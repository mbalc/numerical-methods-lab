function demoLZNK
  SAMPLE_COUNTS=[10, 20, 100];
  TEMPLATE_POLY=[1,-5,2];
  
  printf("Testing %d times with polynomial p:\n", length(SAMPLE_COUNTS))
  disp(TEMPLATE_POLY')
  
  for i=1:length(SAMPLE_COUNTS)
    m=SAMPLE_COUNTS(i);
    printf("\nBeginning test nr %d with m = %d\n", i, m)
    testLZNK(TEMPLATE_POLY, m)
  endfor
endfunction

function testLZNK(polyTemplate, sampleCount)
  [polys,xs,ys]=generateSamples(polyTemplate, sampleCount);
  [A,y]=generateLZNKInput(xs, ys, length(polyTemplate));
  [x,R,B]=Householder(A,y);
  printf("LZNK solution is:\n", sampleCount)
  disp(x)
  restoredQR=rotateMatrix(B, R);
  printf("and ||A - QR|| / ||A|| is %f\n", relativeError(A, restoredQR))
endfunction

function [matrix]=rotateMatrix(householderVectors, matrix)
  if size(householderVectors,2) < 1
    b=[];
    return;
  endif
  matrix=reflectedMatrix(matrix, householderVectors(:,1));
  [matrix(2:end, 2:end)]=rotateMatrix(householderVectors(2:end, 2:end), matrix(2:end, 2:end));
endfunction

function retval=relativeError(original, computed)
  retval=norm(original-computed)/norm(original);
endfunction

function [A, y]=generateLZNKInput(xs, ys, len)
  A = vander(xs, len);
  y = ys;
endfunction

function [polys, xs, ys]=generateSamples(templatePoly, sampleCount)
  xs = rand(sampleCount,1)*10;
  repeatedPoly = repmat(templatePoly, sampleCount, 1);
  polys=[];
  for p=repeatedPoly'
    polys=[polys;perturbPolynomial(p')];
  endfor
  ys = [];
  for i=1:sampleCount
    p=polys(i,:);
    x=xs(i);
    ys=[ys,polyval(p,x)];
  endfor
  ys=ys';
endfunction

function polynom=perturbPolynomial(polynom)
  POLY_PERTURBANCE = 10^-2;
  polynom(end) += (2 * POLY_PERTURBANCE * rand) - POLY_PERTURBANCE;
endfunction

function [x,R,B]=Householder(A, y)
  [R,v,B]=decompose(A,y);
  x = R(1:length(v),:) \ v;
endfunction

function [A,y,b]=decompose(A, y)
  if size(A,2) < 1
    b=[];
    return;
  endif
  v = houseVector(A);
  A = reflectedMatrix(A, v);
  y = reflectedVector(y, v);
  
  [A(2:end, 2:end),y(2:end),b]=decompose(A(2:end, 2:end), y(2:end));
  b = [v';[zeros(1,size(b,2));b]']';
endfunction

function A=reflectedMatrix(A, mirror)
  for i = 1:size(A,2)
    A(:, i) = reflectedVector(A(:, i), mirror);
  endfor
endfunction 

function retval=reflectedVector(vec, mirror)
  lam = 2 / (mirror' * mirror);
  s = mirror' * vec;
  t = lam * s;
  retval=vec - (t * mirror);
endfunction

function retval=houseVector(A)
  firstCol = A(:, 1);
  colNorm = sign(firstCol(1)) * norm(firstCol);
  id = eye(length(firstCol))(:, 1);
  retval=firstCol + (colNorm * id);
endfunction
