function solution
  testPlotFunc(@testFfunc, -2*pi, pi, "cos (3x)");
  testPlotFunc(@testGfunc, -10, 10, "|x|");
endfunction

function retval=testFfunc(x)
  retval = cos(3 * x);
endfunction

function retval=testGfunc(x)
  retval = abs(x);
endfunction

function renderPlots(expTitle, f, a, b, lab, getResults, getNodes)
  RESOLUTION = 1000;
  TEST_RES = 1000;
  
  figure
  rend = linspace (a, b, RESOLUTION);
  plot (rend, f(rend));
  xlabel ("x");
  title (strcat(lab, " - ", expTitle));

  hold on;
  
  for i = [16, 4] # reverse so drawings are clearer
    x = getNodes(a, b, i);
    y = f(x);
    scatter(x, y);
    plot (rend, getResults(x, y, a, b, rend));
  endfor
  
  for i = [4, 16, 64]
    x = getNodes(a, b, i);
    rend = rand(1,TEST_RES);
    y = f(x);
    z = getResults(x, y, a, b, rend);
    v = f(rend);
    disp(cstrcat("norm for ", lab, " - ", expTitle, ", n = ", num2str(i)))
    disp(norm(v - z, Inf))
  endfor
  
  h = legend ("original", "16 nodes", "interpolated from 16 nodes", "4 nodes", "interpolated from 4 nodes");
  
  legend (h, "location", "northeastoutside")

  hold off;
endfunction

function retval=scaleRange(rng, oldA, oldB, newA, newB)
  retval = (((rng - oldA) / (oldB - oldA)) * (newB - newA)) + newA;
endfunction

function retval=getChebyshevNodes(a, b, i)
  normalNodes = cos(((2 * (0:(i-1)) + 1) * pi)/(2*i));
  retval = scaleRange(normalNodes, -1, 1, a, b);
endfunction

function retval=getLagrange(x, y, a, b, rend)
  retval = Horner(interpNewton(x,y), x, rend);
endfunction

function retval=getBSpline(x, y, a, b, rend)
  retval = Bsplval(rend, Bsplnat(y), a, b);
endfunction

function testPlotFunc(f, a, b, lab)
  renderPlots("Lagrange, linspace", f, a, b, lab, @getLagrange, @linspace)
  renderPlots("Lagrange, Chebyshev", f, a, b, lab, @getLagrange, @getChebyshevNodes)
  
  #disabled - broken and slow
  #renderPlots("Bspline, linspace", f, a, b, lab, @getBSpline, @linspace)
endfunction

function c=interpNewton(x,y)
  n = length(x);
  c = y;
  for j = 2:n
    for i = n:-1:j
      c(i) = (c(i) - c(i-1))/(x(i) - x (i-j+1));
    endfor
  endfor
endfunction

function v=Horner(c,x,z)
  n = length(c);
  v = c(n);
  for i = (n-1):-1:1
    v = (v .* (z - x(i))) + c(i);
  endfor
endfunction

function retval=NodeRangeIndex(nodes, x) # TODO binsearch
  for i = 1:(length(nodes)-1)
    if (x < nodes(i+1))
      retval = i;
      return
    endif
    if (1 == nodes(i + 1) && x == 1 && nodes(i) == 1)
      retval = i;
      return
    endif
  endfor
  retval = length(nodes)-1;
endfunction

function retval=buffNodes(nodes, n)
  retval = [zeros(1, n+1), nodes, ones(1, n+5)];
endfunction

function retval=DeBoore(nodes, x, coeffs, n)
  nodes = buffNodes(nodes, n);
  coeffs = buffNodes(coeffs', n);
  k = NodeRangeIndex(nodes, x);
  d = zeros(n + 5, k + 5);
  for i = (k-n):k
    d(1, i+1) = coeffs(i);
  endfor
  
  for j = 1:n
    for i = (k - n + j):k
      alfa = (x - nodes(i+1))/(nodes(i+n+2) - j - nodes(i + 1));
      d(j + 1, i + 1) = ((1 - alfa) * d(j, i)) + (alfa * d(j, i + 1));
    endfor
  endfor
  retval = d(n+1, k + 1);
endfunction

function retval=Mansfield(nodes, x, n)
  nodes = buffNodes(nodes, n);
  f = n + 3 + floor((x - nodes(n+3))/(nodes(n+3)-nodes(n+2))); # assuming equidistance
  k = NodeRangeIndex(nodes, x);
  b = zeros(1, length(nodes) + 3);
  b(k+1) = 1;
  for j = 1:n
    beta = (nodes(k+2) - x)/(nodes(k+2) - nodes(k-j+2));
    b(k - j + 1) = beta * b(k - j + 2);
    for i = (k - j + 1):k
      alfa = 1 - beta;
      beta = (nodes(i+j+2) - x)/(nodes(i+j+2) - nodes(i+2));
      b(i + 1) = alfa * b(i+1) + beta * b(i + 2);
    endfor
    b(k+1) *= (1 - beta);
  endfor
  retval = b;
endfunction

# assuming that nodes are equidistant in range [0,1]
function c = Bsplnat(y)
  DEPTH = 3;
  n = length(y);
  nodes = linspace(0,1,n);
  x = [];
  z = [];
  val = [];
  for i = 3:(n)
    mans = Mansfield(y, nodes(i), DEPTH);
    for j = 1:DEPTH
      x = [x, i + 1 - j];
      z = [z, i - 2];
#      val = [val, mans( 7 + i - i - j)]; # TODO remove hardcoded
      val = [val, mans( DEPTH + 2 + i - i - j)]; # TODO remove hardcoded
    endfor
  endfor
  A = sparse(x,z,val,n,n)';
  c = A \ y';
endfunction

function v=Bsplval(z,c,a,b)
  DEPTH = 3;
  v = []
  n=length(c);
  nodes = linspace(0,1,n);
  z = scaleRange(z, a, b, 0, 1);
  for i = 1:length(z)
    v = [v, DeBoore(nodes, z(i), c, DEPTH)];
  endfor
endfunction
