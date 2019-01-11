function solution
  testPlotFunc(@testFfunc, -2*pi, pi, "cos (3x)");
  testPlotFunc(@testGfunc, -10, 10, "|x|");
  #plotter;
endfunction

function retval=testFfunc(x)
  retval = cos(3 * x);
endfunction

function retval=testGfunc(x)
  retval = abs(x);
endfunction

function renderPlots(expTitle, f, a, b, lab, getResults, getNodes)
  RESOLUTION = 500
  
  figure
  rend = linspace (a, b, RESOLUTION);
  plot (rend, f(rend));
  xlabel ("x");
  #ylabel ("label");
  title (strcat(lab, " - ", expTitle));

  hold on;
  
  for i = [16, 4]
    x = getNodes(a, b, i);
    y = f(x);
    scatter(x, y);
    plot (rend, getResults(x, y, a, b, rend));
  endfor
  h = legend ("original", "16 nodes", "interpolated from 16 nodes", "4 nodes", "interpolated from 4 nodes");
  
  legend (h, "location", "northeastoutside")

  hold off;
endfunction

function retval=getChebyshevNodes(a, b, i)
  normalNodes = cos(((2 * (0:(i-1)) + 1) * pi)/(2*i))
  retval = ((normalNodes + 1) / 2) * (b - a) + a
endfunction

function retval=getLagrange(x, y, a, b, rend)
  retval = Horner(interpNewton(x,y), x, rend)
endfunction

function retval=getBSpline(x, y, a, b, rend)
  retval = Bsplval(rend, Bsplnat(y), a, b)
endfunction

function testPlotFunc(f, a, b, lab)
  renderPlots("Lagrange, linspace", f, a, b, lab, @getLagrange, @linspace)
  renderPlots("Lagrange, Chebyshev", f, a, b, lab, @getLagrange, @getChebyshevNodes)
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

function Bsplnat(y)
endfunction

function Bsplval(z,c,a,b)
endfunction

function plotter
  x = -10:0.1:10;
  plot (x, sin (x));
  xlabel ("x");
  ylabel ("sin (x)");
  title ("Simple 2-D Plot");
endfunction

