function dd = distapprox_poly_mytest(poly,x,y)

a = poly(1);
b = poly(2);
c = poly(3);

k = 3*a*x.^2+2*b*x+c;
m = (a*x.^3+b*x.^2+c*x)-k*x;
dd = (k*x-y+m).^2./(k*k+1);
