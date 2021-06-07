function [R,inliers] = opt_rot_fast(xx,yy,ep)

nn = size(xx,2);

x2 = sum(xx.^2);
y2 = sum(yy.^2);
yx = zeros(9,nn);
for iii = 1:nn
    tmp = -2*yy(:,iii)*xx(:,iii)';
    yx(:,iii)=tmp(:);
end
ss = x2+y2-ep;
[R,inliers] = solver_rotmex(xx,yy,ss,yx);

