function [A,inliers] = affine_colreg_RANSAC_mytest(x,y,ep_aff,iters)

warning('off')
if nargin<4
    iters = 1000;
end
if nargin<3
    ep_aff = 1e-3;
end
n = size(x,2);

if n>5

xe = [x;ones(1,n)];
bestins = 0;

for iii = 1:iters
    idi = randperm(n,4);
    xi = xe(:,idi);
    yi = y(:,idi);
    At = yi/xi;
    inst = sum(sum((At*xe-y).^2)<ep_aff);
    if inst>bestins
        bestA = At;
        bestins = inst;
    end
end
%disp(bestins)

inliers = sum((bestA*xe-y).^2)<ep_aff;
A = y(:,inliers)/xe(:,inliers);

else
    inliers = 1:n;
    A = [eye(3) zeros(3,1)];
end


