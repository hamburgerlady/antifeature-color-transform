function A = affine_colreg(xx,yy,ep)

gids = sum((xx-yy).^2)<=ep;

if sum(gids)>5
xx = [xx(:,gids);ones(1,sum(gids))];
yy = yy(:,gids);

A = yy/xx;
else
    A = [eye(3) zeros(3,1)];
end
