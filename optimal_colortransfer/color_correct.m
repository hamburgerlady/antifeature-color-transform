function [im2t,poly21,R21,A21]=color_correct(im1,im2)


im1 = double(im1);
im2 = double(im2);

m1 = min(im1(:));
im1 = im1-m1;
M1 = max(im1(:));
im1 = im1/M1;
im2 = im2-min(im2(:));
im2 = im2/max(im2(:));

[M,N,~]=size(im1);

ep_poly = 5e-4;
ep_rot = 1e-2;
ep_aff = 2e-2;

maxins = 50;


%% process images

J = im2str_tensor(im1);
ff0 = J2feat_inhib(J);

yy0 = antifeat_decriptor(ff0,im1);
%xx0 = antifeat_decriptor(ff0,im2);

if length(ff0)>maxins
    ff = imsample_antifeat(ff0,M,N,maxins);
else
    ff = ff0;
end

yy = antifeat_decriptor(ff,im1);
xx = antifeat_decriptor(ff,im2);

rx = sqrt(sum(xx.^2,2))/sqrt(3);
ry = sqrt(sum(yy.^2,2))/sqrt(3);


im1g = max(im1,[],3);
meanie1 = mean(im1g(:));

im2g = max(im2,[],3);
meanie2 = mean(im2g(:));


%%


if meanie2<meanie1
    rev = 1;
else
    rev = 0;
end


if rev
    [~,inliers_poly21] = solver_polymex(ry,rx,ep_poly);
    poly21 = fitinliers_poly(ry(inliers_poly21),rx(inliers_poly21));
    
else
    [~,inliers_poly21] = solver_polymex(rx,ry,ep_poly);
    poly21 = fitinliers_poly(rx(inliers_poly21),ry(inliers_poly21));
end

cx = xx./repmat(rx,[1 3])/sqrt(3);
cy = yy./repmat(ry,[1 3])/sqrt(3);

cx = cx(inliers_poly21,:);
cy = cy(inliers_poly21,:);


R21 = opt_rot_fast(cx',cy',ep_rot);


if rev
    im2t = transform_im_full_rev(im2,poly21,R21);    
else
    im2t = transform_im_full(im2,poly21,R21);
end

im2t(im2t>1)=1;
im2t(im2t<0)=0;

yy21 = yy0;
xx21 = antifeat_decriptor(ff0,im2t);
err_t = sum((xx21-yy21).^2,2);
yy21 = yy21(err_t<ep_aff,:);
xx21 = xx21(err_t<ep_aff,:);


A21 = affine_colreg(xx21',yy21',inf);

im2ta = transform_im_affine(im2t,A21);
im2ta(im2ta<0)=0;
im2ta(im2ta>1)=1;


im2t = (im2ta*M1+m1);



