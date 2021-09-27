function [im2t,T]=color_correct_out_mytest(im1,im2,outs)


im1 = double(im1);
im2 = double(im2);

m1 = min(im1(:));
im1 = im1-m1;
M1 = max(im1(:));
im1 = im1/M1;
im2 = im2-min(im2(:));
im2 = im2/max(im2(:));

[M,N,~]=size(im1);

%ep_poly = 5e-4;
ep_poly = 1e-3;

%ep_rot = 1e-2;
%ep_aff = 2e-2;
%ep_aff = 1e-2;
ep_aff = 5e-3;


maxins = 50;


%% process images

J = im2str_tensor(im1);
ff0 = J2feat_inhib(J);

yy0 = antifeat_decriptor(ff0,im1);
ny = size(yy0,1);
outids = randperm(ny);
ony = ceil(outs*ny);
outids = outids(1:ony);
yy0(outids,:)=rand(ony,3);



xx0 = antifeat_decriptor(ff0,im2);

if length(ff0)>maxins
    [ff,ids_f] = imsample_antifeat(ff0,M,N,maxins);
else
    ff = ff0;
    ids_f = 1:length(ff0);
end


%yy = antifeat_decriptor(ff,im1);
%xx = antifeat_decriptor(ff,im2);

yy = yy0(ids_f,:);
xx = xx0(ids_f,:);

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

% cx = xx./repmat(rx,[1 3])/sqrt(3);
% cy = yy./repmat(ry,[1 3])/sqrt(3);
% 
% cx = cx(inliers_poly21,:);
% cy = cy(inliers_poly21,:);
% 
% 
% R21 = opt_rot_fast(cx',cy',ep_rot);
R21 = eye(3);

if rev
    im2t = transform_im_full_rev(im2,poly21,R21);    
else
    im2t = transform_im_full(im2,poly21,R21);
end

im2t(im2t>1)=1;
im2t(im2t<0)=0;

yy21 = yy0;
xx21 = antifeat_decriptor(ff0,im2t);





if 0
err_t = sum((xx21-yy21).^2,2);
yy21 = yy21(err_t<ep_aff,:);
xx21 = xx21(err_t<ep_aff,:);
    
    
A21 = affine_colreg(xx21',yy21',inf);

im2ta = transform_im_affine(im2t,A21);
im2ta(im2ta<0)=0;
im2ta(im2ta>1)=1;


im2t = (im2ta*M1+m1);
else
    %A21 = [eye(3) zeros(3,1)];
    rx = sqrt(sum(xx21.^2,2))/sqrt(3);
    ry = sqrt(sum(yy21.^2,2))/sqrt(3);
    in21 = ((rx-ry).^2/2)<ep_poly;
    [A21,inliers] = affine_colreg_RANSAC_mytest(xx21(in21,:)',yy21(in21,:)',ep_aff);
    xpol = xx21(in21,:);
    xpol = xpol(inliers,:);
    ypol = yy21(in21,:);
    ypol = ypol(inliers,:);
    
    im2tp = my_poly2_transfer(xpol,ypol,im2t,inf);
    %A21 = yy21(in21,:)'/[xx21(in21,:)';ones(1,sum(in21))];
    
    %im2ta = transform_im_affine(im2t,A21);
    %im2ta(im2ta<0)=0;
    %im2ta(im2ta>1)=1;
    %im2t = (im2ta*M1+m1);
    im2tp(im2tp<0)=0;
    im2tp(im2tp>1)=1;
    im2t = (im2tp*M1+m1);  
    
end

   


T.poly = poly21;
T.R = R21;
T.A = A21;
T.M = M1;
T.m = m1;
T.rev = rev;


