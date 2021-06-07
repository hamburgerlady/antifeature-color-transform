function dd = antifeat_decriptor(ff,im)


nn = size(ff,1);
[M,N,K] = size(im);
ff = max(1,ceil(ff));
ff(:,1) = min(N,ff(:,1));
ff(:,2) = min(M,ff(:,2));
indy = sub2ind([M,N],ff(:,2),ff(:,1));
dd = zeros(nn,K);
%sig = 2;
sig = 1.5;

[xx,yy]=meshgrid(-5:5);
gk = exp(-(xx.^2+yy.^2)/2/sig^2);
gk = gk/sum(gk(:));
for iii = 1:K,
    im_sm = conv2(double(im(:,:,iii)),gk,'same');
    dd(:,iii) = im_sm(indy);
end
