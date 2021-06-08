function imt = transform_im(im,T)

im = double(im);
im = im-min(im(:));
im = im/max(im(:));

if T.rev
    imt = transform_im_full_rev(im,T.poly,T.R);    
else
    imt = transform_im_full(im,T.poly,T.R);
end

imt(imt>1)=1;
imt(imt<0)=0;

imt = transform_im_affine(imt,T.A);
imt(imt<0)=0;
imt(imt>1)=1;


imt = (imt*T.M+T.m);
