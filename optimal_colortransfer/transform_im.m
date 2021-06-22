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

if isfield(T,'P')
    tmp = imt;
    a = imt(:,:,1);
    b = imt(:,:,2);
    c = imt(:,:,3);
    
    for iii = 1:3
        tmp(:,:,iii) = a.*b*T.P(1,iii)+a.*c*T.P(2,iii)+b.*c*T.P(3,iii)+a.*a*T.P(4,iii)+b.*b*T.P(5,iii)+c.*c*T.P(6,iii)+a*T.P(7,iii)+b*T.P(8,iii)+c*T.P(9,iii)+T.P(10,iii);
    end
    imt = tmp;
else
    imt = transform_im_affine(imt,T.A);
end

imt(imt<0)=0;
imt(imt>1)=1;

imt = (imt*T.M+T.m);
