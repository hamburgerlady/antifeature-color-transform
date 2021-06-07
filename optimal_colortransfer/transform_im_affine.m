function imout = transform_im_affine(im,A)

imout = im;
imout(:,:,1) = A(1,1)*im(:,:,1)+A(1,2)*im(:,:,2)+A(1,3)*im(:,:,3)+A(1,4);
imout(:,:,2) = A(2,1)*im(:,:,1)+A(2,2)*im(:,:,2)+A(2,3)*im(:,:,3)+A(2,4);
imout(:,:,3) = A(3,1)*im(:,:,1)+A(3,2)*im(:,:,2)+A(3,3)*im(:,:,3)+A(3,4);
