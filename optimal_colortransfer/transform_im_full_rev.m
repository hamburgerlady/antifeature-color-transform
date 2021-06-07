function imout = transform_im_full_rev(im,poly,R)

if max(im(:))>100,
    im = double(im)/255;
end

imr = sqrt(sum((double(im)).^2,3))/sqrt(3);

%imrt = polyval([poly(1:3)' 0],imr);
imrt = polyrealsol_y(poly(1:3)',-imr);

imc = im./repmat(eps+imr,[1 1 3])/sqrt(3);

r0 = imc(:,:,1);
g0 = imc(:,:,2);
b0 = imc(:,:,3);

% r0t = R(1,1)*(r0-m0(1))+R(1,2)*(g0-m0(2))+R(1,3)*(b0-m0(3))+m1(1);
% g0t = R(2,1)*(r0-m0(1))+R(2,2)*(g0-m0(2))+R(2,3)*(b0-m0(3))+m1(2);
% b0t = R(3,1)*(r0-m0(1))+R(3,2)*(g0-m0(2))+R(3,3)*(b0-m0(3))+m1(3);

r0t = R(1,1)*r0+R(1,2)*g0+R(1,3)*b0;
g0t = R(2,1)*r0+R(2,2)*g0+R(2,3)*b0;
b0t = R(3,1)*r0+R(3,2)*g0+R(3,3)*b0;

imct(:,:,1)=r0t;
imct(:,:,2)=g0t;
imct(:,:,3)=b0t;

imout = imct.*repmat(imrt,[1 1 3])*sqrt(3);
