%% load two sample images
im1 = imread(['testims' filesep 'lily1.png']);
im2 = imread(['testims' filesep 'lily2.png']);

%% Test algorithm without outliers

% run with polynomial postprocessing, robust and very accurate
im2t = color_correct_poly(im1,im2);

figure(1)
imshow([im2 im1 im2t])
title('Left: input, Middle: target, Right: output (polynomial)')

disp('Without outliers:');
disp(['PSNR (polynomial): '  num2str(psnr(im1,uint8(im2t)))])
disp(['SSIM (polynomial): '  num2str(ssim(im1,uint8(im2t)))])

%% Test algorithm with outliers


im1c = im1;
for iii = 1:20
    i0 = randi(500);
    j0 = randi(300);
    im1c(i0+(1:50),j0+(1:50),:) = uint8(repmat(rand(1,1,3)*255,[50 50 1]));
end

tmp = abs(double(im1c(:))-double(im1(:)));

out_ratio = sum(tmp>5)/numel(tmp);

disp(['With ' num2str(round(100*out_ratio)) '% outliers:']);


% run with polynomial postprocessing, robust and very accurate
im2t = color_correct_poly(im1c,im2);

figure(2)
imshow([im2 im1c im2t])
title('Left: input, Middle: target, Right: output (polynomial)')

disp(['PSNR (polynomial): '  num2str(psnr(im1,uint8(im2t)))])
disp(['SSIM (polynomial): '  num2str(ssim(im1,uint8(im2t)))])
