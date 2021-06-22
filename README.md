# antifeature-color-transform
A feature-based method, that robustly fits color transforms to data containing gross outliers. 

To run algorithm, add directories 'optimal_colortransfer' and 'antifeatures' to the matlab path.
Run test script 'test_colortransfer.m' to test the method

Main functions:

[im2t,T] = color_correct_poly(im1,im2); 
(Robustly transforms im2 to the colors of im1. Also returns the transformation details in the struct T.)

imt = transform_im(im,T);
(Transforms im to imt using the transform struct T.)


If you use the code please cite
@inproceedings{oskarsson2021robust,
  title={Robust Image-to-Image Color Transfer Using Optimal Inlier Maximization},
  author={Oskarsson, Magnus},
  booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition Workshops},
  year={2021}
}
