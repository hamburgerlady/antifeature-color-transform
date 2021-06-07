function [J,lam1,lam2,E1,E2] = im2str_tensor(im)


[M,N,C]=size(im);
if C>1,
    im = double(max(im,[],3));
end
sigma2 = 1.7;

[imx,imy]=str_diff_nosmooth(im);
[J_xx,J_xy,J_yy]=str_J(imx,imy);
[J_xx_sm,J_xy_sm,J_yy_sm]=str_Jsm(J_xx,J_xy,J_yy,sigma2);

j11 = J_xx_sm;
j12 = J_xy_sm;
j22 = J_yy_sm;

J = zeros(M,N,3);
J(:,:,1)=j11;
J(:,:,2)=j12;
J(:,:,3)=j22;

if nargout>1,
lam2 = j11/2 + j22/2 + (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;
lam1 = j11/2 + j22/2 - (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;

E1_1 = (j11/2 + j22/2 - (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2)./j12 - j22./j12;
E1_2 = 1;
E2_1 = (j11/2 + j22/2 + (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2)./j12 - j22./j12;
E2_2 = 1;

sum1 = sqrt(E1_1.^2+E1_2.^2);
sum2 = sqrt(E2_1.^2+E2_2.^2);

E1_1 = E1_1./sum1;
E1_2 = E1_2./sum1;
E2_1 = E2_1./sum2;
E2_2 = E2_2./sum2;

E1 = zeros(M,N,2);
E2 = zeros(M,N,2);
E1(:,:,1)=E1_1;
E1(:,:,2)=E1_2;
E2(:,:,1)=E2_1;
E2(:,:,2)=E2_2;
end
