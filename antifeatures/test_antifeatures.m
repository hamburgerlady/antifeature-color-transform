nrmax = 200;
im = imread('peppers.png');
im = double(im);
im = im-min(im(:));
im = im/max(im(:));
imgr = repmat(max(im,[],3),[1 1 3]);

[M,N,~] = size(im);
J = im2str_tensor(im);

%% features
ff1 = J2feat(J);
ff2 = J2feat_inhib(J);
ff3 = imsample_antifeat(ff2,M,N,nrmax);

figure(1)
clf
imshow([im im im;imgr imgr imgr]);
hold on

plot(ff1(:,1),ff1(:,2),'.')
plot(ff2(:,1)+N,ff2(:,2),'.')
plot(ff3(:,1)+2*N,ff3(:,2),'.')

%% descriptors
yy1 = antifeat_decriptor(ff1,im);
yy2 = antifeat_decriptor(ff2,im);
yy3 = antifeat_decriptor(ff3,im);



for iii = 1:size(yy1,1)
    ll = plot(ff1(iii,1),ff1(iii,2)+M,'.');
    set(ll,'MarkerSize',20);
    set(ll,'Color',yy1(iii,:));
end

for iii = 1:size(yy2,1)
    ll = plot(ff2(iii,1)+N,ff2(iii,2)+M,'.');
    set(ll,'MarkerSize',20);
    set(ll,'Color',yy2(iii,:));
end

for iii = 1:size(yy3,1)
    ll = plot(ff3(iii,1)+2*N,ff3(iii,2)+M,'.');
    set(ll,'MarkerSize',20);
    set(ll,'Color',yy3(iii,:));
end
