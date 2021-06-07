function ffout = imsample_antifeat(ffin,M,N,nn)

%bordy = 50;
if size(ffin,1)>nn

bordy = ceil(0.013*max(M,N));

ar = M/N;
nd = sqrt(nn/ar);
Mn = ceil(ar*nd);
Nn = ceil(nd);
%[xx,yy]=meshgrid(linspace(1,N,Nn),linspace(1,M,Mn));
[xx,yy]=meshgrid(linspace(bordy,N-bordy,Nn),linspace(bordy,M-bordy,Mn));
nnin = size(ffin,1);
xx=xx(:);
yy=yy(:);
nne = length(xx);
ffout = zeros(nne,2);
for iii = 1:nne
    [~,indy]=min(sum(abs(ffin-repmat([xx(iii) yy(iii)],nnin,1)),2));
    ffout(iii,:)=ffin(indy,:);
end
pp = randperm(nne);
ffout = ffout(pp(1:nn),:);
else
    ffout = ffin;
end