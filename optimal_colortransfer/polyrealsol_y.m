function rr = polyrealsol_y(poly13,d)


lille = 1e-9;

a = poly13(1);
b = poly13(2);
c = poly13(3);

maxr = roots([a b c -1]);
maxr(abs(imag(maxr))>lille)=[];
maxr = real(maxr);
maxr(maxr<0)=[];
[~,id]=min(abs(maxr-1));
maxr = maxr(id);
if isempty(maxr),
    maxr = roots([3*a 2*b c]);
    [~,id]=min(abs(maxr-1));
    maxr = real(maxr(id));
    %keyboard
end

delta = 18*a.*b.*c.*d-4*b.^3.*d+b.^2.*c.^2-4*a.*c.^3-27*a.^2.*d.^2;
delta0 = b.^2-3*a.*c;
delta1 = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
%tmp = sqrt(max(0,-27*a.^2.*delta));
tmp = sqrt(-27*a.^2.*delta);
tmp2 = (delta1+tmp)/2;

aa = angle(tmp2);
nn = (abs(tmp2)).^(1/3);

aa1 = aa/3;
aa2 = aa/3+2*pi/3;
aa3 = aa/3+2*2*pi/3;

C1 = nn.*(cos(aa1)+i*sin(aa1));
C2 = nn.*(cos(aa2)+i*sin(aa2));
C3 = nn.*(cos(aa3)+i*sin(aa3));
%C = ((delta1+tmp)/2).^(1/3);
%C = sign(tmp2).*((abs(tmp2)).^(1/3));
%disp([C.^3; tmp2])
rr1 = -1./a/3.*(b+C1+delta0./C1);
rr2 = -1./a/3.*(b+C2+delta0./C2);
rr3 = -1./a/3.*(b+C3+delta0./C3);


% rr1(abs(imag(rr1))>lille)=0;
% rr2(abs(imag(rr2))>lille)=0;
% rr3(abs(imag(rr3))>lille)=0;

rr1(abs(imag(rr1))>lille)=inf;
rr2(abs(imag(rr2))>lille)=inf;
rr3(abs(imag(rr3))>lille)=inf;

rr1 = real(rr1);
rr2 = real(rr2);
rr3 = real(rr3);

% rr1(rr1<0)=0;
% rr2(rr2<0)=0;
% rr3(rr3<0)=0;
rr1(rr1<-lille)=inf;
rr2(rr2<-lille)=inf;
rr3(rr3<-lille)=inf;

% rr1(rr1>maxr)=0;
% rr2(rr2>maxr)=0;
% rr3(rr3>maxr)=0;
rr = min(rr3,min(rr1,rr2));
%rr = rr1+rr2+rr3;
M = max(rr(isfinite(rr)));
rr(~isfinite(rr))=M;

%rr(d<-1)=1;


