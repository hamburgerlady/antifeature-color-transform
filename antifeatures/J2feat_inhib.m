function [ff,ok,logprodeig] = J2feat_inhib(J)

[M,N,~]=size(J);
logprodeig = log(1+J(:,:,1).*J(:,:,3)-J(:,:,2).^2);
facty = 5;
%inhib = 5; % inhib pixels
%inhib = 3; % inhib pixels
inhib = 10; % inhib pixels

mm = mean(logprodeig(:));
ok1 = logprodeig<(mm*facty);
ok2 = false(size(ok1));
lgp_u = logprodeig(1:end-2,2:end-1);
lgp_d = logprodeig(3:end,2:end-1);
lgp_l = logprodeig(2:end-1,1:end-2);
lgp_r = logprodeig(2:end-1,3:end);

ok2(2:end-1,2:end-1) = (logprodeig(2:end-1,2:end-1)<=lgp_u) & ...
                        (logprodeig(2:end-1,2:end-1)<=lgp_d) & ...
                        (logprodeig(2:end-1,2:end-1)<=lgp_l) & ...
                        (logprodeig(2:end-1,2:end-1)<=lgp_r);
ok = ok1 & ok2;
[yi,xi] = find(ok);

nn = length(yi);
id = randperm(nn);
yi = yi(id);
xi = xi(id);
for iii = 1:nn
    yii = yi(iii);
    xii = xi(iii);
    if ok(yii,xii)
        y1 = max(yii-inhib,1);
        x1 = max(xii-inhib,1);
        y2 = min(yii+inhib,M);
        x2 = min(xii+inhib,N);
        ok(y1:y2,x1:x2)=false;
        ok(yii,xii)=true;
        %disp(sum(ok(:)));
    end
end

[yi,xi] = find(ok);

nn = size(yi,1);
ff = zeros(nn,2);
ff(:,1) = xi;
ff(:,2) = yi;


