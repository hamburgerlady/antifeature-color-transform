function [ff,ok,logprodeig] = J2feat(J)

logprodeig = log(1+J(:,:,1).*J(:,:,3)-J(:,:,2).^2);
%facty = 5;
facty = 2;

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
nn = size(yi,1);
ff = zeros(nn,2);
ff(:,1) = xi;
ff(:,2) = yi;


