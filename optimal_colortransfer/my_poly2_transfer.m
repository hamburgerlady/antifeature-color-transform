function [IR,T] = my_poly2_transfer(X0,X1,I0,ep)


gids = sum((X0-X1).^2,2)<=ep;
nn = sum(gids);

if nn>20
    X0 = X0(gids,:);
    X1 = X1(gids,:);
    a = X0(:,1);
    b = X0(:,2);
    c = X0(:,3);
    
    X0 = [a.*b a.*c b.*c a.*a b.*b c.*c a b c ones(nn,1)];
    T = X0\X1;
    IR = I0;
    a = I0(:,:,1);
    b = I0(:,:,2);
    c = I0(:,:,3);
    
    for iii = 1:3
        IR(:,:,iii) = a.*b*T(1,iii)+a.*c*T(2,iii)+b.*c*T(3,iii)+a.*a*T(4,iii)+b.*b*T(5,iii)+c.*c*T(6,iii)+a*T(7,iii)+b*T(8,iii)+c*T(9,iii)+T(10,iii);
    end
    IR(IR<0)=0;
    IR(IR>1)=1;
else
    IR = I0;
end

