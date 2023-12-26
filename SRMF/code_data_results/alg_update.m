function [U0]=alg_update(U, V, W, R, S, lambda_l, lambda_d)
X = R*V + 2*lambda_d*S*U;
Y = 2*lambda_d*U'*U;
U0 = zeros(size(U));
nu = size(U,2);
D = V'*V;
[m,n] = size(W);
for i = 1:m
    ii = find(W(i,:)>0);
    if size(ii,2) == 0
       B = Y + lambda_l*eye(nu);
    elseif size(ii,2) == n
       B = D + Y + lambda_l*eye(nu);
    else
       A = V(ii, :)'*V(ii, :);
       B = A + Y + lambda_l*eye(nu);
    end
    U0(i, :) = X(i, :)/B;  % see mldivide
end