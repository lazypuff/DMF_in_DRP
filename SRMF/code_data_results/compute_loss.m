function [loss]=compute_loss(U, V, W, lambda_l,lambda_d,lambda_c, intMat, drugMat, cellMat)
loss = sum(sum((W .* (intMat - U*V')).^2));
loss = loss + lambda_l * (sum(sum(U.^2)) + sum(sum(V.^2)));
loss = loss + lambda_d * sum(sum((drugMat-U*U').^2)) + lambda_c * sum(sum((cellMat-V*V').^2));
