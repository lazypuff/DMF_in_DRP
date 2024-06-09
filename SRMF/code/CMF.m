function [bestU,bestV]=CMF(W, intMat, drugMat, cellMat, lambda_l,lambda_d,lambda_c,K,max_iter,seed)
m = size(W,1);
n = size(W,2);
rng('default');
rng(seed);
U0 = sqrt(1/K)*randn(m, K);
V0 = sqrt(1/K)*randn(n, K);
bestU = U0;
bestV = V0;
last_loss = compute_loss(U0, V0, W, lambda_l,lambda_d,lambda_c, intMat, drugMat, cellMat);
bestloss = last_loss;
WR = W .* intMat;
% fid = fopen(['record_K' int2str(K) '_lambda_l=' num2str(lambda_l) '_lambda_d=' num2str(lambda_d) '_lambda_c=' num2str(lambda_c) '.txt'],'wt+');
for t = 1:max_iter
    U = alg_update(U0, V0, W, WR, drugMat, lambda_l, lambda_d);
    V = alg_update(V0, U, W', WR', cellMat, lambda_l, lambda_c);
    curr_loss = compute_loss(U, V, W, lambda_l,lambda_d,lambda_c, intMat, drugMat, cellMat);
%     if curr_loss > last_loss
%         U = U0;
%         V = V0;
% %         disp('Convergence failed.');
%         break
%     end
    if curr_loss < bestloss
        bestU = U;
        bestV = V;
        bestloss = curr_loss;
    end
    delta_loss = (curr_loss-last_loss)/last_loss;
%     fprintf(fid,'%s\n',[sprintf('Iter = \t'),int2str(t),...
%         sprintf('\t curr_loss = \t'), num2str(curr_loss),...
%         sprintf('\t delta_loss = \t'), num2str(delta_loss)]);
    if abs(delta_loss) < 10^(-6)
       break
    end
    last_loss = curr_loss;
    U0 = U;
    V0 = V;
end
% fclose(fid);        
    
    

