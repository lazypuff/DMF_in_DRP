clear;clc;
%%% Part one, read in the splitting information and all other information
Drugsim_fig_mt = readmatrix('identity_gdsc.csv');
Drugsim_fig_mt = Drugsim_fig_mt(:,2:end);
Cellsim_probe = readmatrix('simC_ogdsc_gexp.csv');
Cellsim_probe = Cellsim_probe(:,2:end);

resp_mat = readmatrix('resp_ogdsc.csv');
resp_mat = resp_mat(:, 2:end);
resp_mat = resp_mat';
% second column is the index of cell, third column is the index of drug
fold_info = readmatrix('srmf_valid_info_mask_drug_og.csv');

% create the metric output
RMSE_out = [];
PCC_out =[];
R2_out = [];
results_met = [];
% the output prediction matrix
predMatOut = zeros(size(resp_mat));

% fold iteration
for j = 1:10
    fold_indices = fold_info(:, 1) == j;
    % Extract row indices and column indices
    row_indices = fold_info(fold_indices, 3);
    col_indices = fold_info(fold_indices, 2);



    % Use sub2ind to convert row and column indices to linear indices
    linear_indices = sub2ind(size(resp_mat), row_indices, col_indices);

    % Extract values from resp_mat using linear indices
    scale1 = resp_mat(~isnan(resp_mat));
    num = resp_mat./max(max(scale1),abs(min(scale1)));
    train_df = num;
    valid_df = resp_mat(linear_indices);
    train_df(linear_indices) = [NaN];

    %%% Part three, train the hyper-parameters only use training data
    % 10% testing data in training
    % use drug as row index and cell line as column index
    K = 45; max_iter=50; seed=50;
    index = (1:size(train_df,1));
    index = randsample(index, size(train_df,1));
    results = [];
    for lambda_l = [2^-3, 2^-2, 2^-1, 1, 2, 2^2]
        for lambda_d = [2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 1, 2]
            for lambda_c = [2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 1, 2]
                
                s = round(length(index)/10);
                traintemp = train_df;            
                indextemp = index(1:s);
                    
                testtemp = train_df(indextemp,:);
                traintemp(indextemp,:) = [NaN];
                curnum = traintemp;
                W = ~isnan(curnum);
                curnum(isnan(curnum)) = 0;
                [U,V] = CMF(W,curnum,Drugsim_fig_mt,Cellsim_probe,lambda_l,lambda_d,lambda_c,K,max_iter,seed);
                num = num *max(max(scale1),abs(min(scale1)));
                numpred = U*V'*max(max(scale1),abs(min(scale1)));
                predtemp = numpred(indextemp,:);
                SQE = (testtemp - predtemp).^2;
                RMSE = sqrt(nanmean(SQE(:)));
                
                results = [results;[lambda_l, lambda_d, lambda_c, RMSE]];
            end
        end
    end
    %%% Part Four, get the prediction of test set and calculate the metrics,
    %%% using the hyper-parameters trained in Part Three
    % Extract the hyper-parameters that has the best(lowest) RMSE.

    % Extract the min value of RMSE in training-testing and also the row index
    % of that
    % results = readmatrix('results_grid_search_gdsc_256_raw.csv');
    [minValue, minRowIndex] = min(results(:,4));

    % Extract the values from the first three columns for the row with the lowest fourth column value
    valuesOfInterest = results(minRowIndex, 1:3);
    lambda_l = valuesOfInterest(1);
    lambda_d = valuesOfInterest(2);
    lambda_c = valuesOfInterest(3);

    

    % input the tuned hyperparameters here
    % K = 45; lambda_l = 2^i1; lambda_d = 2^i2; lambda_c = 2^i3; max_iter=50; seed=50;
    K = 45; max_iter=50; seed=50;
    curnum = train_df;
    W = ~isnan(curnum);
    curnum(isnan(curnum)) = 0;
    [U,V] = CMF(W,curnum,Drugsim_fig_mt,Cellsim_probe,lambda_l,lambda_d,lambda_c,K,max_iter,seed);
    num = num *max(max(scale1),abs(min(scale1)));
    numpred = U*V'*max(max(scale1),abs(min(scale1)));
    predtemp = numpred(linear_indices);
    SQE = (valid_df - predtemp).^2;
    RMSE_out = sqrt(nanmean(SQE(:)));
    PCC_out = corr(valid_df, predtemp, 'rows', 'pairwise');
    SST = nansum((valid_df(:) - nanmean(valid_df(:))).^2);
    SSE = nansum((valid_df(:) - predtemp(:)).^2);
    R2_out = (1-SSE/SST);
    results_met = [results_met;[j,lambda_l,lambda_d,lambda_c,RMSE_out, PCC_out, R2_out]];
    
    filterPred = zeros(size(numpred));
    for i = 1:length(row_indices)
    filterPred(row_indices(i), col_indices(i)) = numpred(row_indices(i), col_indices(i));
    end
    
    predMatOut = predMatOut + filterPred;
end

writematrix(results_met,'results_oldgdsc_nd_mask_drug.csv');
writematrix(predMatOut,'predMat_oldgdsc_nd_mask_drug.csv');