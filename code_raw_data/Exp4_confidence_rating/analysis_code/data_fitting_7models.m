
% fit data to different models

% model 1: vision dominance
% model 2: motion dominance
% model 3: bayesian integration
% model 4: bayesian alternation
% model 6: proximity-based integration
% model 7: proximity-based alternation

whichModel=input('which model?');

cd C:\Users\Green\Documents\Documents_thinkpad_05_28_2013\Documents\experiments\CueCombination\2environments_doubledtrials_data_flag

quar=input('boxplot outlier interquartile');
 load(['rich_poor_environments_formodeling' num2str(quar) '.mat']);

%% for rich environment

Nsub=18;
% response {i,n,k}  i, subject ID; n, rich or poor; k, condition
for i=1:Nsub
    for n=1:2
        for k=1:4
            response_dist{i,n,k}=sqrt( response{i,n,k}(:,1) .^2 + response{i,n,k}(:,2) .^2);
        end
    end
end


LH3=0; % likelihood for combination
LH4=0; % likelihood for ponflict
for n=1:2
for i=1:Nsub,
    
    disparity1(i,n)=sqrt( (centers{i,n,1}(1)- 0)^2 + (centers{i,n,1}(2)- 0)^2 );
    disparity2(i,n)=sqrt( (centers{i,n,2}(1)- 0)^2 + (centers{i,n,2}(2)- 0)^2 );
    mu2=0;
    mu1=disparity(i,n);
    std1=std(response{i,n,1});
    std2=std(response{i,n,2});
    
    if whichModel==1  % vision dominance
        std3=std1; % combination condition
        mu3=mu1;
        
        
        
        
        std4=std1; % conflict condition
        mu4=mu1+mean{targ_lm};
        
        
        
        
   