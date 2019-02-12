function [pVal_cor, rValL_cor, rValH_cor] = corr_Bonferroni(rVal,pVal,n,alphaLevel,nComp)
% Bonferroni-corrected p-value and confidence interval for correlational
% coefficient
% Inputs
%   rVal - a correlational coefficient (a single scalar)
%   pVal - a p-value from correlational test (Bonferroni-uncorrected)
%   n - number of observations/subjects in the correlational test
%   alphaLevel - significance alpha level (e.g., 0.05)
%   nComp - number of multiple comparisons
% Outputs
%   pVal_cor - Bonferroni-corrected p-value
%   rValL_cor - lower bound of Bonferroni-corrected confidence interval of
%               the correlational coefficient
%   rValH_cor - upper bound of Bonferroni-corrected confidence interval of
%               the correlational coefficient

% Written by Chi-Hsun Eric Chang, Sep 2018
%%
% convert r to z using fisher-z transformation 
Z_rVal = 0.5*log( (1+rVal)/(1-rVal) );

% Bonferroni-corrected alpha level
alphaLevel_cor = alphaLevel/nComp; 

% Bonferroni-corrected critical z-value
Zcrit = norminv(1-alphaLevel_cor/2);  

% z-value of the corrected confidence interval
Z_rValL = Z_rVal - Zcrit*sqrt(1/(n-3)); % lower bound
Z_rValH = Z_rVal + Zcrit*sqrt(1/(n-3)); % upper bound 

% convert z-value back to r-value
rValL_cor = (exp(2*Z_rValL)-1) / (exp(2*Z_rValL)+1); 
rValH_cor = (exp(2*Z_rValH)-1) / (exp(2*Z_rValH)+1); 

% Bonferroni-corrected p-value
pVal_cor = pVal*nComp; 