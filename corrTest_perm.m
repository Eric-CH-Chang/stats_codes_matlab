function [rVal_act,rVal_perm,p_val_Perm,CI_perm] = corrTest_perm(data1,data2,nPerm,alphaLevel,tail)
% Permutation test for correlational coefficient
% Inputs:
%       data1 - data from group 1. A column vector
%       data2 - data from group 2. A column vector
%       nPerm - number of permutation times
%       alpha level - significance level (alpha)
%       tail - one-tailed or two-tailed p-value
% Outputs:
%       corr_act - actual correlational coefficient between two groups
%       corr_perm - permuted correlational coefficients. A 1 x nPerm vector
%       p_val - p value of the permutation test
%       CI_perm - confidence interval of the permutation test

% Writtien by Chi-Hsun Eric Chang, April 2018

%%
% Check the number of input and set the default
if nargin < 3
    nPerm = 1000;
end

if nargin < 4
    alphaLevel = 0.05;
end

if nargin < 5
    tail = 2;
end

% Check nPerm and alphaLevel
if alphaLevel<=0.05 && nPerm<1000
    disp('Warning: ')
    disp('Too few permutations for alpha level equal or smaller than .05')
elseif alphaLevel<=0.01 && nPerm < 10000
    disp('Warning: ')
    disp('Too few permutations for alpha level equal or smaller than .01')
end

% Check if the number of data is equal in two groups
n1 = length(data1);
n2 = length(data2);
if ~isequal(n1,n2)
    error('Number of data is not equal between two groups')
end

% Check if the inputs of data are column vectors
if isrow(data1)
    data1 = data1';
elseif ~isrow(data1) && size(data1,2)>1
    error('Data1 must be a column vector')
end

if isrow(data2)
    data2 = data2';
elseif ~isrow(data2) && size(data2,2)>1
    error('Data2 must be a column vector')
end

% Check missing data and if so, remove it
% for group1 and group2
ind_missData1 = find(isnan(data1));
ind_missData2 = find(isnan(data2));
ind_missData = unique([ind_missData1 ind_missData2]);
nMissData = length(ind_missData);
if nMissData==0
    n1 = size(data1,1);
    n2 = size(data2,1);
else
    n1 = size(data1,1) - nMissData;
    data1(ind_missData) = [];
    disp('Missing data in data1')
    
    n2 = size(data2,1) - nMissData;
    data2(ind_missData) = [];
    disp('Missing data in data2')
end

% Actual correlational coefficient between group1 and group2
[rVal, pVal] = corrcoef(data1,data2);
rVal_act = rVal(1,2);
% pVal_act = pVal(1,2);

% Permutation
% keep data from group 1 constant, and randomly shuffl data from group 2 in
% each iteration
for i = 1:nPerm
    ind_rand = randperm(n2);
    data2_rand = data2(ind_rand);
    [rVal, pVal] = corrcoef(data1,data2_rand);
    rVal_perm(i) = rVal(1,2);
%     pVal_perm(i) = pVal(1,2);
end

% Compute the p value
% default two-tailed test
tmp_sort = sort(abs(rVal_perm), 'descend');
tmp_act = abs(rVal_act);
rnk = find(tmp_sort>=tmp_act, 1, 'last'); % how many permuted accuracy is greater or equal to the actual accuracy
if isempty(rnk)
    rnk = 0;
end
p_val_Perm = rnk/nPerm; % the proportion of the values that are equal or greater than actual value
if tail==1 || tail==-1
    p_val_Perm = p_val_Perm/2; % 1-tailed;
end

% Compute confidence interval of the permuted values
CI_perm(1) = prctile(tmp_sort,100*alphaLevel/2); %CI lower bound
CI_perm(2) = prctile(tmp_sort,100-100*alphaLevel/2); % CI upper bound


end
