function [meanDiff_act,meanDiff_perm,p_val_Perm,CI_perm] = indepT_perm(data1,data2,nPerm,alphaLevel,tail)
% Permutation test for two-group comparison
% Inputs:
%       data1 - data from group 1
%       data2 - data from group 2
%       nPerm - number of permutation times
%       alpha level - significance level (alpha)
%       tail - one-tailed (1 or -1) or two-tailed p-value (2)
% Outputs:
%       meanDiff_act - actual mean difference between two groups
%       meanDiff_perm - permuted mean difference. A 1 x nPerm vector
%       p_val - p value of the permutation test
%       CI_perm - confidence interval of the permutation test

% Writtien by Chi-Hsun Eric Chang, April, 2018

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

% Check missing data and if so, remove it
% for group1
ind_missData1 = find(isnan(data1));
nMissData1 = sum(isnan(data1));
if nMissData1==0
    n1 = size(data1,1);
else
    n1 = size(data1,1) - nMissData1;
    data1(ind_missData1) = [];
    disp('Missing data in data1')
end
% for group2
ind_missData2 = find(isnan(data2));
nMissData2 = sum(isnan(data2));
if nMissData2==0
    n2 = size(data2,1);
else
    n2 = size(data2,1) - nMissData2;
    data2(ind_missData2) = [];
    disp('Missing data in data2')
end

% Combine two groups
dataAll = [data1; data2];
groupLabel = [repmat(1,n1,1); repmat(2,n2,1)];

% Actual difference between group1 and group2
meanDiff_act = mean(data1) - mean(data2);

% Permutation
for i = 1:nPerm
    groupLabel_rand = groupLabel(randperm(n1+n2));
    ind_group1 = find(groupLabel_rand==1);
    ind_group2 = find(groupLabel_rand==2);
    group1 = dataAll(ind_group1,1);
    group2 = dataAll(ind_group2,1);
    
    meanDiff_perm(i) = mean(group1) - mean(group2); 
end

% Compute the p value
% default two-tailed test
tmp_sort = sort(abs(meanDiff_perm), 'descend');
tmp_act = abs(meanDiff_act);
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

    
    


