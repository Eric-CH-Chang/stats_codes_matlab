function [tVal, df, pVal] = CrawfordHowell(singleCase,controlGroup)
% Crawford-Howell t-test for single case studies
% Inputs
%   singleCase - data from the single case subject (e.g., a patient)
%   controlGroup - data from the control group (a column vector)
% Outputs
%   tVal - t-value of the 2-tailed test
%   df - degrees of freedom
%   pVal - p-value of the 2-tailed test

% Writtien by Chi-Hsun Eric Chang, Oct 2017

%%
nContSubs = length(controlGroup);
tVal = (singleCase - mean(controlGroup)) ./ (std(controlGroup).*sqrt((nContSubs+1) ./ nContSubs));
df = nContSubs-1;
pVal = 2*(1-tcdf(abs(tVal),df));
