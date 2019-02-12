function output = descripStats(dataMat,rowNames,dispTable,tableTitle)
% Conduct descriptive statistics and report them in a table
% Inputs:
%       dataMat - data matrix. Rows are participants and columns are
%                 conditions
%       rowNames - names of conditions used in the table. Need to be cell
%                  vector with strings
%       dispTable - whether displaying the table (1 or 0)
%       tableTitle - the title of the table (optional)
% Output:
%       output - a struct contains descriptive statitics 

% Written by Chi-Hsun Eric Chang, November, 2018

%%
if nargin < 4
    tableTitle = '';
end

if ~strcmp(class(dataMat),'double') || ~strcmp(class(dataMat),'single')
    dataMat = double(dataMat);
end

%%
nSubs = size(dataMat,1); % number of subjects
nConds = size(dataMat,2); % number of conditions
% nGroups = size(dataMat,3); % number of groups
output.nSubs = nSubs;
output.nConds = nConds;

% Mean
averSub = mean(dataMat,1); % average across subjects
output.aver = averSub;

% Standard deviation
sdSub = std(dataMat,0,1);
output.sd = sdSub;

% Standard error
seSub = sdSub/sqrt(nSubs);
output.se = seSub;

% Range
minSub = min(dataMat,[],1);
maxSub = max(dataMat,[],1);
output.minVal = minSub;
output.maxVal = maxSub;

% Table
if dispTable
    disp(tableTitle)
    Headers = {'N','Min','Max','Mean','STD','SE'};
    dscrpStatTable = table(repmat(nSubs,[nConds,1]), minSub', maxSub',...
        averSub', sdSub', seSub',...
        'VariableNames',Headers,'RowNames',rowNames);
    disp(dscrpStatTable)
    fprintf('\n')
end

end