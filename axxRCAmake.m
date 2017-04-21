function axxRCA = axxRCAmake(folderNames,W,condsToUse,c,type,axxRCA)
%Make axxRCA struct that includes the axx Wave data and the W (Ch x
%component) matrix from fullRCA
%INPUTS:
% pathnames (required): cell vector of string directory names housing Axx_c00x.mat
% condsToUse: vector of conditions to use
% 
% OUTPUTS:
% Struct with the following fields
% Wave: Wave data from Axx_c00x.mat
% W: linear transformation matrix to go from sensor space to RC-space

if nargin<6
    axxRCA.Wave = {};
end
%axxRCA(c).W = W;
for f = 1:length(folderNames)
    tempFolders = subfolders(folderNames{f},1);
    %axxpathNames{f} = sprintf('%s/Exp_MATL_HCN_128_Avg',tempFolders{end});
    %axxMatFiles = subfiles(sprintf('%s/Exp_MATL_HCN_128_Avg/Axx*_trials.mat',tempFolders{end}),1);
    axxMatFiles = arrayfun(@(x) ... 
        sprintf('%s/%s_Exp_MATL_HCN_128_Avg/Axx_c%03d.mat',tempFolders{end},type,x),condsToUse,'uni',false);
    tmpStrct = {load(axxMatFiles{(1)}) load(axxMatFiles{(2)})};
    %tmpStrct = load(axxMatFiles{m});
    %readyStrct = mrC.axx(tmpStrct);
    readyStrct = {mrC.axx(tmpStrct{1}),mrC.axx(tmpStrct{2})};
    axxRCA(c).Wave(:,f) = {readyStrct{1}.Wave, readyStrct{2}.Wave};
end

nanDims = [1,2];
axxRCA(c).Wave = cellfun(@(x) Zero2NaN(x,nanDims),axxRCA(c).Wave,'uni',false);
%axxRCA(c).Projected = rcaProject(axxRCA(c).Wave,axxRCA(c).W);
axxRCA(c).Projected = rcaProject(axxRCA(c).Wave,W);
axxRCA(c).Projected = cellfun(@(x) Zero2NaN(x,nanDims),axxRCA(c).Projected,'uni',false);



