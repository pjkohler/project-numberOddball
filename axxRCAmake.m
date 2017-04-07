function axxRCA = axxRCAmake(folderNames,W,condsToUse,c,axxRCA)
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

if nargin<5
    axxRCA.Wave = {};
end
axxRCA(c).W = W;
for f = 1:length(folderNames)
    tempFolders = subfolders(folderNames{f},1);
    %axxpathNames{f} = sprintf('%s/Exp_MATL_HCN_128_Avg',tempFolders{end});
    %axxMatFiles = subfiles(sprintf('%s/Exp_MATL_HCN_128_Avg/Axx*_trials.mat',tempFolders{end}),1);
    axxMatFiles = arrayfun(@(x) ... 
        sprintf('%s/Exp_MATL_HCN_128_Avg/Axx_c%03d.mat',tempFolders{end},x),condsToUse,'uni',false);
    tmpStrct = {load(axxMatFiles{(1)}) load(axxMatFiles{(2)})};
    %tmpStrct = load(axxMatFiles{m});
    %readyStrct = mrC.axx(tmpStrct);
    readyStrct = {mrC.axx(tmpStrct{1}),mrC.axx(tmpStrct{2})};
    axxRCA(c).Wave(:,f) = {readyStrct{1}.Wave, readyStrct{2}.Wave};
end




% %Create AxxRCA struct for time domain analysis of RCA components
%         tmpStrct = {load(axxMatFiles{curCond(1)}) load(axxMatFiles{curCond(2)})};
%         readyStrct = {mrC.axx(tmpStrct{1}),mrC.axx(tmpStrct{2})};
%         axxRCA(c).Wave(:,f) = {readyStrct{1}.Wave, readyStrct{2}.Wave};
%         axxRCA(c).W(:,f) = fullRCA(c).W(curCond,f); %Check if this line runs when data finishes exporting
%         
%         
% for c = 1:(length(condsToUse)/2) %freq pairs
%         curCond = [c-1,c]+c;
%         tmpStrct = {load(axxMatFiles{curCond(1)}) load(axxMatFiles{curCond(2)})};
%         %tmpStrct = load(axxMatFiles{m});
%         %readyStrct = mrC.axx(tmpStrct);
%         readyStrct = {mrC.axx(tmpStrct{1}),mrC.axx(tmpStrct{2})};
%         axxRCA(c).Wave(:,f) = {readyStrct{1}.Wave, readyStrct{2}.Wave};
%     end
