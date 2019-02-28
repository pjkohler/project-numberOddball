function axxRCA = axxRCAmake(path_names,W,condsToUse,type,axxRCA)
    % Make axxRCA struct that includes the axx Wave data and the W (Ch x
    % component) matrix from fullRCA
    % INPUTS:
    %   pathnames (required): cell vector of string directory names housing Axx_c00x.mat
    %   condsToUse: vector of conditions to use
    % 
    % OUTPUTS:
    %   Struct with the following fields
    %   Wave: Wave data from Axx_c00x.mat
    %   W: linear transformation matrix to go from sensor space to RC-space

    if nargin<6
        axxRCA.Wave = {};
    end
    for f = 1:length(path_names)
        axx_matfiles = arrayfun(@(x) ... 
            sprintf('%s/%s_Exp_MATL_HCN_128_Avg/Axx_c%03d.mat',fileparts(path_names{f}),type,x),condsToUse,'uni',false);

        tmpStrct = cellfun(@(x) load(x),axx_matfiles,'uni',false);
        readyStrct = cellfun(@(x) mrC.axx(x),tmpStrct,'uni',false);
        axxRCA.Wave(:,f) = cellfun(@(x) x.Wave, readyStrct,'uni',false);
    end

    nanDims = [1,2];
    axxRCA.Wave = cellfun(@(x) Zero2NaN(x,nanDims),axxRCA.Wave,'uni',false);
    axxRCA.Projected = rcaProject(axxRCA.Wave,W);
    axxRCA.Projected = cellfun(@(x) Zero2NaN(x,nanDims),axxRCA.Projected,'uni',false);
end


