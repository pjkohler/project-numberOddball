function numberOddball_Exp1(varargin)
    % Description:	analyze data from numberOddball Experiment 2
    % 
    % Syntax:	numberOddball_Exp2(<options>)
    % <options>
    %   chanToCompare   - integer indicating the comparison channel [75]
    %    
    %   foldersToUse    - vector of indices of subject folders to analyze
    %                       [analyze all available subject folders]
    %
    %   launchRCA       - if true RCA is run, if false RCA data is loaded
    %                   	[true]/false
    %
    %   trialError      - compute trial-wise error true/[false]
    %
    %   plotSplit       - plot RC run separately on odd and carrier
    %                       true/[false]
    %
    %   forceSourceData - force reading in of raw data when initial running
    %                       RCA [true]/false
    %
    %   axxType         - string indicating the type of axx data to plot
    %                      ['ALL']/'NF1'
    %   plotSNR         - if true use SNR instead of vector based stats
    %                     true/[false]
    
    %% ADD PATHS
    close all;
    if ~exist('code_folder','var')
        code_folder = '/Users/kohler/code';
        rcaCodePath = sprintf('%s/git/rcaBase',code_folder);
        addpath(genpath(rcaCodePath));
        addpath(genpath(sprintf('%s/git/mrC',code_folder)));
        addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',code_folder)));
        addpath(genpath(sprintf('%s/git/export_fig',code_folder)));
    else
    end
    setenv('DYLD_LIBRARY_PATH','')
    
    %% PARSE ARGS
    opt	= ParseArgs(varargin,...
            'chanToCompare'		, 75, ...
            'foldersToUse', [], ...
            'launchRCA', false, ...
            'trialError', false, ...
            'plotSplit', false, ...
            'forceSourceData', false, ...
            'axxType', 'ALL', ...
            'plotSNR', false, ...
            'ampTest', 0, ...
            'behSplit', false ...
            );

    %% VARIABLES THAT CAN BE CHANGED
    topFolder = '/Volumes/Denali_DATA1/kohler/EEG_EXP/DATA/numeroOddball';
    doExp = 1;
    cond_names = {'dist 3','control'};

    %% IDENTIFY DATA LOCATION
    data_location = sprintf('%s/Experiment%.0d',topFolder,doExp);
    fig_location = sprintf('%s/Experiment%.0d/figures',topFolder,doExp);
    if ~exist(fig_location,'dir')
        mkdir(fig_location);
    else
    end

    folder_names=subfolders(sprintf('%s/*20*',data_location),1);
    if ~isempty(opt.foldersToUse)
        folder_names = folder_names(opt.foldersToUse);
    else
    end
    numSubs = size(folder_names,1);
    % and string for saving the data
    saveStr = datestr(clock,26);
    saveStr(strfind(saveStr,'/')) ='';
    saveFileName = sprintf('%s/rcaData_%s.mat',data_location, saveStr); % include the date as a string;

    %% SETUP INPUTS
    use_bins = 0; % use average bin, the zeroeth one
    use_freqs = 1:8; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
    carrier_freq = [6.0, 3.75, 3] ;
    odd_freq = [1.0, 0.75, 0.5 ];
    axx_reps = carrier_freq./odd_freq;
    
    condsToUse = 1:6;
    trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    n_reg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    n_comp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    n_freq = length(use_freqs);
    n_cond = length(condsToUse);
    data_type = 'RLS'; % can also be 'DFT' if you have DFT exports
    rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

    %% call the function --- USER runs this section (no editing necessary) ---
    % generate "filter" for comparison data
    wComparison=zeros(128,1); wComparison(opt.chanToCompare)=1; 

    if ~opt.launchRCA % if users says analysis is not to be launched
        tempFile = dir([data_location,'/rcaData_*.mat']);
        saveFileName = fullfile(data_location,tempFile(end).name(1:(end-4))); % grab newest file and loose .mat suffix;
        load(saveFileName);
    else
        for f = 1:length(folder_names)
            temp_folders = subfolders(folder_names{f},1);
            temp_folders = temp_folders(~ismember(temp_folders, [folder_names{f},'/not_time_corrected']));
            temp_folders = temp_folders(~ismember(temp_folders, [folder_names{f},'/time_corrected']));
            path_names{f} = sprintf('%s/Exp_TEXT_HCN_128_Avg',temp_folders{end});
        end
        if ~exist([saveFileName,'.mat'], 'file')
            save(saveFileName,'*ToUse') % create file if it does not exist
        else
            save(saveFileName,'*ToUse','-append') % otherwise append
        end
        warning('off','all')
        for c = 1:(length(condsToUse)/2)
            curCond = [c-1,c]+c;
            % full RCA, use all harmonics 1F1, 2F1, 1F2 and 2F2
            if c==1
                rca_full(c,:) = rcaSweep(path_names,use_bins,use_freqs,curCond,trialsToUse,n_reg,n_comp,data_type,opt.chanToCompare,[],rcPlotStyle,opt.forceSourceData); % TRUE
            else
                % since this is just a subset of the previous RCA, set forceSourceData to false
                rca_full(c,:) = rcaSweep(path_names,use_bins,use_freqs,curCond,trialsToUse,n_reg,n_comp,data_type,opt.chanToCompare,[],rcPlotStyle,false);
            end
            if c == 3
                for f = 1:size(rca_full,2)
                    rca_full(c,f) = flipSwapRCA(rca_full(c,f),1:5,1:2);
                end
            else
            end
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/rca_full_cond%.0d&%.0d_cov.pdf',fig_location,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
            close all;
            % create axx_rca_full
            fullW = [rca_full(c,1).W, wComparison];

            axx_rca_full(c) = axxRCAmake(path_names,fullW,curCond,opt.axxType); % create struct

            %axx_rca_full(c).Projected = rcaProject(axx_rca_full(c).Wave,axx_rca_full(c).W);
            % oddball RCA, first two frequencies 1F1 and 2F1 
            % since this is just a subset of the previous RCA, set forceSourceData to false
            rca_odd(c, :) = rcaSweep(path_names,use_bins,use_freqs(1:end/2),curCond,trialsToUse,n_reg,n_comp,data_type,opt.chanToCompare,[],rcPlotStyle,false);
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/oddRCA_cond%.0d&%.0d_cov.pdf',fig_location,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
            close all;
            % carrier RCA, last two frequencies 1F2 and 2F2
            rca_carr(c, :) = rcaSweep(path_names,use_bins,use_freqs(end/2+1:end),curCond,trialsToUse,n_reg,n_comp,data_type,opt.chanToCompare,[],rcPlotStyle,false);
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/carrierRCA_cond%.0d&%.0d_cov.pdf',fig_location,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
            oddW = [rca_odd(c).W, wComparison];
            carrierW = [rca_carr(c).W, wComparison];
            
            axx_rca_odd(c) = axxRCAmake(path_names,oddW,curCond,opt.axxType); % create struct
            axx_rca_carr(c) = axxRCAmake(path_names,carrierW,curCond,opt.axxType); % create struct
            close all;
        end
        save(saveFileName,'rca_full', 'rca_odd', 'rca_carr', 'axx_rca_*','-append')
        warning('on','all')
    end
    
    %% COMPUTE VALUES FOR PLOTTING
    % AE RCA
    keepConditions = true;
    errorType = 'SEM';
    doNR = []; % 1 freq, 5 RCs, 6 conditions
    for d = 1:3 % compute values for both full RCA and merge the split oddball/carrier  
        if d == 1
            curRCA = rca_full;
            curAxxRCA = axx_rca_full;
        elseif d == 2
            curRCA = rca_odd;
            curAxxRCA = axx_rca_odd;
        else
            curRCA = rca_carr;
            curAxxRCA = axx_rca_carr;
        end
        for c = 1:size(curRCA,1)
            for f = 1:size(curRCA,2)
                fprintf('\n ... freq no. %0.0f ...\n',f);
                rcStruct = aggregateData(curRCA(c,f),keepConditions,errorType,opt.trialError,doNR);
                % RC
                curRCA(c,f).stats.Amp = squeeze(rcStruct.ampBins(1,:,:,:));
                curRCA(c,f).stats.SubjectAmp = squeeze(rcStruct.subjectAmp(1,:,:,:,:));
                curRCA(c,f).stats.ErrLB = squeeze(rcStruct.ampErrBins(1,:,:,:,1));
                curRCA(c,f).stats.ErrUB = squeeze(rcStruct.ampErrBins(1,:,:,:,2));
                curRCA(c,f).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins(1,:,:,:));
                curRCA(c,f).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise(1,:,:,:,:));
                % within group t-values
                curRCA(c,f).stats.tSqrdP = squeeze(rcStruct.tSqrdP(1,:,:,:));
                curRCA(c,f).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig(1,:,:,:));
                curRCA(c,f).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal(1,:,:,:));
            end
            % axx_rca_full
            n_tps = size(curAxxRCA(c).Projected{1,1},1);
            n_rca = size(curAxxRCA(c).Projected{1,1},2);
            n_subs = size(curAxxRCA(c).Projected,2);
            n_cond = 2;
            curAxxRCA(c).ProjMat = reshape(cell2mat(curAxxRCA(c).Projected),n_tps,n_cond,n_rca,n_subs);
            curAxxRCA(c).avg = nanmean(curAxxRCA(c).ProjMat,4);
            curAxxRCA(c).std = nanstd(curAxxRCA(c).ProjMat,0,4);
            curAxxRCA(c).sem = curAxxRCA(c).std./sqrt(n_subs);
        end
        if d == 1
            rca_full = curRCA;
            axx_rca_full = curAxxRCA;
        elseif d == 2
            rca_split = curRCA;
            axx_rca_odd = curAxxRCA;
        else
            rca_split = cat(2, rca_split, curRCA);
            axx_rca_carr = curAxxRCA;
        end
    end
        
        
%         for c=1:3
%             
%             [tempDataStrct,tempReal,tempImag] = aggregateData(curRCA.data,curRCA.settings,keepConditions,errorType,opt.trialError);
%             realSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempReal,[2,3,5,4,1]); % frequencies x rcs x condition x subjects x freqpair x rcaType
%             imagSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempImag,[2,3,5,4,1]);
%             ampVals(freqIdx,1:5,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
%             tVs0Stat(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
%             tVs0Pval(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdP;
%             tVs0Sig(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
% %             errLB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
% %             errUB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
%             errLB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,1) );
%             errUB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2) );
%             tempNoiseStrct1 = aggregateData(curRCA.noiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
%             tempNoiseStrct2 = aggregateData(curRCA.noiseData.higherSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
%             [snrTmp,noiseTmp] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
%             snrVals(freqIdx,1:5,:,c,rcaIdx) = squeeze(snrTmp);
%             noiseVals(freqIdx,1:5,:,c,rcaIdx) = squeeze(noiseTmp);
%             % COMPARISON
%             [tempDataStrct,tempReal,tempImag] = aggregateData(curRCA.comparisonData,curRCA.settings,keepConditions,errorType,opt.trialError);
%             realSubjs(freqIdx,6,:,:,c,rcaIdx) = permute(tempReal,[2,3,5,4,1]);
%             imagSubjs(freqIdx,6,:,:,c,rcaIdx) = permute(tempImag,[2,3,5,4,1]);
%             ampVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
%             realVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.realBins );
%             imagVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.imagBins );
% %             errLB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
% %             errUB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
%             errLB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,1) );
%             errUB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2) );
%             tVs0Stat(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
%             tVs0Pval(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdP;
%             tVs0Sig(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
%             tempNoiseStrct1 = aggregateData(curRCA.comparisonNoiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
%             tempNoiseStrct2 = aggregateData(curRCA.comparisonNoiseData.higherSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
%             [snrTmp,noiseTmp] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
%             snrVals(freqIdx,6,:,c,rcaIdx) = squeeze(snrTmp);
%             noiseVals(freqIdx,6,:,c,rcaIdx) = squeeze(noiseTmp);
%             % axx_rca_full
%             nTps = size(curAxxRCA.Projected{1,1},1);
%             nRCA = size(curAxxRCA.Projected{1,1},2);
%             nSubjs = size(curAxxRCA.Projected,2);
%             %projmatTest = reshape(cell2mat(cellfun(@ (x,y) cat(3,x,y), axx_rca_full(3).Projected(1,:),axx_rca_full(3).Projected(2,:),'Uni',false)),840,6,15,2)
%             tmp_axxProjMat(c).ProjMat = reshape(cell2mat(cellfun(@(x,y) cat(3,x,y),curAxxRCA.Projected(1,:),curAxxRCA.Projected(2,:),'Uni',false)),nTps,nRCA,nSubjs,2);
%             tmp_axxProjMat(c).avg = nanmean(tmp_axxProjMat(c).ProjMat,3);
%             tmp_axxProjMat(c).std = nanstd(tmp_axxProjMat(c).ProjMat,0,3);
%             tmp_axxProjMat(c).sem = tmp_axxProjMat(c).std./sqrt(nSubjs);
%         end
%         if d == 1
%             axxProjMat = tmp_axxProjMat;
%         elseif d == 2
%             axxProjMatOdd = tmp_axxProjMat;
%         else
%             axxProjMatCarrier = tmp_axxProjMat;
%         end
%     end

    %% STATS
    for f = 1:max(use_freqs)
        for fPair = 1:3
            for rc_type = 1:2
                for c = 1:2
                    if rc_type == 1
                        [rca_data_real,rca_data_imag] = getRealImag(rca_full(fPair,f).data);
                        [comp_data_real,comp_data_imag] = getRealImag(rca_full(fPair,f).comparisonData);
                    else
                        if f > length(rca_odd)
                            [rca_data_real,rca_data_imag] = getRealImag(rca_carr(fPair,f-length(rca_odd)).data);
                            [comp_data_real,comp_data_imag] = getRealImag(rca_carr(fPair,f-length(rca_odd)).comparisonData);
                        else
                            [rca_data_real,rca_data_imag] = getRealImag(rca_odd(fPair,f).data);
                            [comp_data_real,comp_data_imag] = getRealImag(rca_odd(fPair,f).comparisonData);
                        end
                    end
                    % concatenate comparison and rca data
                    rca_data_real = cellfun(@(x,y) cat(2,x,y), rca_data_real, comp_data_real, 'uni', false);
                    rca_data_imag = cellfun(@(x,y) cat(2,x,y), rca_data_imag, comp_data_imag, 'uni', false);

                    for r = 1:6
                        x_data = cell2mat(cellfun(@(x) squeeze(nanmean(x(1,r,:),3)), rca_data_real, 'uni', false));
                        y_data = cell2mat(cellfun(@(x) squeeze(nanmean(x(1,r,:),3)), rca_data_imag, 'uni', false));
                        xy_data = cat(2,permute(x_data,[2,3,1]),permute(y_data,[2,3,1]));
                        % compute vector mean amplitude, to get sign
                        temp_amp = sqrt(nanmean(x_data,2).^2+nanmean(y_data,2).^2);
                        
                        results = tSquaredFourierCoefs(xy_data);
                        between_t2sig(f,r,fPair,rc_type) = results.H;
                        between_t2p(f,r,fPair,rc_type) = results.pVal;
                        between_t2Stat(f,r,fPair,rc_type) = results.tSqrd;
                        between_sign(f,r,fPair,rc_type) = sign(temp_amp(1)-temp_amp(2));
                    end
                end
            end
        end
    end
    
    [ rcaRelExpl, rcaVarExpl ,pcaVarExpl ] = rcaExplained(rca_full(3,1),2);

    
    %% NEW PLOTTING
    close all
    plotSNR = false;
    x_vals = rca_full(1).settings.freqLabels; % frequencies are x-values
    % Which dimension corresponds to what in axxProjMat
    tPtsIdx = 1;
    condsIdx = 2;
    RCAIdx = 3;
    
    % axxProjMat.avg = axxProjMat(:,newOrder,:);
    % axxProjMat.std = axxProjMat(:,newOrder,:);
    % axxProjMat.sem = axxProjMat(:,newOrder,:);
    % axxProjMat.ProjMat = axxProjMat.ProjMat(:,newOrder,:,:);
    axx_ticks = [0 500 1000 1500 2000];
    bgColors = [1,1,1; 0,0,0];
    subColors =  repmat(distinguishable_colors(4,bgColors),2,1);
    subColors = subColors(1:6,:); %Change this from use_conds for consistency in colors between plots
    lWidth = 1;
    fSize = 10;
    gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

    yFigs = n_comp+1;

    if opt.plotSplit == 0    
        xFigs = 5;
    else
        xFigs = 8;
    end
    figHeight = yFigs * 6;
    figWidth = xFigs * 6;

    clear egiH;
    
    c_brewer = load('colorBrewer_new.mat');
    bold_colors = c_brewer.rgb20([9,3,5,13],:);
    weak_colors = c_brewer.rgb20([10,4,6,14],:);
    barH = [0,0];
    yUnit = 1;
    bar_width = 0.3;
    color_idx = [3,1];
    for f = 1:3 % frequency pairs
        % Axx xVals
        nTps = size(axx_rca_full(f).ProjMat,1);
        axx_xvals = linspace(0, 1000/carrier_freq(f) * axx_reps(f) , nTps+1);
        axx_xvals = axx_xvals(2:end);
        stim_onset = axx_xvals(floor(linspace(1,length(axx_xvals),axx_reps(f) + 1)));
        stim_onset = stim_onset(1:end-1);
        figure;
        for r = 1:yFigs 
            if r < 6            
                if opt.plotSplit == 0
                    curRCA = rca_full(f,:);
                    egiH(r,1) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs);
                    hold on
                    rcaColorBar = [min(curRCA(1).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA(1).A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'full';
                else
                    curRCA = rca_split(f,:);
                    egiH(r,1) = subplot(yFigs,xFigs,3+(r-1)*xFigs); %New Candidate
                    hold on
                    rcaColorBar = [min(curRCA(1).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA(1).A(:,r),rcaColorBar);
                    hold off
                    egiH(r,2) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs); %I think this one is correct
                    hold on
                    rcaColorBar = [min(curRCA(5).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA(5).A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'split';
                end
            else
            end

            % axx plot
            for z = 1:2
                if opt.plotSplit == 0
                    % full axx
                    subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
                    subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs, xFigs+(r-1)*xFigs]);
                    cur_axx = axx_rca_full(f);
                    title_str = 'single cycle waveform (all)';
                else
                    if z == 1
                        subplot(yFigs,xFigs,[1+(r-1)*xFigs, 2+(r-1)*xFigs]);
                        cur_axx = axx_rca_odd(f);
                        title_str = 'single cycle waveform (odd)';
                    else
                        subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs,xFigs+(r-1)*xFigs]);
                        cur_axx = axx_rca_carr(f);
                        title_str = 'single cycle waveform (carrier)';
                    end
                end
                dataToPlot = squeeze(cur_axx.avg(:,:,r));
                errorToPLot = squeeze(cur_axx.sem(:,:,r));
                axx_ymin = floor(min((min(cur_axx.avg(:,:,r)-cur_axx.sem(:,:,r)))/2))*2;
                axx_ymax = ceil(max((max(cur_axx.avg(:,:,r)+cur_axx.sem(:,:,r)))/2))*2;
                if ( axx_ymax-axx_ymin ) > 10
                    axx_yunit = 2;
                else
                    axx_yunit = 1;
                end
                hold on
                for c=1:size(dataToPlot,2)
                    plot(axx_xvals,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(color_idx(c),:));
                    ErrorBars(axx_xvals',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(color_idx(c),:));
                end
                ylim([axx_ymin, axx_ymax]);
                xlim([0,axx_xvals(end)]);
                set(gca,gcaOpts{:},'xtick',axx_ticks,'xticklabel',cellfun(@(x) num2str(x),num2cell(axx_ticks),'uni',false),'ytick',axx_ymin:axx_yunit:axx_ymax);
                yLine = repmat(get(gca,'YLim'),axx_reps(f),1)';
                line(repmat(stim_onset,2,1),yLine,'Color','black');
                if r == 1
                    title(title_str,'fontweight','normal');
                elseif r == 6
                    xlabel('time (ms)');
                else
                end        
                s = get(gca, 'Position');
                set(gca, 'Position', [s(1)+0.01, s(2), s(3)*0.9, s(4)]);
                if opt.plotSplit == 0
                    continue
                else
                end
            end

            for t = 1:2
                subplot(yFigs,xFigs,(xFigs-4)+(r-1)*xFigs+(t-1));
                hold on
                freq_set = max(use_freqs)/2;
                cond_set = 2;
                curIdx = (use_freqs(1:end/2))+(t-1)*freq_set;
                bar_spacing = (0:bar_width:(bar_width*(cond_set-1))) - mean(0:bar_width:(bar_width*(cond_set-1)));
                x_vals = repmat((1:freq_set),cond_set,1) + repmat(bar_spacing,freq_set,1)';
                for c=1:cond_set
                    if plotSNR
                        curRange = snrVals(curIdx,:,:,1,opt.plotSplit+1);
                        %valSet = snrVals(curIdx,r,c,f,opt.plotSplit+1)
                    else
                        amp_vals = arrayfun(@(x) curRCA(x).stats.Amp(r,c), curIdx);
                        err_ub = arrayfun(@(x) curRCA(x).stats.ErrUB(r,c), curIdx);
                        err_lb = arrayfun(@(x) curRCA(x).stats.ErrLB(r,c), curIdx);
                        errorbar(x_vals(c,:), amp_vals, err_lb, err_ub, '.k', 'LineWidth',lWidth, 'marker','none');
                        if c == 1
                            y_max = ceil(max(cell2mat(...
                                arrayfun(@(x) max(curRCA(x).stats.Amp(r,:)+curRCA(x).stats.ErrUB(r,:)), curIdx,'uni',false))));
                        else
                        end
                    end
                    within_sig = arrayfun(@(x) curRCA(x).stats.tSqrdP(r,c)<0.05, curIdx);
                    
                    amp_h(c) = bar(x_vals(c,:), amp_vals,'BarWidth',bar_width,'edgecolor','none','facecolor',bold_colors(color_idx(c),:));

                    for z = 1:length(within_sig)
                        if within_sig(z)
                            patch_x = [ x_vals(c,z)-bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)-bar_width/2];
                            patch_y = [0, 0, y_max, y_max];
                            pa_h = patch(patch_x, patch_y, bold_colors(color_idx(c),:),'edgecolor','none','facealpha',.25); 
                            uistack(pa_h,'bottom');
                            %amp_h(c) = bar(x_vals(c,z), amp_vals(z), 'bar_width',bar_width, 'edgecolor', 'none', 'facecolor', weak_colors(c,:));
                        else
                        end
                    end
                end
                % do between group tests
                if opt.ampTest == 0
                    cur_p = squeeze(between_t2p(curIdx,r,f, opt.plotSplit+1));
                elseif opt.ampTest == 1
                    cur_p = squeeze(AmptPval(curIdx,r,f, opt.plotSplit+1));
                elseif opt.ampTest == 2
                    cur_p = squeeze(zSNRtPval(curIdx,r,f, opt.plotSplit+1));
                end
                cur_sign = squeeze(between_sign(curIdx,r,f, opt.plotSplit+1) == 1);

                between_idx = (cur_p < 0.05) & (cur_sign == 1);
                between_idx = between_idx + ((cur_p < 0.005) & (cur_sign == 1));
                if any(between_idx)
                    arrayfun(@(x) text(x,y_max*.95,'*','fontsize',20,'HorizontalAlignment','center'), find(between_idx > 0),'uni',false);
                    arrayfun(@(x) text(x,y_max*.85,'*','fontsize',20,'HorizontalAlignment','center'), find(between_idx == 2),'uni',false);
                else
                end
                
                if y_max > 2
                    y_unit = 1;
                else
                    y_unit = y_max/4;
                end
                ylim([0,y_max]);
                xlim([.5,freq_set+0.5]);
                if r== 3  && t == 1
                    if plotSNR
                        ylabel('SNR','fontname','Helvetica','fontsize',fSize)
                    else
                        ylabel('amplitude (\muVolts)')
                    end
                else
                end

                if t == 1
                    title_str = 'oddball';
                    harmLabels = {'1F1','2F1','3F1','4F1'};
                else
                    title_str = 'carrier';
                    harmLabels = {'1F2','2F2','3F2','4F2'};
                end
                if r==1
                    title(title_str,'fontweight','normal','fontname','Helvetica','fontsize',fSize);
                elseif r==6
                    if t==2
                        lH = legend(amp_h, cond_names);
                        legend boxoff
                        l_pos = get(lH,'position');
                        l_pos(1) = l_pos(1) + .10;
                        l_pos(2) = l_pos(2) + .05;
                        set(lH,'position',l_pos);
                    else
                    end
                else
                end
                set(gca, gcaOpts{:}, 'xtick', [1,2,3,4], 'xticklabel', harmLabels, 'ytick', 0:y_unit:y_max);
                hold off
            end
        end
        drawnow;
        for r = 1:5
            for z = 1:size(egiH,2)
                if opt.plotSplit
                    addVal = 0.72;
                    shiftLeft = 0.02;
                else
                    addVal = 0.8;
                    shiftLeft = 0;
                end
                newPos = get(egiH(r,z),'position');
                newPos(1) = newPos(1)-(newPos(3)*addVal/2) - shiftLeft;
                newPos(2) = newPos(2)-(newPos(4)*addVal/2);
                newPos(3:4) = newPos(3:4)*(1+addVal);
                set(egiH(r,z),'position',newPos);
            end
        end
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = figHeight;
        figPos(3) = figWidth;
        set(gcf,'pos',figPos);
        % test sufix label
        if opt.ampTest == 0
            testplot = 'vec';
        elseif opt.ampTest == 1
            testplot = 'amp';
        elseif opt.ampTest == 2
            testplot = 'zSNR';
        end
        % behSplit label
        if opt.behSplit
            trialLab = 'corTrials';
        else
            trialLab = 'allTrials';
        end
        if opt.plotSNR        
            fig_name = sprintf('%s/%s_rc%d_%s_snr',fig_location,data_type,f,rcaType);
        else
            fig_name = sprintf('%s/%s_rc%d_%s',fig_location,data_type,f,rcaType);
        end
            %if opt.trainData == 0
            %    fig_name = sprintf('%s/%s_car%d_%s_%s_%s',fig_location,data_type,carriers(f),rcaType,testplot,trialLab);
            %elseif opt.trainData == 1
            %    fig_name = sprintf('%s/%s_car%d_%s_%s_%s_%s',fig_location,data_type,carriers(f),rcaType,train_stim,testplot,trialLab);
            %elseif opt.trainData == 2
            %    fig_name = sprintf('%s/%s_car%d_%s_%s_%s_%s',fig_location,data_type,carriers(f),rcaType,train_stim,testplot,trialLab);
            %else
            %end
            export_fig(sprintf('%s.png', fig_name),'-png','-opengl','-m5','-transparent',gcf);
    end
end
                    
%     %% NEW PLOTTING
%     close all
% %     plotSNR = false;
%     xVals = rca_full(1).settings.rcaFreqs; % frequencies are x-values
%     axxTicks = {[0 250 500 750 1000],[0 250 500 750 1000 1250],[0 500 1000 1500 2000]};
%     lWidth = 1;
%     fSize = 10;
%     gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
% 
%     yFigs = n_comp+1;
% 
%     if opt.plotSplit == 0    
%         xFigs = 5;
%     else
%         xFigs = 8;
%     end
%     figHeight = yFigs * 6;
%     figWidth = xFigs * 6;
% 
%     binVals = rca_full(1).settings.freqLabels';
%     clear egiH;
%     
%     cBrewer = load('colorBrewer_new');
%     bold_colors = cBrewer.rgb20([5,9,13,3],:);
%     weakColors = cBrewer.rgb20([6,10,14,4],:);
%     color_idx = [1,3];
% 
%     barH = [0,0];
%     barWidth = .4;
%     yUnit = 1;
%     
%     axx_ticks = [0 500 1000 1500 2000];
% 
%     
%     for f = 1:3 % frequency pairs
%         % Axx xVals
%         nTps = size(axx_rca_full(1).ProjMat,1);
%         axx_xvals = linspace(0, 1000/carrier_freq(f) * axx_reps(f) , nTps+1);
%         axx_xvals = axx_xvals(2:end);
%         stim_onset = axx_xvals(floor(linspace(1,length(axx_xvals),axx_reps(f) + 1)));
%         stim_onset = stim_onset(1:end-1);
%         figure;
%         for r = 1:yFigs 
%             if r<6            
%                 if opt.plotSplit == 0
%                     curRCA = rca_full(f);
%                     egiH(r,1) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs);
%                     hold on
%                     rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
%                     newExtreme = max(abs(rcaColorBar));
%                     rcaColorBar = [-newExtreme,newExtreme];
%                     mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
%                     hold off
%                     rcaType = 'full';
%                 else
%                     curRCA = rca_odd(f);
%                     egiH(r,1) = subplot(yFigs,xFigs,3+(r-1)*xFigs); %New Candidate
%                     hold on
%                     rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
%                     newExtreme = max(abs(rcaColorBar));
%                     rcaColorBar = [-newExtreme,newExtreme];
%                     mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
%                     hold off
%                     curRCA = rca_carr(f);
%                     egiH(r,2) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs); %I think this one is correct
%                     hold on
%                     rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
%                     newExtreme = max(abs(rcaColorBar));
%                     rcaColorBar = [-newExtreme,newExtreme];
%                     mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
%                     hold off
%                     rcaType = 'split';
%                 end
%             else
%             end
%             
%             % axx plot
%             for z = 1:2
%                 if opt.plotSplit == 0
%                     % full axx
%                     subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
%                     subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs, xFigs+(r-1)*xFigs]);
%                     cur_axx = axx_rca_full(f);
%                     title_str = 'single cycle waveform (all)';
%                 else
%                     if z == 1
%                         subplot(yFigs,xFigs,[1+(r-1)*xFigs, 2+(r-1)*xFigs]);
%                         cur_axx = axx_rca_odd(f);
%                         title_str = 'single cycle waveform (odd)';
%                     else
%                         subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs,xFigs+(r-1)*xFigs]);
%                         cur_axx = axx_rca_carr(f);
%                         title_str = 'single cycle waveform (carrier)';
%                     end
%                 end
%                 dataToPlot = squeeze(cur_axx.avg(:,:,r));
%                 errorToPLot = squeeze(cur_axx.sem(:,:,r));
%                 axx_ymin = floor(min((min(cur_axx.avg(:,:,r)-cur_axx.sem(:,:,r)))/2))*2;
%                 axx_ymax = ceil(max((max(cur_axx.avg(:,:,r)+cur_axx.sem(:,:,r)))/2))*2;
%                 if ( axx_ymax-axx_ymin ) > 10
%                     axx_yunit = 2;
%                 else
%                     axx_yunit = 1;
%                 end
%                 hold on
%                 for c=1:size(dataToPlot,2)
%                     plot(axx_xvals,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(color_idx(c),:));
%                     ErrorBars(axx_xvals',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(color_idx(c),:));
%                 end
%                 ylim([axx_ymin, axx_ymax]);
%                 xlim([0,axx_xvals(end)]);
%                 set(gca,gcaOpts{:},'xtick',axx_ticks,'xticklabel',cellfun(@(x) num2str(x),num2cell(axx_ticks),'uni',false),'ytick',axx_ymin:axx_yunit:axx_ymax);
%                 yLine = repmat(get(gca,'YLim'),axx_reps(f),1)';
%                 line(repmat(stim_onset,2,1),yLine,'Color','black');
%                 if r == 1
%                     title(title_str,'fontweight','normal');
%                 elseif r == 6
%                     xlabel('time (ms)');
%                 else
%                 end        
%                 s = get(gca, 'Position');
%                 set(gca, 'Position', [s(1)+0.01, s(2), s(3)*0.9, s(4)]);
%                 if opt.plotSplit == 0
%                     continue
%                 else
%                 end
%             end
%             for t = 1:2
%                 subplot(yFigs,xFigs,(xFigs-4)+(r-1)*xFigs+(t-1));
%                 hold on
%                 freq_set = max(use_freqs)/2;
%                 cond_set = 2;
%                 curIdx = (use_freqs(1:end/2))+(t-1)*freq_set;
%                 bar_spacing = (0:bar_width:(bar_width*(cond_set-1))) - mean(0:bar_width:(bar_width*(cond_set-1)));
%                 x_vals = repmat((1:freq_set),cond_set,1) + repmat(bar_spacing,freq_set,1)';
%                 for c=1:length(cur_order)
%                     if plotSNR
%                         curRange = snrVals(curIdx,:,cur_order,1,opt.plotSplit+1);
%                         %valSet = snrVals(curIdx,r,c,f,opt.plotSplit+1)
%                       
%                     else
%                         amp_vals = arrayfun(@(x) curRCA(x).stats.Amp(r,cur_order(c)), curIdx);
%                         err_ub = arrayfun(@(x) curRCA(x).stats.ErrUB(r,cur_order(c)), curIdx);
%                         err_lb = arrayfun(@(x) curRCA(x).stats.ErrLB(r,cur_order(c)), curIdx);
%                         errorbar(x_vals(c,:), amp_vals, err_lb, err_ub, '.k', 'LineWidth',lWidth, 'marker','none');
%                         if c == 1
%                             y_max = ceil(max(cell2mat(...
%                                 arrayfun(@(x) max(curRCA(x).stats.Amp(r,cur_order)+curRCA(x).stats.ErrUB(r,cur_order)), curIdx,'uni',false))));
%                         else
%                         end
%                     end
%                     within_sig = arrayfun(@(x) curRCA(x).stats.tSqrdP(r,cur_order(c))<0.05, curIdx);
%                     
%                     amp_h(c) = bar(x_vals(c,:), amp_vals,'BarWidth',bar_width,'edgecolor','none','facecolor',bold_colors(color_idx(c),:));
% 
%                     for z = 1:length(within_sig)
%                         if within_sig(z)
%                             patch_x = [ x_vals(c,z)-bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)-bar_width/2];
%                             patch_y = [0, 0, y_max, y_max];
%                             pa_h = patch(patch_x, patch_y, bold_colors(color_idx(c),:),'edgecolor','none','facealpha',.25); 
%                             uistack(pa_h,'bottom');
%                             %amp_h(c) = bar(x_vals(c,z), amp_vals(z), 'bar_width',bar_width, 'edgecolor', 'none', 'facecolor', weak_colors(c,:));
%                         else
%                         end
%                     end
%                 end
%                 % do between group tests
%                 num_tests = size(between_t2p,4)/length(carriers);
%                 test_idx = (1:num_tests)+(f-1)*num_tests;
%                 if opt.ampTest == 0
%                     cur_p = squeeze(between_t2p(curIdx,r,opt.plotSplit+1,test_idx));
%                 elseif opt.ampTest == 1
%                     cur_p = squeeze(AmptPval(curIdx,r,opt.plotSplit+1,test_idx));
%                 elseif opt.ampTest == 2
%                     cur_p = squeeze(zSNRtPval(curIdx,r,opt.plotSplit+1,test_idx));
%                 end
%                 cur_sign = squeeze(between_sign(curIdx,r,opt.plotSplit+1,test_idx) == 1);
%                 if do_exp == 2
%                     for c = 1:3
%                         between_idx = (cur_p(:,c) < 0.05) & (cur_sign(:,c) == 1);
%                         if any(between_idx)
%                             sig_x = [-.2, .2, 0];
%                             sig_y = [1.05,1.05,.95];
%                             sig_labels = {'1','3','N'};
%                             freq = find(between_idx);
%                             arrayfun(@(x) text(x+sig_x(c), y_max*sig_y(c),sig_labels{c},'fontsize',10,'HorizontalAlignment','center'), freq,'uni',false);
%                         end
%                     end
%                 else
%                     between_idx = (cur_p < 0.05) & (cur_sign == 1);
%                     between_idx = between_idx + ((cur_p < 0.005) & (cur_sign == 1));
%                     if any(between_idx)
%                         arrayfun(@(x) text(x,y_max*.95,'*','fontsize',20,'HorizontalAlignment','center'), find(between_idx > 0),'uni',false);
%                         arrayfun(@(x) text(x,y_max*.85,'*','fontsize',20,'HorizontalAlignment','center'), find(between_idx == 2),'uni',false);
%                     else
%                     end
%                 end
%                 
%                 if y_max > 2
%                     y_unit = 1;
%                 else
%                     y_unit = y_max/4;
%                 end
%                 ylim([0,y_max]);
%                 xlim([.5,freq_set+0.5]);
%                 if r== 3  && t == 1
%                     if plotSNR
%                         ylabel('SNR','fontname','Helvetica','fontsize',fSize)
%                     else
%                         ylabel('amplitude (\muVolts)')
%                     end
%                 else
%                 end
% 
%                 if t == 1
%                     title_str = 'oddball';
%                     harmLabels = {'1F1','2F1','3F1','4F1'};
%                 else
%                     title_str = 'carrier';
%                     harmLabels = {'1F2','2F2','3F2','4F2'};
%                 end
%                 if r==1
%                     title(title_str,'fontweight','normal','fontname','Helvetica','fontsize',fSize);
%                 elseif r==6
%                     if t==2
%                         lH = legend(amp_h, cond_names);
%                         legend boxoff
%                         l_pos = get(lH,'position');
%                         l_pos(1) = l_pos(1) + .10;
%                         l_pos(2) = l_pos(2) + .05;
%                         set(lH,'position',l_pos);
%                     else
%                     end
%                 else
%                 end
%                 set(gca, gcaOpts{:}, 'xtick', [1,2,3,4], 'xticklabel', harmLabels, 'ytick', 0:y_unit:y_max);
%                 hold off
%             end
%         end
%         drawnow;
%         for r = 1:5
%             for z = 1:size(egiH,2)
%                 if opt.plotSplit
%                     addVal = 0.72;
%                     shiftLeft = 0.02;
%                 else
%                     addVal = 0.8;
%                     shiftLeft = 0;
%                 end
%                 newPos = get(egiH(r,z),'position');
%                 newPos(1) = newPos(1)-(newPos(3)*addVal/2) - shiftLeft;
%                 newPos(2) = newPos(2)-(newPos(4)*addVal/2);
%                 newPos(3:4) = newPos(3:4)*(1+addVal);
%                 set(egiH(r,z),'position',newPos);
%             end
%         end
%         set(gcf, 'units', 'centimeters');
%         figPos = get(gcf,'pos');
%         figPos(4) = figHeight;
%         figPos(3) = figWidth;
%         set(gcf,'pos',figPos);
%         % test sufix label
%         if opt.ampTest == 0
%             testplot = 'vec';
%         elseif opt.ampTest == 1
%             testplot = 'amp';
%         elseif opt.ampTest == 2
%             testplot = 'zSNR';
%         end
%         % behSplit label
%         if opt.behSplit
%             trialLab = 'corTrials';
%         else
%             trialLab = 'allTrials';
%         end
%         if plotSNR        
%             fig_name = sprintf('%s/%s_car%d_%s_snr',fig_location,data_type,carriers(f),rcaType);
%         else
%             if opt.trainData == 0
%                 fig_name = sprintf('%s/%s_car%d_%s_%s_%s',fig_location,data_type,carriers(f),rcaType,testplot,trialLab);
%             elseif opt.trainData == 1
%                 fig_name = sprintf('%s/%s_car%d_%s_%s_%s_%s',fig_location,data_type,carriers(f),rcaType,train_stim,testplot,trialLab);
%             elseif opt.trainData == 2
%                 fig_name = sprintf('%s/%s_car%d_%s_%s_%s_%s',fig_location,data_type,carriers(f),rcaType,train_stim,testplot,trialLab);
%             else
%             end
%         end
%         export_fig(sprintf('%s.png', fig_name),'-png','-opengl','-m5','-transparent',gcf);
%     end
% end
% 
% %             %Axx plot
% %             if opt.plotSplit == 0
% %                 % Full Axx
% %                 subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
% %                 subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs, xFigs+(r-1)*xFigs]);
% %                 dataToPlot = squeeze(axx_rca_full(f).ProjMat(f).avg(:,r,:));
% %                 errorToPLot = squeeze(axxProjMat(f).sem(:,r,:));
% %                 AxxyMin = floor(min((min(axxProjMat(f).avg(:,r,:)-axxProjMat(f).sem(:,r,:)))/2))*2;
% %                 AxxyMax = ceil(max((max(axxProjMat(f).avg(:,r,:)+axxProjMat(f).sem(:,r,:)))/2))*2;
% %                 if ( AxxyMax-AxxyMin ) > 10
% %                     AxxyUnit = 2;
% %                 else
% %                     AxxyUnit = 1;
% %                 end
% % 
% % 
% %                 hold on
% %                 for c=1:2
% %                     AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(c,:));
% %                     ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(c,:));
% % 
% %                     %fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],bold_colors(c,:),'EdgeColor',bold_colors(c,:),'LineWidth',0.2);
% %                     %alpha(0.2);
% %                 end
% %                 ylim([AxxyMin,AxxyMax]);
% %                 xlim([0,xValsAxx(end)]);
% %                 set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',AxxyMin:AxxyUnit:AxxyMax);
% %                 yLine = repmat(get(gca,'YLim'),nReps(f),1)';
% %                 line(repmat(stimOnset,2,1),yLine,'Color','black');
% %                 if r == 1
% %                     titleStr = 'Waveform';
% %                     title(titleStr);
% %                 end
% %                 s = get(gca, 'Position');
% %                 set(gca, 'Position', [s(1)+0.01, s(2), s(3)*0.9, s(4)]);
% %                 hold off
% % 
% %             else
% %                 % Odd
% %                 cur_axxProjMat = axxProjMatOdd(f);
% %                 subplot(yFigs,xFigs,[1+(r-1)*xFigs, 2+(r-1)*xFigs]);
% %                 dataToPlot = squeeze(cur_axxProjMat.avg(:,r,:));
% %                 errorToPLot = squeeze(cur_axxProjMat.sem(:,r,:));
% %                 AxxyMin = floor(min((min(cur_axxProjMat.avg(:,r,:)-cur_axxProjMat.sem(:,r,:)))/2))*2;
% %                 AxxyMax = ceil(max((max(cur_axxProjMat.avg(:,r,:)+cur_axxProjMat.sem(:,r,:)))/2))*2;
% %                 if ( AxxyMax-AxxyMin ) > 10
% %                     AxxyUnit = 2;
% %                 else
% %                     AxxyUnit = 1;
% %                 end
% %                 hold on
% %                 for c=1:2
% %                     AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(c,:));
% %                     ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(c,:));
% %                     %fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],bold_colors(c,:),'EdgeColor',bold_colors(c,:),'LineWidth',0.2);
% %                     %alpha(0.2);
% %                 end
% %                 ylim([AxxyMin,AxxyMax]);
% %                 xlim([0,xValsAxx(end)]);
% %                 set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',AxxyMin:AxxyUnit:AxxyMax);
% %                 yLine = repmat(get(gca,'YLim'),nReps(f),1)';
% %                 line(repmat(stimOnset,2,1),yLine,'Color','black');
% %                 if r == 1
% %                     titleStr = 'Oddball Waveform';
% %                     title(titleStr);
% %                 end
% %                 s = get(gca, 'Position');
% %                 set(gca, 'Position', [s(1), s(2), s(3)*0.9, s(4)]);
% %                 hold off
% % 
% %                 %Carrier
% %                 cur_axxProjMat = axxProjMatCarrier(f);
% %                 subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs,xFigs+(r-1)*xFigs]);
% %                 dataToPlot = squeeze(cur_axxProjMat.avg(:,r,:));
% %                 errorToPLot = squeeze(cur_axxProjMat.sem(:,r,:));
% %                 AxxyMin = floor(min((min(cur_axxProjMat.avg(:,r,:)-cur_axxProjMat.sem(:,r,:)))/2))*2;
% %                 AxxyMax = ceil(max((max(cur_axxProjMat.avg(:,r,:)+cur_axxProjMat.sem(:,r,:)))/2))*2;
% %                 if ( AxxyMax-AxxyMin ) > 10
% %                     AxxyUnit = 2;
% %                 else
% %                     AxxyUnit = 1;
% %                 end
% %                 hold on
% %                 for c=1:2
% %                     AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(c,:));
% %                     ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(c,:));
% %                     %fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],bold_colors(c,:),'EdgeColor',bold_colors(c,:),'LineWidth',0.2);
% %                     %alpha(0.2);
% %                 end
% %                 ylim([AxxyMin,AxxyMax]);
% %                 xlim([0,xValsAxx(end)]);
% %                 set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',AxxyMin:AxxyUnit:AxxyMax);
% %                 yLine = repmat(get(gca,'YLim'),nReps(f),1)';
% %                 line(repmat(stimOnset,2,1),yLine,'Color','black');
% %                 if r == 1
% %                     titleStr = 'Carrier Waveform';
% %                     title(titleStr);
% %                 end
% %                 s = get(gca, 'Position');
% %                 set(gca, 'Position', [s(1), s(2), s(3)*0.9, s(4)]);
% %                 hold off
% %             end
% 
% 
%             for t = 1:2
%                 subplot(yFigs,xFigs,(xFigs-4)+(r-1)*xFigs+(t-1));
%                 numFreqs = max(use_freqs)/2;
%                 curIdx = (use_freqs(1:end/2))+(t-1)*numFreqs;
%                 xVals = repmat((1:numFreqs),n_cond,1) + repmat(linspace(-barWidth/2,barWidth/2,n_cond),numFreqs,1)';
% 
% 
%                 hold on
%                 for c=1:n_cond;
%                     if opt.plotSNR
%                         yVals = snrVals(curIdx,r,c,f,opt.plotSplit+1);
%                         curRange = snrVals(:,r,:,f,opt.plotSplit+1);
%                     else
%                         yVals = ampVals(curIdx,r,c,f,opt.plotSplit+1);
%                         curRange = ampVals(:,r,:,f,opt.plotSplit+1);
%                         errorbar(xVals(c,:),yVals,errLB(curIdx,r,c,f,opt.plotSplit+1),errUB(curIdx,r,c,f,opt.plotSplit+1),'.','Color','k','LineWidth',lWidth,'marker','none');
%                     end
%                     zeroSig = tVs0Pval(curIdx,r,c,f,opt.plotSplit+1)<0.05;
%                     
%                     barH(c) = bar(xVals(c,:),yVals,'BarWidth',barWidth,'edgecolor','none','facecolor',bold_colors(c,:));
%                     for z = 1:length(zeroSig)
%                         if ~zeroSig(z)
%                             bar(xVals(c,z),yVals(z),'BarWidth',barWidth,'edgecolor','none','facecolor',weakColors(c,:));
%                         else
%                         end
%                     end
%                 end
%                 
%                 yMax = ceil(max(curRange(:))) + yUnit;
%                 ylim([0,yMax]);
%                 xlim([.5,numFreqs+0.5]);
%                 
%                 curSig = between_t2p(curIdx,r,f,opt.plotSplit+1)<0.05;
%                 curSig = curSig+(between_t2p(curIdx,r,f,opt.plotSplit+1)<0.005);
%                 if any(curSig)
%                     arrayfun(@(x) text(x,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center'), find(curSig >0),'uni',false);
%                     arrayfun(@(x) text(x,yMax*.85,'*','fontsize',20,'HorizontalAlignment','center'), find(curSig==2),'uni',false);
%                 else
%                 end
%                 
%                 if r== 3  && t == 1
%                     if opt.plotSNR
%                         ylabel('SNR','fontname','Helvetica','fontsize',fSize)
%                     else
%                         ylabel('Amplitude (\muVolts)','fontname','Helvetica','fontsize',fSize)
%                     end
%                 else
%                 end
% 
%                 if t == 1
%                     titleStr = 'Oddball';
%                     harmLabels = {'1F1','2F1','3F1','4F1'};
%                 else
%                     titleStr = 'Carrier';
%                     harmLabels = {'1F2','2F2','3F2','4F2'};
%                 end
%                 if r==1
%                     title(titleStr,'fontname','Helvetica','fontsize',fSize);
%                 elseif r==6
%                     if t==2
%                         lH = legend(barH,{'number','control'});
%                         legend boxoff
%                         lPos = get(lH,'position');
%                         lPos(1) = lPos(1) + .18;
%                         lPos(2) = lPos(2) + .04;
%                         set(lH,'position',lPos,'fontsize',14,'fontname','Helvetica');
%                     else
%                     end
%                 else
%                 end
%                 set(gca,gcaOpts{:},'xtick',[1,2,3,4],'xticklabel',harmLabels,'ytick',0:yUnit:yMax);
%                 hold off
%             end
%         end
%         drawnow;
%         for r = 1:5;
%             for z = 1:size(egiH,2);
%                 if opt.plotSplit
%                     addVal = 0.8;
%                     shiftLeft = 0.02;
%                 else
%                     addVal = 0.72;
%                     shiftLeft = 0;
%                 end
%                 newPos = get(egiH(r,z),'position');
%                 newPos(1) = newPos(1)-(newPos(3)*addVal/2) - shiftLeft;
%                 newPos(2) = newPos(2)-(newPos(4)*addVal/2);
%                 newPos(3:4) = newPos(3:4)*(1+addVal);
%                 set(egiH(r,z),'position',newPos);
%             end
%         end
%         set(gcf, 'units', 'centimeters');
%         figPos = get(gcf,'pos');
%         figPos(4) = figHeight;
%         figPos(3) = figWidth;
%         set(gcf,'pos',figPos);
%         if opt.plotSNR        
%             export_fig(sprintf('%s/%s_rc%d_%s_snr.png',fig_location,data_type,f,rcaType),'-png','-opengl','-m5','-transparent',gcf);
%         else
%             export_fig(sprintf('%s/%s_rc%d_%s.png',fig_location,data_type,f,rcaType),'-png','-opengl','-m5','-transparent',gcf);
%         end
%     end
