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
    if ~exist('codeFolder','var')
        codeFolder = '/Users/kohler/code';
        rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
        addpath(genpath(rcaCodePath));
        addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
        addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
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
            'plotSNR', false ...
            );

    %% VARIABLES THAT CAN BE CHANGED
    topFolder = '/Volumes/Denali_DATA1/kohler/EEG_EXP/DATA/numeroOddball';
    doExp = 3;

    %% IDENTIFY DATA LOCATION
    dataLocation = sprintf('%s/Experiment%.0d',topFolder,doExp);
    figureLocation = sprintf('%s/Experiment%.0d/figures',topFolder,doExp);
    if ~exist(figureLocation,'dir')
        mkdir(figureLocation);
    else
    end

    folderNames=subfolders(sprintf('%s/*20*',dataLocation),1);
    if ~isempty(opt.foldersToUse)
        folderNames = folderNames(opt.foldersToUse);
    else
    end
    numSubs = size(folderNames,1);
    % and string for saving the data
    saveStr = datestr(clock,26);
    saveStr(strfind(saveStr,'/')) ='';
    saveFileName = sprintf('%s/rcaData_%s.mat',dataLocation, saveStr); % include the date as a string;

    %% SETUP INPUTS
    binsToUse=0; % use average bin, the zeroeth one
    freqsToUse= 1:8; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
    condsToUse = 1:6;
    trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    nFreq = length(freqsToUse);
    nCond = length(condsToUse);
    dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
    rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

    %% call the function --- USER runs this section (no editing necessary) ---
    % generate "filter" for comparison data
    wComparison=zeros(128,1); wComparison(opt.chanToCompare)=1; 

    if ~opt.launchRCA % if users says analysis is not to be launched
        tempFile = dir([dataLocation,'/rcaData_*.mat']);
        saveFileName = fullfile(dataLocation,tempFile(end).name(1:(end-4))); % grab newest file and loose .mat suffix;
        load(saveFileName);
    else
        for f = 1:length(folderNames)
            tempFolders = subfolders(folderNames{f},1);
            pathNames{f} = sprintf('%s/Exp_TEXT_HCN_128_Avg',tempFolders{end});
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
                fullRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse,curCond,trialsToUse,nReg,nComp,dataType,opt.chanToCompare,[],rcPlotStyle,opt.forceSourceData); % TRUE
            else
                % since this is just a subset of the previous RCA, set forceSourceData to false
                fullRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse,curCond,trialsToUse,nReg,nComp,dataType,opt.chanToCompare,[],rcPlotStyle,false);
            end
             if c==3
               fullRCA(c) = flipSwapRCA(fullRCA(c),1:5,1:4);
           else
           end
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/fullRCA_cond%.0d&%.0d_cov.pdf',figureLocation,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
            close all;
            % create axxRCA
            fullW = [fullRCA(c).W, wComparison];
            if c==1
                axxRCA = axxRCAmake(folderNames,fullW,curCond,c,opt.axxType); % create struct
            else
                axxRCA = axxRCAmake(folderNames,fullW,curCond,c,opt.axxType,axxRCA);
            end
            %axxRCA(c).Projected = rcaProject(axxRCA(c).Wave,axxRCA(c).W);
            % oddball RCA, first two frequencies 1F1 and 2F1 
            % since this is just a subset of the previous RCA, set forceSourceData to false
            oddRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse(1:end/2),curCond,trialsToUse,nReg,nComp,dataType,opt.chanToCompare,[],rcPlotStyle,false);
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/oddRCA_cond%.0d&%.0d_cov.pdf',figureLocation,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
            close all;
            % carrier RCA, last two frequencies 1F2 and 2F2
            carrierRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse(end/2+1:end),curCond,trialsToUse,nReg,nComp,dataType,opt.chanToCompare,[],rcPlotStyle,false);
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/carrierRCA_cond%.0d&%.0d_cov.pdf',figureLocation,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
            oddW = [oddRCA(c).W, wComparison];
            carrierW = [carrierRCA(c).W, wComparison];
            if c==1
                axxRCAodd = axxRCAmake(folderNames,oddW,curCond,c,opt.axxType); % create struct
                axxRCAcarrier = axxRCAmake(folderNames,carrierW,curCond,c,opt.axxType); % create struct
            else
                axxRCAodd = axxRCAmake(folderNames,oddW,curCond,c,opt.axxType,axxRCAodd);
                axxRCAcarrier = axxRCAmake(folderNames,carrierW,curCond,c,opt.axxType,axxRCAcarrier);
            end
            close all;
        end
        save(saveFileName,'fullRCA','oddRCA','carrierRCA','axxRCA*','-append')
        warning('on','all')
    end
    %% rca replaces NaNs with zeroes, correct this
    nanDims = [1,2]; % if all time points are zero, or all channels are zero
    structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
    noiseVars = {'lowerSideBand','higherSideBand'};
    for z=1:length(structVars)
        if strfind(lower(structVars{z}),'noise')
            for n = 1:length(noiseVars)
                for f=1:length(fullRCA)
                    fullRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),fullRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
                    oddRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),oddRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
                    carrierRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),carrierRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
                end
            end
        else
            for f=1:length(fullRCA)
                fullRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),fullRCA(f).(structVars{z}),'uni',false);
                oddRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),oddRCA(f).(structVars{z}),'uni',false);
                carrierRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),carrierRCA(f).(structVars{z}),'uni',false);
            end
        end
    end

    %% COMPUTE VALUES FOR PLOTTING
    % AE RCA
    keepConditions = true;
    errorType = 'SEM';
    for d = 1:3 % compute values for both full RCA and merge the split oddball/carrier
        for c=1:3
            if d == 1
                curRCA = fullRCA(c);
                freqIdx = freqsToUse;
                rcaIdx = 1;
                curAxxRCA = axxRCA(c);
            elseif d == 2
                curRCA = oddRCA(c);
                freqIdx = freqsToUse(1:end/2);
                rcaIdx = 2;
                curAxxRCA = axxRCAodd(c);
            else
                curRCA = carrierRCA(c);
                freqIdx = freqsToUse(end/2+1:end);
                rcaIdx = 2; % we want to store odd and carrier together
                curAxxRCA = axxRCAcarrier(c);
            end
            [tempDataStrct,tempReal,tempImag] = aggregateData(curRCA.data,curRCA.settings,keepConditions,errorType,opt.trialError);
            realSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempReal,[2,3,5,4,1]); % frequencies x rcs x condition x subjects x freqpair x rcaType
            imagSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempImag,[2,3,5,4,1]);
            ampVals(freqIdx,1:5,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
            tVs0Stat(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
            tVs0Pval(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdP;
            tVs0Sig(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
%             errLB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
%             errUB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
            errLB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,1) );
            errUB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2) );
            tempNoiseStrct1 = aggregateData(curRCA.noiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
            tempNoiseStrct2 = aggregateData(curRCA.noiseData.higherSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
            [snrTmp,noiseTmp] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
            snrVals(freqIdx,1:5,:,c,rcaIdx) = squeeze(snrTmp);
            noiseVals(freqIdx,1:5,:,c,rcaIdx) = squeeze(noiseTmp);
            % COMPARISON
            [tempDataStrct,tempReal,tempImag] = aggregateData(curRCA.comparisonData,curRCA.settings,keepConditions,errorType,opt.trialError);
            realSubjs(freqIdx,6,:,:,c,rcaIdx) = permute(tempReal,[2,3,5,4,1]);
            imagSubjs(freqIdx,6,:,:,c,rcaIdx) = permute(tempImag,[2,3,5,4,1]);
            ampVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
            realVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.realBins );
            imagVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.imagBins );
%             errLB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
%             errUB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
            errLB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,1) );
            errUB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2) );
            tVs0Stat(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
            tVs0Pval(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdP;
            tVs0Sig(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
            tempNoiseStrct1 = aggregateData(curRCA.comparisonNoiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
            tempNoiseStrct2 = aggregateData(curRCA.comparisonNoiseData.higherSideBand,curRCA.settings,keepConditions,errorType,opt.trialError);
            [snrTmp,noiseTmp] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
            snrVals(freqIdx,6,:,c,rcaIdx) = squeeze(snrTmp);
            noiseVals(freqIdx,6,:,c,rcaIdx) = squeeze(noiseTmp);
            % AxxRCA
            nTps = size(curAxxRCA.Projected{1,1},1);
            nRCA = size(curAxxRCA.Projected{1,1},2);
            nSubjs = size(curAxxRCA.Projected,2);
            %projmatTest = reshape(cell2mat(cellfun(@ (x,y) cat(3,x,y), axxRCA(3).Projected(1,:),axxRCA(3).Projected(2,:),'Uni',false)),840,6,15,2)
            tmp_axxProjMat(c).ProjMat = reshape(cell2mat(cellfun(@(x,y) cat(3,x,y),curAxxRCA.Projected(1,:),curAxxRCA.Projected(2,:),'Uni',false)),nTps,nRCA,nSubjs,2);
            tmp_axxProjMat(c).avg = nanmean(tmp_axxProjMat(c).ProjMat,3);
            tmp_axxProjMat(c).std = nanstd(tmp_axxProjMat(c).ProjMat,0,3);
            tmp_axxProjMat(c).sem = tmp_axxProjMat(c).std./sqrt(nSubjs);
        end
        if d == 1
            axxProjMat = tmp_axxProjMat;
        elseif d == 2
            axxProjMatOdd = tmp_axxProjMat;
        else
            axxProjMatCarrier = tmp_axxProjMat;
        end
    end

    %% STATS
    for f = 1:max(freqsToUse)
        for r = 1:6
            for fPair = 1:3
                for rcType = 1:2
                    for c = 1:2
                        tempReal = squeeze(realSubjs(f,r,c,:,fPair,rcType));
                        tempImag = squeeze(imagSubjs(f,r,c,:,fPair,rcType));
                        if opt.trialError
                            tempReal = nanmean(reshape(tempReal,[],numSubs),1)';
                            tempImag = nanmean(reshape(tempImag,[],numSubs),1)';
                        else
                        end
                        xyData(:,:,c) = [tempReal,tempImag];
                    end
                    % if r == 1 && fPair == 3 && f == 5 && rcType == 1
                    %     disp('debug');
                    % else
                    % end 
                    results = tSquaredFourierCoefs(xyData);
                    tSig(f,r,fPair,rcType) = results.H;
                    tPval(f,r,fPair,rcType) = results.pVal;
                    tStat(f,r,fPair,rcType) = results.tSqrd;
                end
            end
        end
    end

    %% NEW PLOTTING
    close all
%     plotSNR = false;
    xVals = fullRCA(1).settings.freqsToUse; % frequencies are x-values
    axxTicks = {[0 250 500 750 1000],[0 250 500 750 1000 1250],[0 500 1000 1500 2000]};
    lWidth = 1;
    fSize = 10;
    gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

    yFigs = nComp+1;

    if opt.plotSplit == 0    
        xFigs = 5;
    else
        xFigs = 8;
    end
    figHeight = yFigs * 6;
    figWidth = xFigs * 6;

    binVals = fullRCA(1).settings.freqLabels';
    clear egiH;
    
    numConds = size(ampVals,3);
    cBrewer = load('colorBrewer');
    boldColors = cBrewer.rgb20([5,9,13,3],:);
    weakColors = cBrewer.rgb20([6,10,14,4],:);
    barH = [0,0];
    barWidth = .4;
    yUnit = 1;
    
    for f = 1:3 % frequency pairs
        % Axx xVals
        nTps = size(axxProjMat(f).ProjMat,1);
        carrierFreqs = [6.0 3.75 3.0];
        nReps = [6 5 6];
        xValsAxx = linspace(0,1000/carrierFreqs(f)*nReps(f),nTps+1);
        xValsAxx = xValsAxx(2:end);
        stimOnset = xValsAxx(floor(linspace(1,length(xValsAxx),nReps(f)+1)));
        stimOnset = stimOnset(1:end-1);
        figure;
        for r = 1:yFigs 
            if r<6            
                if opt.plotSplit == 0
                    curRCA = fullRCA(f);
                    egiH(r,1) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs);
                    hold on
                    rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'full';
                else
                    curRCA = oddRCA(f);
                    egiH(r,1) = subplot(yFigs,xFigs,3+(r-1)*xFigs); %New Candidate
                    hold on
                    rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                    hold off
                    curRCA = carrierRCA(f);
                    egiH(r,2) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs); %I think this one is correct
                    hold on
                    rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'split';
                end
            else
            end

            %Axx plot
            if opt.plotSplit == 0
                % Full Axx
                subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
                subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs, xFigs+(r-1)*xFigs]);
                dataToPlot = squeeze(axxProjMat(f).avg(:,r,:));
                errorToPLot = squeeze(axxProjMat(f).sem(:,r,:));
                AxxyMin = floor(min((min(axxProjMat(f).avg(:,r,:)-axxProjMat(f).sem(:,r,:)))/2))*2;
                AxxyMax = ceil(max((max(axxProjMat(f).avg(:,r,:)+axxProjMat(f).sem(:,r,:)))/2))*2;
                if ( AxxyMax-AxxyMin ) > 10
                    AxxyUnit = 2;
                else
                    AxxyUnit = 1;
                end


                hold on
                for c=1:2
                    AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',boldColors(c,:));
                    ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',boldColors(c,:));

                    %fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],boldColors(c,:),'EdgeColor',boldColors(c,:),'LineWidth',0.2);
                    %alpha(0.2);
                end
                ylim([AxxyMin,AxxyMax]);
                xlim([0,xValsAxx(end)]);
                set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',AxxyMin:AxxyUnit:AxxyMax);
                yLine = repmat(get(gca,'YLim'),nReps(f),1)';
                line(repmat(stimOnset,2,1),yLine,'Color','black');
                if r == 1
                    titleStr = 'Waveform';
                    title(titleStr);
                end
                s = get(gca, 'Position');
                set(gca, 'Position', [s(1)+0.01, s(2), s(3)*0.9, s(4)]);
                hold off

            else
                % Odd
                cur_axxProjMat = axxProjMatOdd(f);
                subplot(yFigs,xFigs,[1+(r-1)*xFigs, 2+(r-1)*xFigs]);
                dataToPlot = squeeze(cur_axxProjMat.avg(:,r,:));
                errorToPLot = squeeze(cur_axxProjMat.sem(:,r,:));
                AxxyMin = floor(min((min(cur_axxProjMat.avg(:,r,:)-cur_axxProjMat.sem(:,r,:)))/2))*2;
                AxxyMax = ceil(max((max(cur_axxProjMat.avg(:,r,:)+cur_axxProjMat.sem(:,r,:)))/2))*2;
                if ( AxxyMax-AxxyMin ) > 10
                    AxxyUnit = 2;
                else
                    AxxyUnit = 1;
                end
                hold on
                for c=1:2
                    AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',boldColors(c,:));
                    ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',boldColors(c,:));
                    %fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],boldColors(c,:),'EdgeColor',boldColors(c,:),'LineWidth',0.2);
                    %alpha(0.2);
                end
                ylim([AxxyMin,AxxyMax]);
                xlim([0,xValsAxx(end)]);
                set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',AxxyMin:AxxyUnit:AxxyMax);
                yLine = repmat(get(gca,'YLim'),nReps(f),1)';
                line(repmat(stimOnset,2,1),yLine,'Color','black');
                if r == 1
                    titleStr = 'Oddball Waveform';
                    title(titleStr);
                end
                s = get(gca, 'Position');
                set(gca, 'Position', [s(1), s(2), s(3)*0.9, s(4)]);
                hold off

                %Carrier
                cur_axxProjMat = axxProjMatCarrier(f);
                subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs,xFigs+(r-1)*xFigs]);
                dataToPlot = squeeze(cur_axxProjMat.avg(:,r,:));
                errorToPLot = squeeze(cur_axxProjMat.sem(:,r,:));
                AxxyMin = floor(min((min(cur_axxProjMat.avg(:,r,:)-cur_axxProjMat.sem(:,r,:)))/2))*2;
                AxxyMax = ceil(max((max(cur_axxProjMat.avg(:,r,:)+cur_axxProjMat.sem(:,r,:)))/2))*2;
                if ( AxxyMax-AxxyMin ) > 10
                    AxxyUnit = 2;
                else
                    AxxyUnit = 1;
                end
                hold on
                for c=1:2
                    AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',boldColors(c,:));
                    ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',boldColors(c,:));
                    %fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],boldColors(c,:),'EdgeColor',boldColors(c,:),'LineWidth',0.2);
                    %alpha(0.2);
                end
                ylim([AxxyMin,AxxyMax]);
                xlim([0,xValsAxx(end)]);
                set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',AxxyMin:AxxyUnit:AxxyMax);
                yLine = repmat(get(gca,'YLim'),nReps(f),1)';
                line(repmat(stimOnset,2,1),yLine,'Color','black');
                if r == 1
                    titleStr = 'Carrier Waveform';
                    title(titleStr);
                end
                s = get(gca, 'Position');
                set(gca, 'Position', [s(1), s(2), s(3)*0.9, s(4)]);
                hold off
            end


            for t = 1:2
                subplot(yFigs,xFigs,(xFigs-4)+(r-1)*xFigs+(t-1));
                numFreqs = max(freqsToUse)/2;
                curIdx = (freqsToUse(1:end/2))+(t-1)*numFreqs;
                xVals = repmat((1:numFreqs),numConds,1) + repmat(linspace(-barWidth/2,barWidth/2,numConds),numFreqs,1)';


                hold on
                for c=1:numConds;
                    if opt.plotSNR
                        yVals = snrVals(curIdx,r,c,f,opt.plotSplit+1);
                        curRange = snrVals(:,r,:,f,opt.plotSplit+1);
                    else
                        yVals = ampVals(curIdx,r,c,f,opt.plotSplit+1);
                        curRange = ampVals(:,r,:,f,opt.plotSplit+1);
                        errorbar(xVals(c,:),yVals,errLB(curIdx,r,c,f,opt.plotSplit+1),errUB(curIdx,r,c,f,opt.plotSplit+1),'.','Color','k','LineWidth',lWidth,'marker','none');
                    end
                    zeroSig = tVs0Pval(curIdx,r,c,f,opt.plotSplit+1)<0.05;
                    
                    barH(c) = bar(xVals(c,:),yVals,'BarWidth',barWidth,'edgecolor','none','facecolor',boldColors(c,:));
                    for z = 1:length(zeroSig)
                        if ~zeroSig(z)
                            bar(xVals(c,z),yVals(z),'BarWidth',barWidth,'edgecolor','none','facecolor',weakColors(c,:));
                        else
                        end
                    end
                end
                
                yMax = ceil(max(curRange(:))) + yUnit;
                ylim([0,yMax]);
                xlim([.5,numFreqs+0.5]);
                
                curSig = tPval(curIdx,r,f,opt.plotSplit+1)<0.05;
                curSig = curSig+(tPval(curIdx,r,f,opt.plotSplit+1)<0.005);
                if any(curSig)
                    arrayfun(@(x) text(x,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center'), find(curSig >0),'uni',false);
                    arrayfun(@(x) text(x,yMax*.85,'*','fontsize',20,'HorizontalAlignment','center'), find(curSig==2),'uni',false);
                else
                end
                
                if r== 3  && t == 1
                    if opt.plotSNR
                        ylabel('SNR','fontname','Helvetica','fontsize',fSize)
                    else
                        ylabel('Amplitude (\muVolts)','fontname','Helvetica','fontsize',fSize)
                    end
                else
                end

                if t == 1
                    titleStr = 'Oddball';
                    harmLabels = {'1F1','2F1','3F1','4F1'};
                else
                    titleStr = 'Carrier';
                    harmLabels = {'1F2','2F2','3F2','4F2'};
                end
                if r==1
                    title(titleStr,'fontname','Helvetica','fontsize',fSize);
                elseif r==6
                    if t==2
                        lH = legend(barH,{'number','control'});
                        legend boxoff
                        lPos = get(lH,'position');
                        lPos(1) = lPos(1) + .18;
                        lPos(2) = lPos(2) + .04;
                        set(lH,'position',lPos,'fontsize',14,'fontname','Helvetica');
                    else
                    end
                else
                end
                set(gca,gcaOpts{:},'xtick',[1,2,3,4],'xticklabel',harmLabels,'ytick',0:yUnit:yMax);
                hold off
            end
        end
        drawnow;
        for r = 1:5;
            for z = 1:size(egiH,2);
                if opt.plotSplit
                    addVal = 0.8;
                    shiftLeft = 0.02;
                else
                    addVal = 0.72;
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
        if opt.plotSNR        
            export_fig(sprintf('%s/%s_rc%d_%s_snr.pdf',figureLocation,dataType,f,rcaType),'-pdf','-transparent',gcf);
        else
            export_fig(sprintf('%s/%s_rc%d_%s.pdf',figureLocation,dataType,f,rcaType),'-pdf','-transparent',gcf);
        end
    end
