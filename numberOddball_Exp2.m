function numberOddball_Exp2(varargin)
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
    %                   	true/[false]
    %
    %   trialError      - compute trial-wise error true/[false]
    %
    %   plotSplit       - plot RC run separately on odd and carrier
    %                       true/[false]
    %
    %   forceSourceData - force reading in of raw data when initial running
    %                       RCA true/[false]
    %
    %   axxType         - string indicating the type of axx data to plot
    %                      ['ALL']/'NF1'
    %
    %   trainData       - integer indicating whether to train RCA on all data 
    %                     or only on conditions with carrier = 6 or carrier = 8
    %                     [0](all)/1(6)/2(8)
    %
    %   ampTest         - integer indicating how to test amplitudes
    %                     between conditions: 
    %                     [0](vector)/1(amplitude only)/2(sSNR)
    %
    %   behSplit        - only use correct trials
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
            'trainData', 0, ...
            'ampTest', 0, ...
            'behSplit', false ...
            );

    %% VARIABLES THAT CAN BE CHANGED
    topFolder = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/numeroOddball';
    doExp = 2;

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
    trainData = opt.trainData;
    if trainData == 1
        condsToUse = 4:6;
    elseif trainData == 2
        condsToUse = 1:3;
    end
    trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    nFreq = length(freqsToUse);
    nCond = length(condsToUse);
    dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
    rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

    %% LAUNCH RCA OR LOAD DATA
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
        for c = 1:(length(condsToUse)/length(condsToUse)) % This is not doing anything but stays in case we want to loop and separate data differently
            curCond = condsToUse;
            % full RCA, use all harmonics 1F1, 2F1, 1F2 and 2F2
            if c==1
                fullRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse,curCond,trialsToUse,nReg,nComp,dataType,opt.chanToCompare,[],rcPlotStyle,opt.forceSourceData); % TRUE
            else
                % since this is just a subset of the previous RCA, set forceSourceData to false
                fullRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse,curCond,trialsToUse,nReg,nComp,dataType,opt.chanToCompare,[],rcPlotStyle,false);
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
    
    %% NOW SPLIT BASED ON BEHAVIOR
    if opt.behSplit
        
        % load trial idx, and convert it to a more manageable variable
        load(sprintf('%s/trialIdx.mat',dataLocation));
        if ~iscell(TrialIdx)
            TrialIdx = permute(TrialIdx,[1,3,2]);
            tSize = size(TrialIdx);
            trialIdx = cell(size(fullRCA.data));
            if any(size(trialIdx) ~= tSize(2:end));
                error('mismatch between index and data');
            else
            end
            for q = 1 : prod(tSize(2:end))
                trialIdx{q} = TrialIdx(:,q);
            end
            clear TrialIdx;
        else
        end
    
        [fullRCA,fullRCA_incorrect] = splitRCA(fullRCA, trialIdx);
        [oddRCA,oddRCA_incorrect] = splitRCA(oddRCA, trialIdx);
        [carrierRCA,carrierRCA_incorrect] = splitRCA(carrierRCA, trialIdx);
    end

    %% COMPUTE VALUES FOR PLOTTING
    % AE RCA
    keepConditions = true;
    errorType = 'SEM';
    for d = 1:3 % compute values for both full RCA and merge the split oddball/carrier
        for c=length(fullRCA)
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
            realSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempReal,[2,3,5,4,1]);
            imagSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempImag,[2,3,5,4,1]);
            ampVals(freqIdx,1:5,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
            zSNRVals(freqIdx,1:5,:,:,c,rcaIdx) = permute(squeeze(tempDataStrct.zSNR.subj),[1,2,4,3]);
            tVs0Stat(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
            tVs0Pval(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdP;
            tVs0Sig(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
            errLB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
            errUB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
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
            zSNRVals(freqIdx,6,:,:,c,rcaIdx) = permute(squeeze(tempDataStrct.zSNR.subj),[1,3,2]);
            realVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.realBins );
            imagVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.imagBins );
            errLB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
            errUB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
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
            nConds = size(curAxxRCA.Projected,1);
            %projmatTest = reshape(cell2mat(cellfun(@ (x,y) cat(3,x,y), axxRCA(3).Projected(1,:),axxRCA(3).Projected(2,:),'Uni',false)),840,6,15,2)
    %         tmp_axxProjMat(c).ProjMat = reshape(cell2mat(cellfun(@(x,y) cat(3,x,y),curAxxRCA.Projected(1,:),curAxxRCA.Projected(2,:),'Uni',false)),nTps,nRCA,nSubjs,2);
            tmp_axxProjMat(c).ProjMat = reshape(cell2mat(curAxxRCA.Projected),nTps,nCond,nRCA,nSubjs);
            tmp_axxProjMat(c).avg = nanmean(tmp_axxProjMat(c).ProjMat,4);
            tmp_axxProjMat(c).std = nanstd(tmp_axxProjMat(c).ProjMat,0,4);
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
            for fPair = 1 %was 1:3 before
                for rcType = 1:2
                    for c = 1:nCond %was 2 before
                        tempReal = squeeze(realSubjs(f,r,c,:,fPair,rcType));
                        tempImag = squeeze(imagSubjs(f,r,c,:,fPair,rcType));
                        if opt.trialError
                            tempReal = nanmean(reshape(tempReal,[],numSubs),1)';
                            tempImag = nanmean(reshape(tempImag,[],numSubs),1)';
                        else
                        end
                        xyData(:,:,c) = [tempReal,tempImag];
                        ampData(:,c) = sqrt(tempReal.^2+tempImag.^2);
                        zSNRData(:,c) = squeeze(zSNRVals(f,r,c,:,fPair,rcType));
                    end
                    % think how to compute the appropriate pairs to do the
                    % ttests and how they are stored
                    contrastOrder = [[4 5];[6 5];[6 4];[3 2];[1 2];[1 3]];
                    storeOrder = [[1 1];[2 1];[3 1];[1 2];[2 2];[3 2]];
                    if trainData == 1
                        contrastOrder = [[1 2];[3 2];[3 1]];
                    elseif trainData == 2
                        contrastOrder = [[3 2];[1 2];[1 3]];
                    else
                    end
                    for t=1:nCond % number of tTests
                        results = tSquaredFourierCoefs(xyData(:,:,contrastOrder(t,:)));
                        tSig(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = results.H;
                        tPval(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = results.pVal;
                        tStat(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = results.tSqrd;
                        % amplitude paired ttest
                        [h,p,ci,ampResults] = ttest(ampData(:,contrastOrder(t,1)),ampData(:,contrastOrder(t,2)));
                        AmptSig(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = h;
                        AmptPval(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = p;
                        AmptStat(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = ampResults.tstat;
                        %zSNR paired ttest
                        [h,p,ci,zSNRResults] = ttest(zSNRData(:,contrastOrder(t,1)),zSNRData(:,contrastOrder(t,2)));
                        zSNRtSig(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = h;
                        zSNRtPval(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = p;
                        zSNRtStat(f,r,fPair,rcType,storeOrder(t,1),storeOrder(t,2)) = zSNRResults.tstat;
                    end

                end
            end
        end
    end

    %% NEW PLOTTING
    close all
    plotSNR = false;
    xVals = fullRCA(1).settings.freqsToUse; % frequencies are x-values
    % Which dimension corresponds to what in axxProjMat
    tPtsIdx = 1;
    condsIdx = 2;
    RCAIdx = 3;
    %Reorder to make plotting easier
    carriers = [6, 8];
    newOrders = [[5 4 6]; [2 3 1]]; %6v6,6v5,6v9; 8v8,8v9,8v5
    if trainData == 1
        newOrders = [[2 1 3]];
        carriers = [6];
        trainStim = 'train6';
    elseif trainData == 2
        newOrders = [[2 3 1]];
        carriers = [8];
        trainStim = 'train8';
    else
    end
    % axxProjMat.avg = axxProjMat(:,newOrder,:);
    % axxProjMat.std = axxProjMat(:,newOrder,:);
    % axxProjMat.sem = axxProjMat(:,newOrder,:);
    % axxProjMat.ProjMat = axxProjMat.ProjMat(:,newOrder,:,:);
    axxTicks = {[0 500 1000 1500 2000], [0 500 1000 1500 2000]};
    bgColors = [1,1,1; 0,0,0];
    subColors =  repmat(distinguishable_colors(4,bgColors),2,1);
    subColors = subColors(1:6,:); %Change this from condsToUse for consistency in colors between plots
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
    
    numConds = size(newOrders(1,:),2);
    cBrewer = load('colorBrewer');
    boldColors = cBrewer.rgb20([9,3,5,13],:);
    weakColors = cBrewer.rgb20([10,4,6,14],:);
    barH = [0,0];
    barWidth = .3;
    yUnit = 1;
    
    for f = 1:length(carriers) % carriers (i.e. 6 & 8)
        % Axx xVals
        nTps = size(axxProjMat.ProjMat,tPtsIdx);
        carrierFreq = [3.0 3.0];
        nReps = [6 6];
        xValsAxx = linspace(0,1000/carrierFreq(f)*nReps(f),nTps+1);
        xValsAxx = xValsAxx(2:end);
        stimOnset = xValsAxx(floor(linspace(1,length(xValsAxx),nReps(f)+1)));
        stimOnset = stimOnset(1:end-1);
        figure;
        for r = 1:yFigs 
            if r<6            
                if opt.plotSplit == 0
                    curRCA = fullRCA;
                    egiH(r,1) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs);
                    hold on
                    rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'full';



                else
                    curRCA = oddRCA;
                    egiH(r,1) = subplot(yFigs,xFigs,3+(r-1)*xFigs); %New Candidate
                    hold on
                    rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                    hold off
                    curRCA = carrierRCA;
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
                dataToPlot = squeeze(axxProjMat.avg(:,newOrders(f,:),r));
                errorToPLot = squeeze(axxProjMat.sem(:,newOrders(f,:),r));
                AxxyMin = floor(min((min(axxProjMat.avg(:,newOrders(f,:),r)-axxProjMat.sem(:,newOrders(f,:),r)))/2))*2;
                AxxyMax = ceil(max((max(axxProjMat.avg(:,newOrders(f,:),r)+axxProjMat.sem(:,newOrders(f,:),r)))/2))*2;
                if ( AxxyMax-AxxyMin ) > 10
                    AxxyUnit = 2;
                else
                    AxxyUnit = 1;
                end


                hold on
                for c=1:3
                    AxxH(r) = plot(xValsAxx,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',boldColors(c,:));
                    ErrorBars(xValsAxx',dataToPlot(:,c),errorToPLot(:,c),'color',boldColors(c,:));

%                     fill([(xValsAxx)';flipud((xValsAxx)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],boldColors(c,:),'EdgeColor',boldColors(c,:),'LineWidth',0.2);
%                     alpha(0.2);
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
                cur_axxProjMat = axxProjMatOdd;
                subplot(yFigs,xFigs,[1+(r-1)*xFigs, 2+(r-1)*xFigs]);
                dataToPlot = squeeze(cur_axxProjMat.avg(:,newOrders(f,:),r));
                errorToPLot = squeeze(cur_axxProjMat.sem(:,newOrders(f,:),r));
                AxxyMin = floor(min((min(cur_axxProjMat.avg(:,newOrders(f,:),r)-cur_axxProjMat.sem(:,newOrders(f,:),r)))/2))*2;
                AxxyMax = ceil(max((max(cur_axxProjMat.avg(:,newOrders(f,:),r)+cur_axxProjMat.sem(:,newOrders(f,:),r)))/2))*2;
                if ( AxxyMax-AxxyMin ) > 10
                    AxxyUnit = 2;
                else
                    AxxyUnit = 1;
                end
                hold on
                for c=1:3
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
                cur_axxProjMat = axxProjMatCarrier;
                subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs,xFigs+(r-1)*xFigs]);
                dataToPlot = squeeze(cur_axxProjMat.avg(:,newOrders(f,:),r));
                errorToPLot = squeeze(cur_axxProjMat.sem(:,newOrders(f,:),r));
                AxxyMin = floor(min((min(cur_axxProjMat.avg(:,newOrders(f,:),r)-cur_axxProjMat.sem(:,newOrders(f,:),r)))/2))*2;
                AxxyMax = ceil(max((max(cur_axxProjMat.avg(:,newOrders(f,:),r)+cur_axxProjMat.sem(:,newOrders(f,:),r)))/2))*2;
                if ( AxxyMax-AxxyMin ) > 10
                    AxxyUnit = 2;
                else
                    AxxyUnit = 1;
                end
                hold on
                for c=1:3
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
                xVals = repmat((1:numFreqs),numConds,1) + repmat(linspace(-barWidth,barWidth,numConds),numFreqs,1)';

                hold on
                for c=1:numConds
                    curOrder = newOrders(f,:);
                    if plotSNR
                        curRange = snrVals(curIdx,:,curOrder,1,opt.plotSplit+1);
                        %valSet = snrVals(curIdx,r,c,f,opt.plotSplit+1);
                        yVals = snrVals(curIdx,r,curOrder(c),1,opt.plotSplit+1);
                    else
                        curRange = ampVals(curIdx,r,curOrder,1,opt.plotSplit+1) + errUB(curIdx,r,curOrder,1,opt.plotSplit+1);
                        %valSet = ampVals(:,r,c,f,opt.plotSplit+1);
                        yVals = ampVals(curIdx,r,curOrder(c),1,opt.plotSplit+1);
                        %plot(1:2,noiseVals(curIdx,r,c,f,opt.plotSplit+1),'sq','Color',subColors(c,:),'MarkerSize',5);
                        errorbar(xVals(c,:),yVals,errLB(curIdx,r,curOrder(c),1,opt.plotSplit+1),errUB(curIdx,r,curOrder(c),1,opt.plotSplit+1),'.','Color','k','LineWidth',lWidth,'marker','none');
                    end
                    
                    yMax = ceil(max(curRange(:)));
                    zeroSig = tVs0Pval(curIdx,r,curOrder(c),1,opt.plotSplit+1)<0.05;
                    
                    barH(c) = bar(xVals(c,:),yVals,'BarWidth',barWidth,'edgecolor','none','facecolor',boldColors(c,:));
                    
                    for z = 1:length(zeroSig)
                        if ~zeroSig(z)
                            bar(xVals(c,z),yVals(z),'BarWidth',barWidth,'edgecolor','none','facecolor',weakColors(c,:));
                        else
                        end
                    end
                    
                    if opt.ampTest == 0
                        curSig = tPval(curIdx,r,1,opt.plotSplit+1,:,f)<0.05;
                    elseif opt.ampTest == 1
                        curSig = AmptPval(curIdx,r,1,opt.plotSplit+1,:,f)<0.05;
                    elseif opt.ampTest == 2
                        curSig = zSNRtPval(curIdx,r,1,opt.plotSplit+1,:,f)<0.05;
                    end
                    
                    %curSig = curSig+(tPval(curIdx,r,1,opt.plotSplit+1,:,f)<0.005);
                    curSig = squeeze(curSig);
%                     if any(any(curSig))
%                         [freq,test] = find(curSig >0);
%                         arrayfun(@(x) text(x,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center'), freq,'uni',false);
%                         [freq,test] = find(curSig==2);
%                         arrayfun(@(x) text(x,yMax*.85,'*','fontsize',20,'HorizontalAlignment','center'), freq,'uni',false);
%                     else
%                     end
                    if any(curSig(:,c))
                        if c==1
                            freq = find(curSig(:,c) >0);
                            arrayfun(@(x) text(x-.3,yMax*.85,'1','fontsize',10,'HorizontalAlignment','center'), freq,'uni',false);
                        elseif c==2
                            freq = find(curSig(:,c) >0);
                            arrayfun(@(x) text(x,yMax*.85,'3','fontsize',10,'HorizontalAlignment','center'), freq,'uni',false);
                        elseif c==3
                            freq = find(curSig(:,c) >0);
                            arrayfun(@(x) text(x+.3,yMax*.85,'N','fontsize',10,'HorizontalAlignment','center'), freq,'uni',false);
                        end
                    end
                end
                
                yUnit = 1;
                ylim([0,yMax]);
                xlim([.5,numFreqs+0.5]);
                if r== 3  && t == 1
                    if plotSNR
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
                        lH = legend(barH,{'control','dist 1','dist 3'});
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
                    addVal = 0.72;
                    shiftLeft = 0.02;
                else
                    addVal = 1;
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
        if plotSNR        
            export_fig(sprintf('%s/%s_car%d_%s_snr.pdf',figureLocation,dataType,carriers(f),rcaType),'-pdf','-transparent',gcf);
        else
            if trainData == 0
                export_fig(sprintf('%s/%s_car%d_%s_%s_%s.pdf',figureLocation,dataType,carriers(f),rcaType,testplot,trialLab),'-pdf','-transparent',gcf);
            elseif trainData == 1
                export_fig(sprintf('%s/%s_car%d_%s_%s_%s_%s.pdf',figureLocation,dataType,carriers(f),rcaType,trainStim,testplot,trialLab),'-pdf','-transparent',gcf);
            elseif trainData == 2
                export_fig(sprintf('%s/%s_car%d_%s_%s_%s_%s.pdf',figureLocation,dataType,carriers(f),rcaType,trainStim,testplot,trialLab),'-pdf','-transparent',gcf);
            else
            end
        end
        end
    end
%% Test plot just one condition
% close all
% freqsToUse= 1:8;
% t = 1; % 1 = odd; 2 = carrier 
% numFreqs = max(freqsToUse)/2;
% curIdx = (freqsToUse(1:end/2))+(t-1)*numFreqs;
% f = 1; % carrier either 1 = 6 or 2 = 8
% rc = 2; %rca component number
% newOrders = [[5 4 6]; [2 3 1]]; %6v6,6v5,6v9; 8v8,8v9,8v5
% curOrder = newOrders(f,:);
% hold on
% for c = 1:3
%     curOrder = newOrders(f,:);
%     curRange = ampVals(curIdx,rc,curOrder,1,1) + errUB(curIdx,rc,curOrder,1,1);
%     ampH(c)=plot(1:numFreqs,ampVals(curIdx,rc,curOrder(c),1,opt.plotSplit+1),'-','LineWidth',lWidth,'Color',subColors(c,:));
%     errorbar(1:numFreqs,ampVals(curIdx,rc,curOrder(c),1,1),errLB(curIdx,rc,curOrder(c),1,1),errUB(curIdx,rc,curOrder(c),1,1),'Color',subColors(c,:),'LineWidth',lWidth);
%     yMax = ceil(max(curRange(:)));
%     zeroSig = tVs0Pval(curIdx,rc,curOrder(c),1,1)<0.05;
%     if any(zeroSig)
%         if c==1
%             arrayfun(@(x) ...
%                 text(x-.3,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center','color',subColors(c,:)),...
%                 find(zeroSig==1),'uni',false);
%         elseif c==2
%             arrayfun(@(x) ...
%                 text(x,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center','color',subColors(c,:)),...
%                 find(zeroSig==1),'uni',false);
%         elseif c ==3
%             arrayfun(@(x) ...
%                 text(x+.3,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center','color',subColors(c,:)),...
%                 find(zeroSig==1),'uni',false);
%         end
%     else
%     end
%     curSig = tPval(curIdx,rc,1,1,:,f)<0.05;
%     %curSig = curSig+(tPval(curIdx,r,1,opt.plotSplit+1,:,f)<0.005);
%     curSig = squeeze(curSig);
%     if any(curSig(:,c))
%         if c==1
%             freq = find(curSig(:,c) >0);
%             arrayfun(@(x) text(x-.3,yMax*.85,'A','fontsize',20,'HorizontalAlignment','center'), freq,'uni',false);
%         elseif c==2
%             freq = find(curSig(:,c) >0);
%             arrayfun(@(x) text(x,yMax*.85,'B','fontsize',20,'HorizontalAlignment','center'), freq,'uni',false);
%         elseif c==3
%             freq = find(curSig(:,c) >0);
%             arrayfun(@(x) text(x+.3,yMax*.85,'C','fontsize',20,'HorizontalAlignment','center'), freq,'uni',false);
%         end
%     else
%     end
% end
% hold off
end

function [ trueStruct, falseStruct ] = splitRCA(rcaStruct, splitIdx)
    structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
    noiseVars = {'lowerSideBand','higherSideBand'};
    trueStruct = rcaStruct; falseStruct = rcaStruct;
    for z=1:length(structVars)
        if strfind(lower(structVars{z}),'noise')
            for n = 1:length(noiseVars)
                for f=1:length(rcaStruct)
                    trueStruct(f).(structVars{z}).(noiseVars{n}) = ...
                        cellfun(@(x,y) x(:,:,y), rcaStruct(f).(structVars{z}).(noiseVars{n}), splitIdx,'uni',false);
                    falseStruct(f).(structVars{z}).(noiseVars{n}) = ...
                        cellfun(@(x,y) x(:,:,~y), rcaStruct(f).(structVars{z}).(noiseVars{n}), splitIdx,'uni',false);
                end
            end
        else
            for f=1:length(rcaStruct)
                trueStruct(f).(structVars{z}) = ...
                    cellfun(@(x,y) x(:,:,y), rcaStruct(f).(structVars{z}), splitIdx,'uni',false);
                falseStruct(f).(structVars{z}) = ...
                    cellfun(@(x,y) x(:,:,~y), rcaStruct(f).(structVars{z}), splitIdx,'uni',false);
            end
        end
    end
end
