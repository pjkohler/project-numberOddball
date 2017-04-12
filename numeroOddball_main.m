%% add paths
clear all except codeFolder; 
close all;
if ~exist('codeFolder','var')
    codeFolder = '/Users/kohler/code';
    rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
    addpath(genpath(rcaCodePath));
    addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
else
end
setenv('DYLD_LIBRARY_PATH','')

%% VARIABLES THAT CAN BE CHANGED
topFolder = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/numeroOddball';
printFigures=true; % set to true if you want to automatically print the figures and save them into the high-level directory where the data are
launchAnalysis = true;
trialError = false;
forceSourceData = false; % generate source data for first instance of rca?
doExp = 1;

%% IDENTIFY DATA LOCATION
dataLocation = sprintf('%s/Experiment%.0d',topFolder,doExp);
figureLocation = sprintf('%s/Experiment%.0d/figures',topFolder,doExp);
if ~exist(figureLocation,'dir')
    mkdir(figureLocation);
else
end

folderNames=subfolders(sprintf('%s/*20*',dataLocation),1);
folderNames = folderNames(1:10);
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
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

%% call the function --- USER runs this section (no editing necessary) ---
% generate "filter" for comparison data
wComparison=zeros(128,1); wComparison(chanToCompare)=1; 

if ~launchAnalysis % if users says analysis is not to be launched
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
        % (always forceSourceData)
        if c==1
            fullRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse,curCond,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData); % TRUE
        else
            fullRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse,curCond,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData);
        end
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/fullRCA_cond%.0d&%.0d_cov.pdf',figureLocation,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
        close all;
        % create axxRCA
        if c==1
            axxRCA = axxRCAmake(folderNames,fullRCA(c).W,curCond,c); % create struct
        else
            axxRCA = axxRCAmake(folderNames,fullRCA(c).W,curCond,c,axxRCA);
        end
        %axxRCA(c).Projected = rcaProject(axxRCA(c).Wave,axxRCA(c).W);
        % oddball RCA, first two frequencies 1F1 and 2F1 
        % (always forceSourceData)
        oddRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse(1:end/2),curCond,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData);
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/oddRCA_cond%.0d&%.0d_cov.pdf',figureLocation,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
        close all;
        % carrier RCA, last two frequencies 1F2 and 2F2
        carrierRCA(c) = rcaSweep(pathNames,binsToUse,freqsToUse(end/2+1:end),curCond,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData);
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/carrierRCA_cond%.0d&%.0d_cov.pdf',figureLocation,curCond(1),curCond(2)),'-pdf','-transparent',rcaH);
        if c==1
            axxRCAodd = axxRCAmake(folderNames,oddRCA(c).W,curCond,c); % create struct
            axxRCAcarrier = axxRCAmake(folderNames,carrierRCA(c).W,curCond,c); % create struct
        else
            axxRCAodd = axxRCAmake(folderNames,oddRCA(c).W,curCond,c,axxRCAodd);
            axxRCAcarrier = axxRCAmake(folderNames,carrierRCA(c).W,curCond,c,axxRCAcarrier);
        end
        close all;
    end
    save(saveFileName,'fullRCA','oddRCA','carrierRCA','axxRCA','-append')
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
        [tempDataStrct,tempReal,tempImag] = aggregateData(curRCA.data,curRCA.settings,keepConditions,errorType,trialError);
        realSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempReal,[3,4,5,2,1]);
        imagSubjs(freqIdx,1:5,:,:,c,rcaIdx) = permute(tempImag,[3,4,5,2,1]);
        ampVals(freqIdx,1:5,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
        tVs0Stat(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
        tVs0Pval(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdP;
        tVs0Sig(freqIdx,1:5,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
        errLB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
        errUB(freqIdx,1:5,:,c,rcaIdx)   = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
        tempNoiseStrct1 = aggregateData(curRCA.noiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError);
        tempNoiseStrct2 = aggregateData(curRCA.noiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError);
        [snrTmp,noiseTmp] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
        snrVals(freqIdx,1:5,:,c,rcaIdx) = squeeze(snrTmp);
        noiseVals(freqIdx,1:5,:,c,rcaIdx) = squeeze(noiseTmp);
        % COMPARISON
        [tempDataStrct,tempReal,tempImag] = aggregateData(curRCA.comparisonData,curRCA.settings,keepConditions,errorType,trialError);
        realSubjs(freqIdx,6,:,:,c,rcaIdx) = permute(tempReal,[3,4,5,2,1]);
        imagSubjs(freqIdx,6,:,:,c,rcaIdx) = permute(tempImag,[3,4,5,2,1]);
        ampVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins );
        realVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.realBins );
        imagVals(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.imagBins );
        errLB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1) );
        errUB(freqIdx,6,:,c,rcaIdx) = squeeze( tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins );
        tVs0Stat(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdVal;
        tVs0Pval(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdP;
        tVs0Sig(freqIdx,6,:,c,rcaIdx) = tempDataStrct.tSqrdSig;
        tempNoiseStrct1 = aggregateData(curRCA.comparisonNoiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError);
        tempNoiseStrct2 = aggregateData(curRCA.comparisonNoiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError);
        [snrTmp,noiseTmp] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
        snrVals(freqIdx,6,:,c,rcaIdx) = squeeze(snrTmp);
        noiseVals(freqIdx,6,:,c,rcaIdx) = squeeze(noiseTmp);
        % AxxRCA
        nTps = size(curAxxRCA.Projected{1,1},1);
        nRCA = size(curAxxRCA.Projected{1,1},2);
        nSubjs = size(curAxxRCA.Projected,2);
        tmp_axxProjMat(c).ProjMat = reshape(cell2mat(arrayfun(@(x) cat(2,x{:}),curAxxRCA.Projected,'Uni',false)),[nTps,nRCA,2,nSubjs]);
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
        for fPair = 1:3
            for rcType = 1:2
                for c = 1:2
                    tempReal = squeeze(realSubjs(f,r,c,:,fPair,rcType));
                    tempImag = squeeze(imagSubjs(f,r,c,:,fPair,rcType));
                    if trialError
                        tempReal = nanmean(reshape(tempReal,[],numSubs),1)';
                        tempImag = nanmean(reshape(tempImag,[],numSubs),1)';
                    else
                    end
                    xyData(:,:,c) = [tempReal,tempImag];
                end
                if fPair == 1 && rcType == 1 && r == 6 && ( f == 3 ) 
                    results = tSquaredFourierCoefs(xyData(:,:,1),xyData(:,:,2));
                else
                    results = tSquaredFourierCoefs(xyData(:,:,1),xyData(:,:,2));
                end
                tSig(f,r,fPair,rcType) = results.H;
                tPval(f,r,fPair,rcType) = results.pVal;
                tStat(f,r,fPair,rcType) = results.tSqrd;
            end
        end
    end
end

%% NEW PLOTTING
close all
plotSNR = false;
plotSplit = 1; % 0 plot full RCA, 1 plot split RCA 
xVals = fullRCA(1).settings.freqsToUse; % frequencies are x-values
bgColors = [1,1,1; 0,0,0];
subColors =  repmat(distinguishable_colors(4,bgColors),2,1);
subColors = subColors(condsToUse,:);
lWidth = 1;
fSize = 8;
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Arial','linewidth',lWidth};

yFigs = nComp+1;

if plotSplit == 1    
    xFigs = 6;
else
    xFigs = 4;
end
figHeight = yFigs * 6;
figWidth = xFigs * 6;

binVals = fullRCA(1).settings.freqLabels';
clear egiH;
for f = 1:3 % frequency pairs
    figure;
    for r = 1:yFigs 
        if r<6            
            if plotSplit == 0
                curRCA = fullRCA(f);
                egiH(r,1) = subplot(yFigs,xFigs,(xFigs-1)+(r-1)*xFigs);
                hold on
                rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                hold off
                rcaType = 'full';
                
                %Axx plot
                nTps = size(axxProjMat(f).ProjMat,1);
                subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
                dataToPlot = squeeze(axxProjMat(f).avg(:,r,:));
                errorToPLot = squeeze(axxProjMat(f).sem(:,r,:));
                AxxyMin = floor(min(axxProjMat(f).avg(:)-axxProjMat(f).sem(:)));
                AxxyMax = ceil(max(axxProjMat(f).avg(:)+axxProjMat(f).sem(:)));
                AxxyUnit = 2;
                hold on
                for c=1:2
                    AxxH(r) = plot(1:nTps,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',subColors(c,:));
                    fill([(1:nTps)';flipud((1:nTps)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],subColors(c,:),'EdgeColor',subColors(c,:),'LineWidth',0.2);
                    alpha(0.2);
                end
                ylim([AxxyMin,AxxyMax]);
                xlim([-5,nTps+5]);
                set(gca,gcaOpts{:},'xtick',[1,nTps],'xticklabel',{'1',num2str(nTps)},'ytick',AxxyMin:AxxyUnit:AxxyMax);
                if r == 1
                    titleStr = 'Axx Wave';
                    title(titleStr);
                end
                hold off
                
            else
                curRCA = oddRCA(f);
                %egiH(r,1) = subplot(yFigs,xFigs,1+(r-1)*xFigs); %Probably need to change this
                egiH(r,1) = subplot(yFigs,xFigs,2+(r-1)*xFigs); %New Candidate
                hold on
                rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                hold off
                curRCA = carrierRCA(f);
                egiH(r,2) = subplot(yFigs,xFigs,(xFigs-1)+(r-1)*xFigs); %I think this one is correct
                hold on
                rcaColorBar = [min(curRCA.A(:,r)),max(curRCA.A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(curRCA.A(:,r),rcaColorBar);
                hold off
                rcaType = 'split';
                
                %Axx plot
                % Odd
                cur_axxProjMat = axxProjMatOdd(f);
                nTps = size(cur_axxProjMat.ProjMat,1);
                subplot(yFigs,xFigs,1+(r-1)*xFigs);
                dataToPlot = squeeze(cur_axxProjMat.avg(:,r,:));
                errorToPLot = squeeze(cur_axxProjMat.sem(:,r,:));
                AxxyMin = floor(min(cur_axxProjMat.avg(:)-cur_axxProjMat.sem(:)));
                AxxyMax = ceil(max(cur_axxProjMat.avg(:)+cur_axxProjMat.sem(:)));
                AxxyUnit = 2;
                hold on
                for c=1:2
                    AxxH(r) = plot(1:nTps,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',subColors(c,:));
                    fill([(1:nTps)';flipud((1:nTps)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],subColors(c,:),'EdgeColor',subColors(c,:),'LineWidth',0.2);
                    alpha(0.2);
                end
                ylim([AxxyMin,AxxyMax]);
                xlim([-5,nTps+5]);
                set(gca,gcaOpts{:},'xtick',[1,nTps],'xticklabel',{'1',num2str(nTps)},'ytick',AxxyMin:AxxyUnit:AxxyMax);
                if r == 1
                    titleStr = 'Oddball Axx Wave';
                    title(titleStr);
                end
                hold off
                
                %Carrier
                cur_axxProjMat = axxProjMatCarrier(f);
                nTps = size(cur_axxProjMat.ProjMat,1);
                subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
                dataToPlot = squeeze(cur_axxProjMat.avg(:,r,:));
                errorToPLot = squeeze(cur_axxProjMat.sem(:,r,:));
                AxxyMin = floor(min(cur_axxProjMat.avg(:)-cur_axxProjMat.sem(:)));
                AxxyMax = ceil(max(cur_axxProjMat.avg(:)+cur_axxProjMat.sem(:)));
                hold on
                for c=1:2
                    AxxH(r) = plot(1:nTps,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',subColors(c,:));
                    fill([(1:nTps)';flipud((1:nTps)')],[dataToPlot(:,c)-errorToPLot(:,c);flipud(dataToPlot(:,c)+errorToPLot(:,c))],subColors(c,:),'EdgeColor',subColors(c,:),'LineWidth',0.2);
                    alpha(0.2);
                end
                ylim([AxxyMin,AxxyMax]);
                xlim([-5,nTps+5]);
                set(gca,gcaOpts{:},'xtick',[1,nTps],'xticklabel',{'1',num2str(nTps)},'ytick',AxxyMin:AxxyUnit:AxxyMax);
                if r == 1
                    titleStr = 'Carrier Axx Wave';
                    title(titleStr);
                end
                hold off
            end
        else
        end
        for t = 1:2
            subplot(yFigs,xFigs,(xFigs-3)+(r-1)*xFigs+(t-1));
            numFreqs = max(freqsToUse)/2;
            curIdx = (freqsToUse(1:end/2))+(t-1)*numFreqs;
             
            hold on
            for c=1:2
                if plotSNR
                    curRange = snrVals(:,:,:,f,plotSplit+1);
                    %valSet = snrVals(curIdx,r,c,f,plotSplit+1);
                    ampH(c)=plot(1:numFreqs,snrVals(curIdx,r,c,f,plotSplit+1),'-','LineWidth',lWidth,'Color',subColors(c,:));
                else
                    curRange = ampVals(:,:,:,f,plotSplit+1);
                    %valSet = ampVals(:,r,c,f,plotSplit+1);
                    ampH(c)=plot(1:numFreqs,ampVals(curIdx,r,c,f,plotSplit+1),'-','LineWidth',lWidth,'Color',subColors(c,:));
                    %plot(1:2,noiseVals(curIdx,r,c,f,plotSplit+1),'sq','Color',subColors(c,:),'MarkerSize',5);
                    errorbar(1:numFreqs,ampVals(curIdx,r,c,f,plotSplit+1),errLB(curIdx,r,c,f,plotSplit+1),errUB(curIdx,r,c,f,plotSplit+1),'Color',subColors(c,:),'LineWidth',lWidth);
                end
                yMax = ceil(max(curRange(:)));
                zeroSig = tVs0Pval(curIdx,r,c,f,plotSplit+1)<0.05;
                if any(zeroSig)
                    if c==1
                        arrayfun(@(x) ...
                            text(x-.15,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center','color',subColors(c,:)),...
                            find(zeroSig==1),'uni',false);
                    else
                        arrayfun(@(x) ...
                            text(x+.15,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center','color',subColors(c,:)),...
                            find(zeroSig==1),'uni',false);
                    end
                else
                end
            end
            
            curSig = tPval(curIdx,r,f,plotSplit+1)<0.05;
            curSig = curSig+(tPval(curIdx,r,f,plotSplit+1)<0.005);
            if any(curSig)
                arrayfun(@(x) text(x,yMax*.95,'*','fontsize',20,'HorizontalAlignment','center'), find(curSig >0),'uni',false);
                arrayfun(@(x) text(x,yMax*.85,'*','fontsize',20,'HorizontalAlignment','center'), find(curSig==2),'uni',false);
            else
            end
            
            
            yUnit = 1;
            ylim([0,yMax]);
            xlim([.5,numFreqs+0.5]);
            if r== 3  && t == 1
                if plotSNR
                    ylabel('SNR')
                else
                    ylabel('Amplitude (\muVolts)')
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
                title(titleStr);
            elseif r==6
                if t==2
                    lH = legend(ampH,{'number','control'});
                    legend boxoff
                    lPos = get(lH,'position');
                    lPos(1) = lPos(1) + .10;
                    lPos(2) = lPos(2) + .05;
                    set(lH,'position',lPos);
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
            addVal = 1;
            newPos = get(egiH(r,z),'position');
            newPos(1) = newPos(1)-(newPos(3)*addVal/2);
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
    if plotSNR        
        export_fig(sprintf('%s/%s_rc%d_%s_snr.pdf',figureLocation,dataType,f,rcaType),'-pdf','-transparent',gcf);
    else
        export_fig(sprintf('%s/%s_rc%d_%s.pdf',figureLocation,dataType,f,rcaType),'-pdf','-transparent',gcf);
    end
end
