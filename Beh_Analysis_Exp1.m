close all;
clear all;
dataLocation = '/Volumes/Denali_DATA1/kohler/EEG_EXP/DATA/numeroOddball/Experiment1';
%dataLocation = '~/Documents/Research/Stanford/SSVEP';
folderNames=subfolders(sprintf('%s/20*',dataLocation),1);
for s = 1:length(folderNames)
   tempFolders = subfolders(sprintf('%s/',folderNames{s}),1);
   tempFolders = tempFolders(~ismember(tempFolders, [folderNames{s},'/not_time_corrected']));
   tempFolders = tempFolders(~ismember(tempFolders, [folderNames{s},'/time_corrected']));
   matFiles = subfiles(sprintf('%s/ALL_Exp_MATL_HCN_128_Avg/RT*',tempFolders{end}),1);
   blockNum = 0;
   for m = 1:length(matFiles)
       tmpData = load(matFiles{m});
       if ~isempty(tmpData.CndTiming)
           blockNum = blockNum + 1;
           if s == 1 && blockNum == 1
               numTrials = size(tmpData.TimeLine,1); % trials per block, assume same for all
               conditions = unique(cat(1,tmpData.TimeLine.cndNmb));
           end
           trialIdx = (1:numTrials)+(blockNum-1)*numTrials;
           respData(trialIdx,1,s) = cat(1,tmpData.TimeLine.cndNmb); % condition label
           respData(trialIdx,2,s) = cell2mat(cellfun(@(x) find(ismember({'Mis','Ra','La'},x)),{tmpData.TimeLine.respString},'uni',false))-1; % response (0 = mis, 1 = Ra, 2 = La )
           respData(trialIdx,3,s) = cat(1,tmpData.TimeLine.respTimeSec); % response time
           clear tmpData;
       else
       end
   end
   IDs{s} = folderNames{s}(end-6:end);
   misIdx(:,s) = respData(:,2,s)==0;
   % Change correct response mapping for control and experimental
   % condiitons
   for c=1:length(conditions)
       if mod(c,2) == 1
           corrResp = 1;
       elseif mod(c,2) == 0
           corrResp = 2;
       end
       curIdx = respData(:,1,s) == conditions(c);
       percMis(c,s) = length(find(misIdx(curIdx,s)))./length(find(curIdx));
       aveAcc(c,s) = (sum(respData(~misIdx(:,s) & curIdx,2,s) == corrResp)+0.5)./(length(respData(~misIdx(:,s) & curIdx,2,s) == corrResp)+1);
   end
end
%Put into variables, calculate dprime, FA, TP

%Look at data
mean(aveAcc,1)
mean(aveAcc,2)

FAIdx = [2,4,6]; % False alarm
HrFa = aveAcc;
HrFa(FAIdx,:) = 1 - HrFa(FAIdx,:); % FA = 1-accuracy
Zsc = norminv(HrFa);
dPr = zeros(length(conditions)/2,length(folderNames));
bias = zeros(length(conditions)/2,length(folderNames));
for c=1:length(conditions)/2 % assume freq pairs
    cIdx = (c-1)*2 + 1;
    dPr(c,:) = Zsc(cIdx,:) - Zsc(cIdx+1,:);
    bias(c,:) = -(Zsc(cIdx,:) + Zsc(cIdx+1,:))/2;
end

%For plotting
figureLocation = sprintf('%s/figures',dataLocation);
lWidth = 2;
fSize = 12;
freqLabels = {'6 & 1Hz','3.75 & 0.75Hz','3 & 0.5Hz'};
avgdPr = mean(dPr,2);
avgBias = mean(bias,2);
errdPr = std(dPr,[],2)./sqrt(length(folderNames)-1);
errBias = std(bias,[],2)./sqrt(length(folderNames)-1);

% Individual subject plot
cBrewer = load('colorBrewer_new');
colors = cBrewer.rgb20(round(linspace(1,20,length(IDs))),:);
gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

%dPrime
for z = 1:2
    subplot(2,1,z)
    if z == 1
        errorb(avgdPr,errdPr,'Color',colors(1,:), 'barwidth', 0.25);
        hold on
        pH = plot([1 2 3], avgdPr,'o','markersize',10,'LineWidth',2,'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',[1,1,1]);
        y_lims = [0,4];
        ylabel('dPrime');
    else
        errorb(avgBias,errBias,'Color',colors(1,:), 'barwidth', 0.25);
        hold on
        pH = plot([1 2 3], avgBias,'o','markersize',10,'LineWidth',2,'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',[1,1,1]);
        plot([1 2 3], avgBias,'.','markersize',20,'Color','b');
        y_lims = [-0.5,1];
        ylabel('bias');
    end
    uistack(pH, 'top');
    xlabel('freq pairs');
    set(gca,gcaOpts{:},'xtick',[1,2,3],'xticklabel',freqLabels,'xlim',[0.5,3.5],'ylim',y_lims,'box','off');
    hold off
end
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = 20;
figPos(3) = 20;
set(gcf,'pos',figPos);
export_fig(sprintf('%s/avg_dPrime.pdf',figureLocation),'-pdf','-transparent',gcf);

%% 
xFigs = 2;% dPrime or Bias
yFigs = 3;% Freq pairs
titleStr = {'dPrime', 'bias'};

for c=1:2 %dPrime or Bias
    for f=1:3 % Freq pairs
        subplot(yFigs,xFigs, c+(f-1)*xFigs)
        if c == 1
            dataToPlot = dPr;
            y_min = 0;
            y_max = 5;
        else
            dataToPlot = bias;
            y_min = -1;
            y_max = 2;
        end
        x_min = 0.5;
        x_max = 15.5;
        scatter(1:length(IDs),dataToPlot(f,:),100,colors);
        %plot([1:length(IDs)],dataToPlot,'.','markersize',20);
        hold on
        if c==1
            %ref_h = plot([x_min, x_max],ones(1,2)*errdPr(f)*2,'k-','linewidth',2);
            % uistack(ref_h,'bottom');
        else
            ref_h = plot([x_min, x_max], zeros(1,2),'k-','linewidth',2);
            uistack(ref_h,'bottom');
        end
        
        if f == 1
            title(titleStr{c}, 'fontweight','normal')
        elseif f == 3
            xlabel('participants');
        end
        if c == 1
            ylabel(freqLabels{f})
        end
        yMax = max(max(dataToPlot));
        set(gca,gcaOpts{:},'xtick',[1:length(IDs)], 'xticklabel', {''}, 'xlim',[x_min, x_max],'ylim',[y_min, y_max],'box','off','clipping','off');
        %xtickangle(90);
        hold off
    end
    
end
tightfig;
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = 15;
figPos(3) = 20;
set(gcf,'pos',figPos);
export_fig(sprintf('%s/SubjPerformance.pdf',figureLocation),'-pdf','-transparent',gcf);

% Paired tTests dPr
% dPr(1,:); %6Hz-1Hz
% dPr(2,:); %3.75Hz-1.75Hz
% dPr(3,:); %3Hz-0.5Hz
test_order = [[1 2];[1,3];[2,3]];
for freq=1:3
    x = test_order(freq,1);
    y = test_order(freq,2);
    [h,p,ci,stats] = ttest(dPr(x,:),dPr(y,:));
    dPr_res.h(freq) = h;
    dPr_res.p(freq) = p;
    dPr_res.stats(freq) = stats;
end

% Paired tTests Bias
% bias(1,:); %6Hz-1Hz
% bias(2,:); %3.75Hz-1.75Hz
% bias(3,:); %3Hz-0.5Hz
test_order = [[1 2];[1,3];[2,3]];
for freq=1:3
    x = test_order(freq,1);
    y = test_order(freq,2);
    [h,p,ci,stats] = ttest(bias(x,:),bias(y,:));
    bias_res.h(freq) = h;
    bias_res.p(freq) = p;
    bias_res.stats(freq) = stats;
end