close all;
clear all;
dataLocation = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/numeroOddball/Experiment1';
%dataLocation = '~/Documents/Research/Stanford/SSVEP';
folderNames=subfolders(sprintf('%s/20*',dataLocation),1);
for s = 1:length(folderNames)
   tempFolders = subfolders(sprintf('%s/',folderNames{s}),1);
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
mean(aveAcc,1);
mean(aveAcc,2);
median(aveAcc,1);
median(aveAcc,2);
FAIdx = [2,4,6];
HrFa = aveAcc;
HrFa(FAIdx,:) = 1 - HrFa(FAIdx,:);
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
freqLabels = {'6:1Hz','3.75:0.75Hz','3:0.5Hz'};
gcaOpts = {'box','off','fontname','Helvetica','linewidth',lWidth,'box','off'};
avgdPr = mean(dPr,2);
avgBias = mean(bias,2);
errdPr = std(dPr,[],2)./sqrt(length(folderNames)-1);
errBias = std(bias,[],2)./sqrt(length(folderNames)-1);

%dPrime
plot([1 2 3], avgdPr,'.','markersize',20,'Color','b');
hold on
errorbar(avgdPr,errdPr,'Color','b');
ylabel('dPrime');
xlabel('Freq pairs');
set(gca,gcaOpts{:},'xtick',[1,2,3],'xticklabel',freqLabels, 'xlim',[0.5,3.5],'box','off');
hold off
export_fig(sprintf('%s/avg_dPrime.pdf',figureLocation),'-pdf','-transparent',gcf);

%Bias
plot([1 2 3], avgBias,'.','markersize',20,'Color','b');
hold on
errorbar(avgBias,errBias,'Color','b');
ylabel('Bias');
xlabel('Freq pairs');
set(gca,gcaOpts{:},'xtick',[1,2,3],'xticklabel',freqLabels, 'xlim',[0.5,3.5],'box','off');
hold off
export_fig(sprintf('%s/avg_Bias.pdf',figureLocation),'-pdf','-transparent',gcf);

% Individual subject plot
cBrewer = load('colorBrewer');
colors = cBrewer.rgb20(round(linspace(1,20,length(IDs))),:);
yFigs = 2;% dPrime or Bias
xFigs = 3;% Freq pairs
titleStr = {'6:1Hz','3.75:0.75Hz','3:0.5Hz'};
titleStr2 = {'dPrime', 'Bias'};

for f=1:3 % Freq pairs
    for c=1:2 %dPrime or Bias
        subplot(yFigs,xFigs,xFigs*(c-1)+f)
        if c == 1
            dataToPlot = dPr(f,:);
        else
            dataToPlot = bias(f,:);
        end
        scatter(1:length(IDs),dataToPlot,60,colors,'filled');
        %plot([1:length(IDs)],dataToPlot,'.','markersize',20);
        hold on
        if c==1
            hline = refline([0,errdPr(f)*2]); % 2* std error of dprime sample distribution
%             hline = refline([0,2]); %dPrime = 2 ~ p<0.05
            
        else
            hline = refline([0,0]); %bias = 0
        end
        hline.Color = 'k';
        hline.LineWidth = 2;
        yMin = min(dataToPlot);
        yMax = max(dataToPlot);
        set(gca,gcaOpts{:},'xtick',[1:length(IDs)],'xticklabel',IDs, 'xlim',[0.5,length(IDs)+0.5],'ylim',[0-0.5,yMax+0.5],'box','off');
        xtickangle(90);
        if c == 1
            title(titleStr{f})
            if f == 1
                ylabel(titleStr2{c})
            end
        else
            if f == 1
                ylabel(titleStr2{c})
            end
        end
        hold off
    end
end
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