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
lWidth = 1;
freqLabels = {'1:6Hz','0.75:3.75Hz','0.5:3Hz'};
gcaOpts = {'box','off','fontname','Arial','linewidth',lWidth};
avgdPr = mean(dPr,2);
avgBias = mean(bias,2);
errdPr = std(dPr,[],2)./sqrt(length(folderNames)-1);
errBias = std(bias,[],2)./sqrt(length(folderNames)-1);

plot([1 2 3], avgdPr,'-','LineWidth',lWidth);
errorbar(avgdPr,errdPr);
ylabel('dPrime');
xlabel('Freq pairs');
set(gca,gcaOpts{:},'xtick',[1,2,3],'xticklabel',freqLabels, 'xlim',[0.5,3.5]);
export_fig(sprintf('%s/avg_dPrime.pdf',figureLocation));

plot([1 2 3], avgBias,'-','LineWidth',lWidth);
errorbar(avgBias,errBias);
ylabel('Bias');
xlabel('Freq pairs');
set(gca,gcaOpts{:},'xtick',[1,2,3],'xticklabel',freqLabels, 'xlim',[0.5,3.5]);
export_fig(sprintf('%s/avg_Bias.pdf',figureLocation));

xFigs = 2;% dPrime or Bias
yFigs = 3;% Freq pairs
titleStr = {'dPrime', 'Bias'};
for f=1:3 % Freq pairs
    for c=1:2 %dPrime or Bias
        subplot(yFigs,xFigs,xFigs*(f-1)+c)
        if c == 1
            dataToPlot = dPr(f,:);
        else
            dataToPlot = bias(f,:);
        end
        plot([1:length(IDs)],dataToPlot,'.','markersize',20);
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
        set(gca,gcaOpts{:},'xtick',[1:length(IDs)],'xticklabel',IDs, 'xlim',[0.5,length(IDs)+0.5],'ylim',[0-0.5,yMax+0.5]);
        xtickangle(90);
        if f == 1
            title(titleStr{c})
        end
        hold off
    end
end
export_fig(sprintf('%s/SubjPerformance.pdf',figureLocation));

    
       
