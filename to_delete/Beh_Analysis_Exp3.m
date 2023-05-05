close all;
clear all;
dataLocation = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/numeroOddball/Experiment3';
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
   % Change correct response mapping for control and experimental
   % condiitons
   cntrlcond = [1,3,5,7];
   for c=1:length(conditions)
       if any(c==cntrlcond)
           corrResp = 2;
       else
           corrResp = 1;
       end
       curIdx = respData(:,1,s) == conditions(c);
       TrialIdx(:,s,c) = respData(curIdx,2,s) == corrResp;
       respData(~misIdx(:,s) & curIdx,4,s) = respData(~misIdx(:,s) & curIdx,2,s) == corrResp; % Calculate accuracy
       percMis(c,s) = length(find(misIdx(curIdx,s)))./length(find(curIdx)); % Calculate percent missing
       aveAcc(c,s) = (sum(respData(~misIdx(:,s) & curIdx,2,s) == corrResp)+0.5)./(length(respData(~misIdx(:,s) & curIdx,2,s) == corrResp)+1); % average accuracy across trials
   end
end
%Put into variables, calculate dprime, FA, TP



Carriers = [5 6 8 9];
Oddballs = [8 9 5 6];
FAIdx = cntrlcond; % false alarm
HRmask = ~ismember(conditions,FAIdx); % hit rate

%Look at data
mean(aveAcc,1)
mean(aveAcc,2)




HrFa = aveAcc;
HrFa(FAIdx,:) = 1 - HrFa(FAIdx,:); % FA = 1-accuracy
Zsc = norminv(HrFa);
dPr = zeros(sum(HRmask),length(folderNames));
bias = zeros(sum(HRmask),length(folderNames));

for c=1:length(Carriers) % number of carriers
    car = FAIdx(c);
    odd = car+1;
    dPr(c,:) = Zsc(odd,:) - Zsc(car,:);
    bias(c,:) = -(Zsc(odd,:) + Zsc(car,:))/2;
end


% % Identify accurate subjects to analyze separately using 6v5 condition (the
% % most even condition)
% accurateS = dPr(1,:)>=median(dPr(1,:));
% inaccurateS = ~accurateS;
% filename= sprintf('%s/mediansplit.mat',dataLocation);
% save(filename,'accurateS','inaccurateS');
% 
% % Identify correct and incorrect trials by subject to indes RCA analysis
% % only with correct trials)
% filename = sprintf('%s/TrialIdx.mat',dataLocation);
% save(filename,'TrialIdx');



%For plotting
figureLocation = sprintf('%s/figures',dataLocation);
lWidth = 2;
condLabels = {'5v8','6v9','8v5','9v6'};
gcaOpts = {'box','off','fontname','Helvetica','linewidth',lWidth};
avgdPr = mean(dPr,2);
avgBias = mean(bias,2);
errdPr = std(dPr,[],2)./sqrt(length(folderNames)-1);
errBias = std(bias,[],2)./sqrt(length(folderNames)-1);


%Increasing Vs Decreasing
plot([1 2], avgdPr(1:2),'.','markersize',20,'Color','b');
hold on
errorbar([1 2],avgdPr(1:2),errdPr(1:2),'Color','b');
plot([3 4], avgdPr(3:4),'.','markersize',20,'Color','b');
errorbar([3 4],avgdPr(3:4),errdPr(3:4),'Color','b');
ylabel('dPrime');
xlabel('Stim Pairs');
yMin = min(min(dPr));
yMax = max(max(dPr));
set(gca,gcaOpts{:},'xtick',[1,2,3,4],'xticklabel',condLabels, 'xlim',[0.5,4.5],'ylim',[0,yMax+0.5],'box','off');
hold off
export_fig(sprintf('%s/Exp3_avg_dPrime_dir.pdf',figureLocation),'-pdf','-transparent',gcf);

%5-8 Vs 6-9
plot([1 2], avgdPr([1 3]),'.','markersize',20,'Color','b');
hold on
errorbar([1 2],avgdPr([1 3]),errdPr([1 3]),'Color','b');
plot([3 4], avgdPr([2 4]),'.','markersize',20,'Color','b');
errorbar([3 4],avgdPr([2 4]),errdPr([2 4]),'Color','b');
ylabel('dPrime');
xlabel('Stim Pairs');
yMin = min(min(dPr));
yMax = max(max(dPr));
set(gca,gcaOpts{:},'xtick',[1,2,3,4],'xticklabel',condLabels([1 3 2 4]), 'xlim',[0.5,4.5],'ylim',[0,yMax+0.5],'box','off');
hold off
export_fig(sprintf('%s/Exp3_avg_dPrime_num.pdf',figureLocation),'-pdf','-transparent',gcf);


% Bias Increasing Vs Decreasing
plot([1 2], avgBias(1:2),'.','markersize',20,'Color','b');
hold on
errorbar([1 2],avgBias(1:2),errBias(1:2),'Color','b');
plot([3 4], avgBias(3:4),'.','markersize',20,'Color','b');
errorbar([3 4],avgBias(3:4),errBias(3:4),'Color','b');
ylabel('Bias');
xlabel('Stim Pairs');
yMin = min(min(bias));
yMax = max(max(bias));
set(gca,gcaOpts{:},'xtick',[1,2,3,4],'xticklabel',condLabels, 'xlim',[0.5,4.5],'ylim',[yMin,yMax],'box','off');
hold off
export_fig(sprintf('%s/Exp3_avg_Bias_dir.pdf',figureLocation),'-pdf','-transparent',gcf);

% Bias 5-8 Vs 6-9
plot([1 2], avgBias([1 3]),'.','markersize',20,'Color','b');
hold on
errorbar([1 2],avgBias([1 3]),errBias([1 3]),'Color','b');
plot([3 4], avgBias([2 4]),'.','markersize',20,'Color','b');
errorbar([3 4],avgBias([2 4]),errBias([2 4]),'Color','b');
ylabel('Bias');
xlabel('Stim Pairs');
yMin = min(min(bias));
yMax = max(max(bias));
set(gca,gcaOpts{:},'xtick',[1,2,3,4],'xticklabel',condLabels([1 3 2 4]), 'xlim',[0.5,4.5],'ylim',[yMin,yMax],'box','off');
hold off
export_fig(sprintf('%s/Exp3_avg_Bias_num.pdf',figureLocation),'-pdf','-transparent',gcf);

%Individual subject plot
cBrewer = load('colorBrewer');
colors = cBrewer.rgb20(round(linspace(1,20,length(IDs))),:);
yFigs = 2;% dPrime or Bias
xFigs = 4;% Stim pairs
titleStr = {'5v8','6v9','8v5','9v6'};
titleStr2 = {'dPrime', 'Bias'};
for f=1:4 % Stim pairs
    for c=1:2 %dPrime or Bias
        subplot(yFigs,xFigs,xFigs*(c-1)+f)
        if c == 1
            yMin = min(min(dPr));
            yMax = max(max(dPr));
            dataToPlot = dPr(f,:);
        else
            dataToPlot = bias(f,:);
            yMin = min(min(bias));
            yMax = max(max(bias));
        end
        scatter(1:length(IDs),dataToPlot,60,colors,'filled');
        %plot([1:length(IDs)],dataToPlot,'.','markersize',20,'Color',colors);
        hold on
        if c==1
            hline = refline([0,errdPr(f)*2]); % 2* std error of dprime sample distribution
%             hline = refline([0,2]);
            
        else
            hline = refline([0,0]); %bias = 0
        end
        hline.Color = 'k';
        hline.LineWidth = 2;
        
        set(gca,gcaOpts{:},'xtick',[1:length(IDs)],'xticklabel',IDs, 'xlim',[0.5,length(IDs)+0.5],'ylim',[yMin-0.5,yMax+0.5],'box','off');
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
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = 14;
figPos(3) = 28;
set(gcf,'pos',figPos);
export_fig(sprintf('%s/Exp3_SubjPerformance.pdf',figureLocation),'-pdf','-transparent',gcf);

% 2-Way repeated measures ANOVA dPR
n = numel(dPr);
Y = reshape(dPr,n,1);
S = reshape(repmat(1:size(dPr,2),size(dPr,1),1),n,1);
F1 = reshape(repmat([1 1 2 2]',1,size(dPr,2)),n,1); % Number pair (carrier 6 vs carrier 8)
F2 = reshape(repmat([1 2 1 2]',1,size(dPr,2)),n,1); % Distance (1 vs 3)
FACTNAMES = {'carrier','distance'};

res_dP = rm_anova2(Y,S,F1,F2,FACTNAMES);

% 2-Way repeated measures ANOVA dPR
n = numel(bias);
Y = reshape(bias,n,1);
S = reshape(repmat(1:size(bias,2),size(bias,1),1),n,1);
F1 = reshape(repmat([1 1 2 2]',1,size(bias,2)),n,1); % Number pair (carrier 6 vs carrier 8)
F2 = reshape(repmat([1 2 1 2]',1,size(bias,2)),n,1); % Distance (1 vs 3)
FACTNAMES = {'carrier','distance'};

res_bias = rm_anova2(Y,S,F1,F2,FACTNAMES);






    
       
