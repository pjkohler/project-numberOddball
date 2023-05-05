close all;
clear all;

code_folder = '/Users/kohler/code';
addpath(genpath(sprintf('%s/git/export_fig',code_folder)),'-end');
addpath(genpath(sprintf('%s/git/mrC',code_folder)),'-end');

opt.data_type = 'RLS';
do_exp = 1;
data_location = sprintf('/Volumes/GoogleDrive/My Drive/numeroOddball/Experiment%.0d', do_exp);
figure_location = "/Users/kohler/Google Drive/WRITING/Articles/2019_KohlerNumerositySSVEP/figures/finished_figures";

folder_names=subfolders(sprintf('%s/20*',data_location),1);
for s = 1:length(folder_names)
    temp_folders = subfolders(folder_names{s},1);
    % remove folders that are used to store whatever data (corrected or uncorrected) 
    % that are not to be used
    temp_folders = temp_folders(~ismember(temp_folders, [folder_names{s},'/not_time_corrected']));
    temp_folders = temp_folders(~ismember(temp_folders, [folder_names{s},'/time_corrected']));
    mat_files = subfiles(sprintf('%s/ALL_Exp_MATL_HCN_128_Avg/RT*',temp_folders{end}),1);
    blockNum = 0;
    for m = 1:length(mat_files)
        tmp_data = load(mat_files{m});
        if ~isempty(tmp_data.CndTiming)
           blockNum = blockNum + 1;
           if s == 1 && blockNum == 1
               numTrials = size(tmp_data.TimeLine,1); % trials per block, assume same for all
               conditions = unique(cat(1,tmp_data.TimeLine.cndNmb));
           end
           trialIdx = (1:numTrials)+(blockNum-1)*numTrials;
           resp_data(trialIdx,1,s) = cat(1,tmp_data.TimeLine.cndNmb); % condition label
           resp_data(trialIdx,2,s) = cell2mat(cellfun(@(x) find(ismember({'Mis','Ra','La'},x)),{tmp_data.TimeLine.respString},'uni',false))-1; % response (0 = mis, 1 = Ra, 2 = La )
           resp_data(trialIdx,3,s) = cat(1,tmp_data.TimeLine.respTimeSec); % response time
           clear tmp_data;
        else
        end
    end
    IDs{s} = folder_names{s}(end-6:end);
    mis_idx(:,s) = resp_data(:,2,s)==0;
   % Change correct response mapping for control and experimental
   % condiitons
   if do_exp == 1
       corr_resp = [1, 2, 1, 2, 1, 2];
   elseif do_exp == 2
      corr_resp = [1, 2, 1, 1, 2, 1];
   else
       corr_resp = [2, 1, 2, 1, 2, 1, 2, 1];
   end
   for c = 1:length(conditions)
       curIdx = resp_data(:,1,s) == conditions(c);
       perc_mis(c,s) = length(find(mis_idx(curIdx,s)))./length(find(curIdx));
       ave_acc(c,s) = (sum(resp_data(~mis_idx(:,s) & curIdx,2,s) == corr_resp(c))+0.5)./(length(resp_data(~mis_idx(:,s) & curIdx,2,s) == corr_resp(c))+1);
   end
end

%% put into variables, calculate dprime, FA, TP
% look at data
mean(ave_acc,1)
mean(ave_acc,2)

fa_idx = find(corr_resp == 2); % False alarm
hr_fa = ave_acc;
hr_fa(fa_idx,:) = 1 - hr_fa(fa_idx,:); % FA = 1-accuracy
Zsc = norminv(hr_fa);
dPr = zeros(length(conditions)/2,length(folder_names));
bias = zeros(length(conditions)/2,length(folder_names));

hit_idx = find(corr_resp == 1); % hits
if do_exp == 1
    which_fa = [1, 2, 3];
elseif do_exp == 2
    which_fa = [1, 1, 2, 2];
else
    which_fa = [1, 2, 3, 4];
end

for c=1:length(hit_idx)
    disp([hit_idx(c), fa_idx(which_fa(c))])
    dPr(c,:) = Zsc(hit_idx(c),:) - Zsc(fa_idx(which_fa(c)),:);
    bias(c,:) = -.5 * (Zsc(hit_idx(c),:) + Zsc(fa_idx(which_fa(c)),:));
end

%% prepare plotting
lWidth = 1.5;
fSize = 12;
freq_labels = {'6 & 1Hz','3.75 & 0.75Hz','3 & 0.5Hz'};
avgdPr = mean(dPr,2);
avgBias = mean(bias,2);
%errdPr = std(dPr,[],2)./sqrt(length(folder_names)-1);
%errBias = std(bias,[],2)./sqrt(length(folder_names)-1);
errdPr = std(dPr,[],2);
errBias = std(bias,[],2);

% Individual subject plot
cBrewer = load('colorBrewer_new');
colors = cBrewer.rgb20;
gray_color = mean(colors(15:16,:),1);
colors = colors(~ismember(1:length(colors),[15,16]),:);
colors = [colors(1:2:end,:); colors(2:2:end,:)];
colors = colors(1:length(IDs),:);
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
text_opts = {'fontweight','normal','fontname','Helvetica','fontsize',fSize};
%% averages
if do_exp == 1
    x_vals = 1:3;
    x_lims = [0.5, 3.5];
    cond_labels = freq_labels;
    cond_order = 1:3;
elseif do_exp == 2
    x_vals = 1:4;
    x_lims = [0.5, 4.5];
    cond_labels = {'ref 6 odd 5','ref 6 odd 9','ref 8 odd 9','ref 8 odd 5'};
    cond_order = [3,4,2,1];
else
    x_vals = 1:4;
    x_lims = [0.5, 4.5];
    cond_labels = {'ref 5 odd 8','ref 6 odd 9','ref 8 odd 5','ref 9 odd 6'};
    cond_order = 1:4;
end

for z = 1:2
    subplot(2,1,z)
    if z == 1
        errorb(avgdPr(cond_order), errdPr(cond_order),'Color',colors(1,:), 'barwidth', 0.25);
        hold on
        pH = plot(x_vals, avgdPr(cond_order),'o','markersize',10,'LineWidth',lWidth,'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',[1,1,1]);
        y_lims = [-1,4];
        ylabel('dPrime');
    else
        errorb(avgBias(cond_order), errBias(cond_order),'Color',colors(1,:), 'barwidth', 0.25);
        hold on
        pH = plot(x_vals, avgBias(cond_order),'o','markersize',10,'LineWidth',lWidth,'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',[1,1,1]);
        plot(x_vals, avgBias(cond_order),'.','markersize',20,'Color','b');
        y_lims = [-1,1.5];
        ylabel('bias (c)');
    end
    uistack(pH, 'top');
    xlabel('freq pairs');
    set(gca,gcaOpts{:},'xtick',x_vals,'xticklabel',cond_labels,'xlim',x_lims,'ylim',y_lims,'box','off');
    hold off
end
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = 20;
figPos(3) = 20;
set(gcf,'pos',figPos);
%export_fig(sprintf('%s/avg_dPrime.pdf',figure_location),'-pdf','-transparent',gcf);

%% individual plots
figure;
xFigs = 2; % dPrime or Bias
yFigs = 4; %length(cond_order);% Freq pairs
titleStr = {'d''', 'bias (c)'};

for c=1:2 %dPrime or Bias
    for f=1:yFigs
        subplot(xFigs, yFigs, f+(c-1)*yFigs)
        if f > length(cond_order)
            set(gca,'visible','off');
            continue;
        else
        end 
        if c == 1
            dataToPlot = dPr;
            ave_data = [avgdPr(cond_order(f))-errdPr(cond_order(f)), avgdPr(cond_order(f))+errdPr(cond_order(f))];
            y_min = -2;
            y_max = 6;
            y_units = 2;
        else
            dataToPlot = bias;
            ave_data = [avgBias(cond_order(f))-errBias(cond_order(f)), avgBias(cond_order(f))+errBias(cond_order(f))];
            y_min = -2;
            y_max = 2;
            y_units = 1;
        end
        x_min = 0.5;
        x_max = 15.5;
        pa_h(2) = patch([x_min,x_min,x_max,x_max], [ave_data, fliplr(ave_data)],gray_color);
        set(pa_h(2),'linestyle','none');
        hold on;
        pa_h(1) = plot([x_min, x_max], ones(2,1).*mean(ave_data), 'w--', 'linewidth',lWidth);
        scatter(1:length(IDs),dataToPlot(cond_order(f),:),100,colors,'linewidth',lWidth);
        %plot([1:length(IDs)],dataToPlot,'.','markersize',20);
        if c==1
            ref_h(1) = plot([x_min, x_min], [y_min, y_max],'k-','linewidth',lWidth);
        else
            ref_h(1) = plot([x_min, x_min], [y_min, y_max],'k-','linewidth',lWidth);
            ref_h(2) = plot([x_min, x_max], zeros(1,2),'k-','linewidth', lWidth);
        end
        arrayfun(@(x) uistack(x,'bottom'), [ref_h, pa_h])
        if f == 1
            ylabel(['\it',titleStr{c}], text_opts{:})
        end
        
        if c == 1
            title(['\it', cond_labels{f}], text_opts{:})
        else
            if f == 1
                %xlabel('participants', text_opts{:});
                ax_pos = get(gca,'position');
                ax_pos(4) = ax_pos(4)+ax_pos(4)*.2;
                set(gca,'position', ax_pos);
            else
                cur_pos = get(gca,'position');
                cur_pos(4) = ax_pos(4);
                set(gca,'position', cur_pos);
            end
        end
        yMax = max(max(dataToPlot));
        set(gca,gcaOpts{:},'xtick',[], 'ytick',y_min:y_units:y_max, 'xticklabel', {''}, 'xlim',[x_min, x_max],'ylim',[y_min, y_max],'box','off','clipping','off');
        %xtickangle(90);
        hold off
    end
    
end
tightfig;
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = 10;
figPos(3) = 20;
set(gcf,'pos',figPos);
export_fig(sprintf('%s/exp%d_behavior.pdf',figure_location, do_exp),'-pdf','-transparent','-painters',gcf);

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