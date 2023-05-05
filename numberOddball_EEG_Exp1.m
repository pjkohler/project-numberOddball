function numberOddball_EEG_Exp1(varargin)
    % Description:	analyze data from numberOddball Experiment 2
    % 
    % Syntax:	numberOddball_Exp2(<options>)
    % <options>
    %    
    %   folders_to_use    - vector of indices of subject folders to analyze
    %                       [analyze all available subject folders]
    %
    %   launch_rca       - if true RCA is run, if false RCA data is loaded
    %                   	[true]/false
    %
    %   trial_error      - compute trial-wise error true/[false]
    %
    %   plot_split       - plot RC run separately on odd and carrier
    %                       true/[false]
    %
    %   plot_axx         - plot axx waveforms passed through the rc
    %                       components [true]/false
    %
    %   force_source - force reading in of raw data when initial running
    %                       RCA [true]/false
    %
    %   axx_type         - string indicating the type of axx data to plot
    %                      ['ALL']/'NF1'
    %   use_projected   - if true use projected amplitudes instead of vector based stats
    %                     true/[false]
    %
    %   beh_split      - logical indicating whether to split the data on
    %                       behavior and only use correct trials true/[false]
    %   data_type      - string, ['RLS'] or 'DFT'
    
    %% ADD PATHS
    close all;
    if ~exist('code_folder','var')
        code_folder = '/Users/kohler/code';
        addpath(genpath(sprintf('%s/git/rcaBase',code_folder)),'-end');
        addpath(genpath(sprintf('%s/git/mrC',code_folder)),'-end');
        addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/misc',code_folder)),'-end');
        addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/figure',code_folder)),'-end');
        addpath(genpath('~/code/git/export_fig'));
    else
    end
    setenv('DYLD_LIBRARY_PATH','')
    
    %% PARSE ARGS
    opt	= Parse_Args(varargin,...
            'folders_to_use', [], ...
            'launch_rca', false, ...
            'trial_error', false, ...
            'plot_split', false, ...
            'plot_axx', true, ...
            'force_source', false, ...
            'axx_type', 'ALL', ...
            'use_projected', true, ...
            'beh_split', false, ...
            'data_type', 'RLS' ...
            );

    %% VARIABLES THAT CAN BE CHANGED
    do_exp = 1;
    top_folder = sprintf('/Volumes/MyGoogleDrive/ONGOING/numeroOddball/Experiment%.0d',do_exp);
    cond_names = {'dist 3','control'};

    %% IDENTIFY DATA LOCATION
    data_location = sprintf('/Volumes/MyGoogleDrive/WRITING/Articles/2019_KohlerNumerositySSVEP/figures/results/experiment%.0d/%s', do_exp, opt.data_type);
    fig_location = sprintf('/Volumes/MyGoogleDrive/WRITING/Articles/2019_KohlerNumerositySSVEP/figures/results/experiment%.0d', do_exp);
    if ~exist(data_location,'dir')
        mkdir(data_location);
    else
    end
    if ~exist(fig_location,'dir')
        mkdir(fig_location);
    else
    end

    folder_names=subfolders(sprintf('%s/*20*',top_folder),1);
    if ~isempty(opt.folders_to_use)
        folder_names = folder_names(opt.folders_to_use);
    else
    end
    % and string for saving the data
    saveStr = datestr(clock,26);
    saveStr(strfind(saveStr,'/')) ='';
    save_name = sprintf('%s/rca_data%s.mat',data_location, saveStr); % include the date as a string;

    %% SETUP INPUTS
    use_bins = 0; % use average bin, the zeroeth one
    use_freqs = 1:8; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
    use_conds = 1:6;
    use_trials = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    use_subs = [];   % subset of subjects to use for analysis (if set to false or empty, all trials will be used)
        
    compare_chan = 75; % comparison channel
    
    carrier_freq = [6.0, 3.75, 3] ;
    odd_freq = [1.0, 0.75, 0.5 ];
    axx_reps = carrier_freq./odd_freq;
    
    n_reg = 10; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    n_comp = 5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

    %% call the function --- USER runs this section (no editing necessary) ---
    % generate "filter" for comparison data
    w_comparison=zeros(128,1); w_comparison(compare_chan)=1; 

    if ~opt.launch_rca % if users says analysis is not to be launched
        temp_name = subfiles(sprintf('%s/rca_data*.mat', data_location),1);
        temp_name = temp_name{end}; % grab newest file and loose .mat suffix;
        load(temp_name);
    else
        for f = 1:length(folder_names)
            temp_folders = subfolders(folder_names{f},1);
            temp_folders = temp_folders(~ismember(temp_folders, [folder_names{f},'/not_time_corrected']));
            temp_folders = temp_folders(~ismember(temp_folders, [folder_names{f},'/time_corrected']));
            path_names{f} = sprintf('%s/Exp_TEXT_HCN_128_Avg',temp_folders{end});
        end
        if ~exist([save_name,'.mat'], 'file')
            save(save_name,'use_*') % create file if it does not exist
        else
            save(save_name,'use_*','-append') % otherwise append
        end
        warning('off','all')
        for c = 1:(length(use_conds)/2)
            cur_cond = [c-1,c]+c;
            % full RCA, use all harmonics 1F1, 2F1, 1F2 and 2F2
            if c==1
                rca_full(c,:) = rcaSweep(path_names, use_bins, use_freqs, cur_cond, use_trials, use_subs, n_reg, n_comp, opt.data_type, compare_chan, false, false); % TRUE
            else
                % since this is just a subset of the previous RCA, set force_source to false
                rca_full(c,:) = rcaSweep(path_names, use_bins, use_freqs, cur_cond, use_trials, use_subs, n_reg, n_comp, opt.data_type, compare_chan, false, false);
            end
            if c == 2
                for f = 1:size(rca_full,2)
                    rca_full(c,f) = flipSwapRCA(rca_full(c,f),1:5,2);
                end
            else
            end
            %elseif c == 1
            %    for f = 1:size(rca_full,2)
            %        rca_full(c,f) = flipSwapRCA(rca_full(c,f),1:5,2);
            %    end
            %end
            rcaH = grabCovFig(gcf);
            %export_fig(sprintf('%s/rca_full_cond%.0d&%.0d_cov.pdf', data_location, cur_cond(1), cur_cond(2)),'-pdf','-transparent',rcaH);
            close all;
            
            % create axx_rca_full
            full_w = [rca_full(c,1).W, w_comparison];
            axx_rca_full(c) = rcaWaveProject(path_names,full_w,cur_cond,opt.axx_type); % create struct

            % oddball RCA, first two frequencies 1F1 and 2F1 
            % since this is just a subset of the previous RCA, set force_source to false
            rca_odd(c, :) = rcaSweep(path_names, use_bins, use_freqs(1:end/2), cur_cond, use_trials, use_subs, n_reg, n_comp, opt.data_type, compare_chan, false, false);
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/rca_odd_cond%.0d&%.0d_cov.pdf', data_location, cur_cond(1), cur_cond(2)),'-pdf','-transparent',rcaH);
            close all;
            
            % carrier RCA, last two frequencies 1F2 and 2F2
            rca_carr(c, :) = rcaSweep(path_names, use_bins, use_freqs(end/2+1:end), cur_cond, use_trials, use_subs, n_reg, n_comp, opt.data_type, compare_chan, false, false);
            rcaH = grabCovFig(gcf);
            export_fig(sprintf('%s/rca_carr_cond%.0d&%.0d_cov.pdf', data_location, cur_cond(1), cur_cond(2)),'-pdf','-transparent',rcaH);
            oddW = [rca_odd(c).W, w_comparison];
            carrierW = [rca_carr(c).W, w_comparison];
            
            axx_rca_odd(c) = rcaWaveProject(path_names,oddW,cur_cond,opt.axx_type); % create struct
            axx_rca_carr(c) = rcaWaveProject(path_names,carrierW,cur_cond,opt.axx_type); % create struct
            close all;
        end
        save(save_name,'rca_full', 'rca_odd', 'rca_carr', 'axx_rca_*','-append')
        warning('on','all')
    end
    
    % beh_split label
    if opt.beh_split
        trial_type = 'cor_trials';
    else
        trial_type = 'all_trials';
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
                curRCA(c,f) = aggregateData(curRCA(c,f), keepConditions, errorType, opt.trial_error, doNR);
            end
            % axx_rca_full
            n_tps = size(curAxxRCA(c).rcaWave{1,1},1);
            n_rca = size(curAxxRCA(c).rcaWave{1,1},2);
            n_subs = size(curAxxRCA(c).rcaWave,2);
            n_cond = 2;
            curAxxRCA(c).ProjMat = reshape(cell2mat(curAxxRCA(c).rcaWave),n_tps,n_cond,n_rca,n_subs);
            curAxxRCA(c).avg = nanmean(curAxxRCA(c).ProjMat,4);
            curAxxRCA(c).std = nanstd(curAxxRCA(c).ProjMat,0,4);
            curAxxRCA(c).sem = curAxxRCA(c).std./sqrt(n_subs);
            for r = 1:size(curAxxRCA(c).ProjMat,3)
                [curAxxRCA(c).realT(:,r), curAxxRCA(c).realP(:,r), ~, curAxxRCA(c).corrT(:,r)] = ...
                    ttest_permute_sstats(squeeze(curAxxRCA(c).ProjMat(:,1,r,:) - curAxxRCA(c).ProjMat(:,2,r,:)), 5000, 'mass');
            end      
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
     % shut down parallel pool, which was used for fitting Naka-Rushton
    delete(gcp('nocreate'));
    clc;

    %% STATS
    for f = 1:max(use_freqs)
        for fPair = 1:3
            for rc_type = 1:2
                for c = 1:2
                    if rc_type == 1
                        [rca_data_real,rca_data_imag] = getRealImag(rca_full(fPair,f).rca_data);
                        rca_proj = squeeze(rca_full(fPair,f).subjects.proj_amp_signal);
                    else
                        if f > length(rca_odd)
                            [rca_data_real,rca_data_imag] = getRealImag(rca_carr(fPair,f-length(rca_odd)).rca_data);
                            rca_proj = squeeze(rca_split(fPair,f-length(rca_odd)).subjects.proj_amp_signal);
                        else
                            [rca_data_real,rca_data_imag] = getRealImag(rca_odd(fPair,f).rca_data);
                            rca_proj = squeeze(rca_split(fPair,f).subjects.proj_amp_signal);
                        end
                    end
                    
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
                        
                        % do the same thing for projected amplitudes
                        [between_tstatSig(f,r,fPair,rc_type), between_tstatP(f,r,fPair,rc_type), ci, sts] = ... 
                            ttest(squeeze(rca_proj(r,:,1)), squeeze(rca_proj(r,:,2)), 'alpha', 0.05, 'dim', 2, 'tail','right');
                        between_tstatVal(f,r,fPair,rc_type) = sts.tstat;
                    end
                end
            end
        end
    end
    
    for z = 1:3    
        [ rca_rel_expl(z), rca_var_expl(z) ,pca_var_expl ] = rcaExplained(rca_full(z,1),2);
    end

    
    %% NEW PLOTTING
    % note: plotting split data does not really work in the current version
    % could be built in fairly easily, but would require some work
    % rca_full seems cooler, simpler, and more informative anyway
    close all
    
    % significance colors
    p_colormap = jmaColors('pval');
    p_colormap(end,:) = [1 1 1];
    
    lWidth = 1;
    fSize = 12;
    text_opts = {'fontweight','normal','fontname','Helvetica','fontsize',fSize};
    gcaOpts = [{'tickdir','out','ticklength',[0.0200,0.0200],'box','off','linewidth',lWidth, 'clipping', 'off'}, text_opts];

    yFigs = 2;
    color_max = 0.2;

    if opt.plot_split == 0
        if opt.plot_axx 
            xFigs = 5;
        else
            xFigs = 3;
        end
    else
        if opt.plot_axx 
            xFigs = 10;
        else
            xFigs = 6;
        end
    end

    clear egiH;
    
    c_brewer = load('colorBrewer_new.mat');
    bold_colors = c_brewer.rgb20([9,3,5,13],:);
    weak_colors = c_brewer.rgb20([10,4,6,14],:);
    bar_width = 0.3;
    color_idx = [3,1];
    
    fig_label = {'A','B','C','D','E','F','G','H','I','J'};
    
    num_pairs = 3;
    
    figHeight = yFigs * num_pairs * 6;
    figWidth = xFigs * 6;
    
    figure;

    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = figHeight;
    figPos(3) = figWidth;
    set(gcf,'pos',figPos);
    
    pair_labels = {'6 & 1 Hz', '3.75 & 0.75 Hz', '3 & 0.5 Hz'};

    for f = 1:num_pairs % frequency pairs
        if f == 3
             axx_ticks = 0:500:2000;
        elseif f == 2
            axx_ticks = 0:200:1200;
        elseif  f == 1
            axx_ticks = 0:200:1000;
        else
        end
        % Axx xVals
        nTps = size(axx_rca_full(f).ProjMat,1);
        axx_xvals = linspace(0, 1000/carrier_freq(f) * axx_reps(f) , nTps+1);
        axx_xvals = axx_xvals(2:end);
        stim_onset = axx_xvals(floor(linspace(1,length(axx_xvals),axx_reps(f) + 1)));
        stim_onset = stim_onset(1:end-1);
        for r = 1:yFigs
            row_no = (r-1)*num_pairs+f;   
            if r < 6            
                if opt.plot_split == 0
                    curRCA = rca_full(f,:);
                    egiH(r,1) = subplot(yFigs*num_pairs,xFigs,3+(row_no-1)*xFigs);
                    hold on
                    [~, colorH] = mrC.plotOnEgi(curRCA(1).A(:,r),[-color_max, color_max], r==yFigs);
                    rcaType = 'full';
                else
                    curRCA = rca_split(f,:);
                    egiH(r,1) = subplot(yFigs*num_pairs,xFigs,3+(row_no-1)*xFigs);
                    hold on
                    [~, colorH(1)] = mrC.plotOnEgi(curRCA(1).A(:,r), [-color_max, color_max], r==yFigs);
                    hold off
                    egiH(r,2) = subplot(yFigs*num_pairs,xFigs,(xFigs-2)+(row_no-1)*xFigs);
                    hold on
                    [~, colorH(2)] = mrC.plotOnEgi(curRCA(5).A(:,r),[-color_max, color_max], r==yFigs);
                    rcaType = 'split';
                end
                if r == yFigs
                    for z = 1:length(colorH)
                        set(colorH(z),'units','centimeters','location','south', text_opts{:}, 'xtick', [-color_max:.1:color_max]);
                        color_pos = get(colorH(z),'position');
                        color_pos(1) = color_pos(1) - color_pos(3)*.4; 
                        color_pos(2) = color_pos(2) - color_pos(4)*4.3;
                        color_pos(3) = color_pos(3) * 1.8;
                        set(colorH(z),'position', color_pos, 'AxisLocation','out'); 
                        if f ~= num_pairs
                            set(colorH(z),'visible','off');
                        else
                        end
                    end
                else
                end
                arrayfun(@(x) uistack(x,'bottom'), egiH(r,:));
                hold off
            else
            end
            if opt.plot_axx
                % axx plot
                for z = 1:2
                    if opt.plot_split == 0
                        % full axx
                        subplot(yFigs*num_pairs, xFigs,[xFigs-1,xFigs]+(row_no-1)*xFigs);
                        cur_axx = axx_rca_full(f);
                        title_str = 'single cycle waveform';
                    else
                        if z == 1
                            subplot(yFigs*num_pairs,xFigs,[xFigs/2-1,xFigs/2]+(row_no-1)*xFigs);
                            cur_axx = axx_rca_odd(f);
                            title_str = 'single cycle waveform (odd)';
                        else
                            subplot(yFigs*num_pairs,xFigs,[xFigs-1,xFigs]+(row_no-1)*xFigs);
                            cur_axx = axx_rca_carr(f);
                            title_str = 'single cycle waveform (carrier)';
                        end
                    end
                    dataToPlot = squeeze(cur_axx.avg(:,:,r));
                    errorToPLot = squeeze(cur_axx.sem(:,:,r));
                    %axx_ymin = floor(min((min(cur_axx.avg(:,:,r)-cur_axx.sem(:,:,r)))/2))*2;
                    %axx_ymax = ceil(max((max(cur_axx.avg(:,:,r)+cur_axx.sem(:,:,r)))/2))*2;
                    %axx_ymax = max(abs([axx_ymin, axx_ymax]));
                    %axx_ymin = -axx_ymax;

                    if row_no > 3
                        axx_yunit = 4; axx_ymin = -8; axx_ymax = 8;

                    else
                        axx_yunit = 2; axx_ymin = -6; axx_ymax = 6;
                    end
                    axx_xmin = 0;
                    axx_xmax = max(axx_xvals);
                    sig_pos(1) = axx_ymax;
                    sig_pos(2) = axx_ymax-( axx_ymax-axx_ymin) *.1; 

                    hold on
                    for c=1:size(dataToPlot,2)
                        plot(axx_xvals,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(color_idx(c),:));
                        ErrorBars(axx_xvals',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(color_idx(c),:));
                    end
                    ylim([axx_ymin, axx_ymax]);
                    xlim([axx_xmin, axx_xmax]);
                    set(gca,gcaOpts{:},'xtick',axx_ticks,'xticklabel',cellfun(@(x) num2str(x),num2cell(axx_ticks),'uni',false),'ytick',axx_ymin:axx_yunit:axx_ymax);
                    yLine = repmat(get(gca,'YLim'),axx_reps(f),1)';
                    line(repmat(stim_onset,2,1),yLine,'Color','black');
                    if r == 1 && f == 1
                        title(title_str, text_opts{:}, 'fontangle','italic');
                    elseif r == yFigs && f == num_pairs
                        xlabel('time (ms)', text_opts{:}, 'fontangle','italic');
                    else
                    end

                    % plot corrected t-values
                    regionIdx = bwlabel(cur_axx.corrT(:,r));
                    for m=1:max(regionIdx)
                        tmp = regionprops(regionIdx == m,'centroid');
                        idx = round(tmp.Centroid(2));
                        hTxt = text(axx_xvals(idx),sig_pos(1),'*','fontsize',18,'fontname','Helvetica','horizontalalignment','center','verticalalignment','cap');
                    end

                    % plot uncorrected t-values
                    cur_p = repmat( cur_axx.realP(:,r)',20,1 );
                    h_img = image([min(axx_xvals),max(axx_xvals)],[sig_pos(1),sig_pos(2)], cur_p, 'CDataMapping', 'scaled','Parent',gca);
                    colormap( gca, p_colormap );   
                    c_mapmax = .05+2*.05/(size(p_colormap,1));
                    set( gca, 'CLim', [ 0 c_mapmax ] ); % set range for color scale
                    set(gca, gcaOpts{:});
                    uistack(h_img,'bottom')

                    s = get(gca, 'Position');
                    set(gca, 'Position', [s(1)+0.01, s(2), s(3)*0.9, s(4)]);
                    if opt.plot_split == 0
                        continue
                    else
                    end  
                end
                hold off;
            end
            for t = 1:2
                subplot(yFigs*num_pairs,xFigs,1+(row_no-1)*xFigs+(t-1));
                hold on
                freq_set = max(use_freqs)/2;
                cond_set = 2;
                curIdx = (use_freqs(1:end/2))+(t-1)*freq_set;
                bar_spacing = (0:bar_width:(bar_width*(cond_set-1))) - mean(0:bar_width:(bar_width*(cond_set-1)));
                x_vals = repmat((1:freq_set),cond_set,1) + repmat(bar_spacing,freq_set,1)';
                for c=1:cond_set
                    if opt.use_projected
                        amp_vals = arrayfun(@(x) curRCA(x).mean.amp_signal(:,:,r,c), curIdx);
                        amp_temp = arrayfun(@(x) nanmean(curRCA(x).subjects.proj_amp_signal(:,:,r,:,c), 4), curIdx);
                        df = arrayfun(@(x) sum(~isnan(curRCA(x).subjects.proj_amp_signal(:,:,r,:,c))), curIdx);
                        if sum(amp_vals - amp_temp) > 1e-10
                            error('average based on projected values is different');
                        else
                            amp_vals = amp_temp;
                        end
                        err_ub = arrayfun(@(x) nanstd(curRCA(x).subjects.proj_amp_signal(:,:,r,:,c), 0, 4), curIdx)./sqrt(df);
                        err_lb = arrayfun(@(x) nanstd(curRCA(x).subjects.proj_amp_signal(:,:,r,:,c), 0, 4), curIdx)./sqrt(df);
                        within_sig = arrayfun(@(x) curRCA(x).stats.t_p(:,:,r,c)<0.05, curIdx);
                        
                        if c == 1
                            % DO LME
                            anova_c1 = cell2mat(arrayfun(@(x) squeeze(curRCA(x).subjects.proj_amp_signal(:,:,r,:,1))', curIdx, 'uni', false));
                            anova_c2 = cell2mat(arrayfun(@(x) squeeze(curRCA(x).subjects.proj_amp_signal(:,:,r,:,2))', curIdx, 'uni', false));
                            anova_ready = [anova_c1, anova_c2];
                            anova_names = {'dist3','dist0'};
                            cond_lbl = arrayfun(@(x) anova_names(x), [ones(size(anova_c1)), 2*ones(size(anova_c2))]);
                            if t == 1
                                if row_no > 3
                                    y_max = 3;
                                else
                                    y_max = 1;
                                end
                                title_str = 'oddball';
                                harm_labels = {'1F1','2F1','3F1','4F1'};
                                r_table = sprintf('%s/%s_oddball_rc%d_freq%d_%s_projected_%s.csv', fig_location, opt.data_type, r, f, rcaType, trial_type);
                            else
                                if row_no > 3
                                    y_max = 2;
                                else
                                    y_max = 4;
                                end
                                title_str = 'carrier';
                                harm_labels = {'1F2','2F2','3F2','4F2'};
                                r_table = sprintf('%s/%s_carrier_rc%d_freq%d_%s_projected_%s.csv', fig_location, opt.data_type, r, f, rcaType, trial_type);
                            end
                            harm_lbl = arrayfun(@(x) harm_labels(x), repmat(1:4, length(cond_lbl)/8, 1));
                            harm_lbl = repmat(harm_lbl(:),2,1);
                            subj_lbl = arrayfun(@(x) num2str(x, 's%.02d'), repmat(1:length(anova_c1)/4,1,8), 'uni', false);
                            tbl = table(cond_lbl', subj_lbl', harm_lbl, anova_ready','VariableNames', {'condition','subject','harmonic','data'});
                            tbl.subject = categorical(tbl.subject);
                            tbl.condition = categorical(tbl.condition);
                            tbl.harmonic = categorical(tbl.harmonic);
                            writetable(tbl, r_table, 'Delimiter',',','QuoteStrings',true);
                            anova_file = strrep(r_table,"trials","results");
                            if exist(anova_file,'file')
                                anova_out = readtable(anova_file,'Delimiter',',','ReadRowNames',true);
                            else
                            end 
                        else
                        end
                        
                    else
                        amp_vals = arrayfun(@(x) curRCA(x).mean.amp_signal(:,:,r,c), curIdx);
                        err_ub = arrayfun(@(x) curRCA(x).stats.amp_up_err(:,:,r,c), curIdx);
                        err_lb = arrayfun(@(x) curRCA(x).stats.amp_lo_err(:,:,r,c), curIdx);
                        within_sig = arrayfun(@(x) curRCA(x).stats.t2_p(:,:,r,c)<0.05, curIdx);
                    end
                    errorbar(x_vals(c,:), amp_vals, err_lb, err_ub, '.k', 'LineWidth',lWidth, 'marker','none');
                    amp_h(c) = bar(x_vals(c,:), amp_vals,'BarWidth',bar_width,'edgecolor','none','facecolor',bold_colors(color_idx(c),:));
                    for z = 1:length(within_sig)
                        if within_sig(z)
                            patch_x = [ x_vals(c,z)-bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)-bar_width/2];
                            patch_y = [0, 0, y_max, y_max];
                            pa_h = patch(patch_x, patch_y, bold_colors(color_idx(c),:),'edgecolor','none','facealpha',.25); 
                            uistack(pa_h,'bottom');
                        else
                        end
                    end
                end
                % do between group tests
                if opt.use_projected
                    cur_p = squeeze(between_tstatP(curIdx,r,f, opt.plot_split+1));
                    between_idx = (cur_p < 0.05);
                    between_idx = between_idx + (cur_p < 0.005);
                else
                    cur_p = squeeze(between_t2p(curIdx, r, f, opt.plot_split+1));
                    cur_sign = squeeze(between_sign(curIdx, r,f, opt.plot_split+1) == 1);
                    between_idx = (cur_p < 0.05) & (cur_sign == 1);
                    between_idx = between_idx + ((cur_p < 0.005) & (cur_sign == 1));
                end
               
                %if any(between_idx)
                %    arrayfun(@(x) text(x,y_max*.95,'*', text_opts{:}, 'fontsize',20,'HorizontalAlignment','center'), find(between_idx > 0),'uni',false);
                %    arrayfun(@(x) text(x,y_max*.85,'*', text_opts{:}, 'fontsize',20,'HorizontalAlignment','center'), find(between_idx == 2),'uni',false);
                %else
                %end
                p_opts = {'fontsize', 10, 'Interpreter', 'tex', 'HorizontalAlignment','left'};
                if exist('anova_out','var');
                    if table2array(anova_out({'cond'},{'Pr__F_'})) < 0.05
                        main_p = table2array(anova_out({'cond'},{'Pr__F_'}));
                        if main_p < 0.0001
                            text(1,y_max*.95, 'N_{main} {\itp} < 0.0001', text_opts{:}, p_opts{:});
                        elseif main_p < 0.001
                            text(1,y_max*.95, 'N_{main} {\itp} < 0.001', text_opts{:}, p_opts{:});
                        else
                            text(1,y_max*.95, strcat("N_{main} {\itp} = ", sprintf("%.3f", main_p)), text_opts{:}, p_opts{:});
                        end
                    else
                    end
                    if table2array(anova_out({'cond:harm'},{'Pr__F_'})) < 0.05
                        inter_p = table2array(anova_out({'cond:harm'},{'Pr__F_'}));
                        if inter_p < 0.0001
                            text(1,y_max*.8, "N\timesH {\itp} < 0.0001", text_opts{:}, p_opts{:});
                        elseif inter_p < 0.001
                            text(1,y_max*.8, "N\timesH {\itp} < 0.001", text_opts{:}, p_opts{:});
                        else
                            text(1,y_max*.8, strcat("N\timesH {\itp} = ", sprintf("%.3f", inter_p)), text_opts{:}, p_opts{:});
                        end
                    else
                    end
                else
                end
                if y_max >= 2
                    y_unit = 1;
                else
                    y_unit = .2;
                end
                ylim([0,y_max]);
                xlim([.5,freq_set+0.5]);
                if r == yFigs  && t == 1 && f == num_pairs
                    ylabel('amplitude (\muVolts)', text_opts{:}, 'fontangle','italic');
                    xlabel('harmonics', text_opts{:}, text_opts{:}, 'fontangle','italic');
                else
                end
                if r==1 && f == 1
                    title(title_str, text_opts{:}, 'fontangle','italic');
                elseif r==yFigs && f == num_pairs
                    if t==2
                        lH = legend(amp_h, cond_names, text_opts{:});
                        legend boxoff
                        l_pos = get(lH,'position');
                        l_pos(1) = l_pos(1) - l_pos(3)*.5;
                        l_pos(2) = l_pos(2) - l_pos(4)*3.5;
                        set(lH,'position',l_pos);
                    else
                    end
                else
                end
                if t == 1
                    text(-1.5, y_max, fig_label{row_no}, text_opts{:}, 'fontsize', 20, 'fontangle','italic');
                    textH = text(-1.5, y_max*.35, pair_labels{f}, text_opts{:}, 'fontangle','italic',...
                            'HorizontalAlignment','center','verticalalignment','top');
                    set(textH,'Rotation',90);
                else
                end
                set(gca, gcaOpts{:}, 'xtick', [1,2,3,4], 'ytick', 0:y_unit:y_max);
                set(gca, 'xticklabel', harm_labels);
                hold off
            end
        end
        drawnow;
        for r = 1:yFigs
            for z = 1:size(egiH,2)
                if opt.plot_split
                    addVal = 0.2;
                    shiftLeft = 0.02;
                else
                    addVal = 0.4;
                    shiftLeft = 0.01;
                end
                newPos = get(egiH(r,z),'position');
                newPos(1) = newPos(1)-(newPos(3)*addVal/2) - shiftLeft;
                newPos(2) = newPos(2)-(newPos(4)*addVal/2)+(newPos(4)*addVal/30)*(r);
                newPos(3:4) = newPos(3:4)*(1+addVal);
                set(egiH(r,z),'position',newPos);
            end
        end
    end
    
    if opt.use_projected    
        fig_name = sprintf('%s/%s_all_%s_projected_%s', fig_location, opt.data_type, rcaType, trial_type);
    else
        fig_name = sprintf('%s/%s_all_%s_%s',fig_location, opt.data_type, rcaType, trial_type);
    end
    export_fig(sprintf('%s.png', fig_name),'-png','-opengl','-m5','-transparent',gcf);
end
