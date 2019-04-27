function numberOddball_Exp2(varargin)
    % Description:	analyze data from numberOddball Experiment 2
    % 
    % Syntax:	numberOddball_Exp2(<options>)
    % <options>
    %    
    %   folders_to_use    - vector of indices of subject folders to analyze
    %                       [analyze all available subject folders]
    %
    %   launch_rca       - if true RCA is run, if false RCA data is loaded
    %                   	true/[false]
    %
    %   trial_error      - compute trial-wise error true/[false]
    %
    %   plot_split       - plot RC run separately on odd and carrier
    %                       true/[false]
    %
    %   force_source - force reading in of raw data when initial running
    %                       RCA true/[false]
    %
    %   axx_type         - string indicating the type of axx data to plot
    %                      ['ALL']/'NF1'
    %
    %   train_data       - integer indicating whether to train RCA on all data 
    %                     or only on conditions with carrier = 6 or carrier = 8
    %                     [0](all)/1(6)/2(8)
    %
    %   use_projected   - if true use projected amplitudes instead of vector based stats
    %                     true/[false]
    %
    %   beh_split      - logical indicating whether to split the data on
    %                       behavior and only use correct trials true/[false]
    %   data_type      - string, ['RLS'] or 'DFT'
    
    %% ADD PATHS
    close all;
    if ~exist('codeFolder','var')
        code_folder = '/Users/kohler/code';
        addpath(genpath(sprintf('%s/git/rcaBase',code_folder)));
        addpath(genpath(sprintf('%s/git/mrC',code_folder)));
        addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',code_folder)));
        addpath(genpath(sprintf('%s/git/export_fig',code_folder)));
    else
    end
    setenv('DYLD_LIBRARY_PATH','')
    
    %% PARSE ARGS
    opt	= ParseArgs(varargin,...
            'folders_to_use', [], ...
            'launch_rca', true, ...
            'trial_error', false, ...
            'plot_split', false, ...
            'force_source', false, ...
            'axx_type', 'ALL', ...
            'train_data', 0, ...
            'use_projected', true, ...
            'beh_split', false, ...
            'data_type', 'RLS' ...
            );

    %% VARIABLES THAT CAN BE CHANGED
    do_exp = 2;
    top_folder = sprintf('/Volumes/Denali_DATA1/kohler/EEG_EXP/DATA/numeroOddball/Experiment%.0d',do_exp);

    %% IDENTIFY DATA LOCATION
    data_location = sprintf('/Users/kohler/Dropbox/WRITING/Articles/2019_KohlerNumerositySSVEP/figures/results/experiment%.0d/%s', do_exp, opt.data_type);
    fig_location = sprintf('/Users/kohler/Dropbox/WRITING/Articles/2019_KohlerNumerositySSVEP/figures/results/experiment%.0d', do_exp);
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
    carrier_freq = 3.0;
    axx_reps = 6;
    
    compare_chan = 75; % comparison channel
    
    % think how to compute the appropriate pairs to do the
    % ttests and how they are stored
    if do_exp == 2
        color_idx = [3,1,2];
        if opt.train_data == 0
            use_conds = [6 5 4, 1 2 3]; %6v9,6v6,6v5; 8v5,8v8,8v9
            %Reorder to make plotting easier
            carriers = [6, 8];
            contrast_order = [[1 2]; [3 2]; [1 3]; [4 5]; [6 5]; [4 6]];
        elseif opt.train_data == 1
            use_conds = [6 5 4];
            carriers = 6;
            train_stim = 'train6';
            contrast_order = [[1 2];[3 2];[1 3]];
        elseif opt.train_data == 2
            use_conds = [1 2 3];
            carriers = 8;
            train_stim = 'train8';
            contrast_order = [[1 2];[3 2];[1 3]];
        end
        cond_names = {'dist 3','dist 0','dist 1'};
        sig_labels = {'3','1','N'};
    elseif do_exp == 3
        color_idx = [3,1];
        if opt.train_data == 0
            use_conds = [2 1 4 3 6 5 8 7]; %5v8, 5v5; 6v9, 6v6; 8v5, 8v8; 9v6, 9v9
            carriers = [5 6 8 9];
            contrast_order = [[1 2];[3 4];[5 6];[7 8]];
        elseif opt.train_data == 1
            use_conds = [2 1 6 5];
            carriers = [5 8];
            train_stim = 'train5_8';
            contrast_order = [[1 2];[3 4]];
        elseif opt.train_data == 2
            use_conds = [4 3 8 7];
            carriers = [6 9];
            train_stim = 'train6_9';
            contrast_order = [[1 2];[3 4]];
        end
        cond_names = {'dist 3','dist 0'};
    end

    use_trials = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    n_reg = 10; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    n_comp = 5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    rc_plotstyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

    %% LAUNCH RCA OR LOAD DATA

    if ~opt.launch_rca % if users says analysis is not to be launched
        temp_name = subfiles(sprintf('%s/rca_data*.mat', data_location),1);
        temp_name = temp_name{end}; % grab newest file and loose .mat suffix;
        load(temp_name);
    else
        % generate "filter" for comparison data
        w_comparison = zeros(128,1); w_comparison(compare_chan) = 1; 
        % find text files
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
        % full RCA, use all harmonics 1F1, 2F1, 1F2 and 2F2
        rca_full = rcaSweep(path_names, use_bins, use_freqs, use_conds, use_trials, n_reg, n_comp, opt.data_type, compare_chan, [], rc_plotstyle, opt.force_source); % TRUE
        
        if do_exp == 3
            for f = 1:numel(rca_full)
                rca_full(1,f) = flipSwapRCA(rca_full(1,f),1:5,2);
            end
        else
        end
        
        rca_h = grabCovFig(gcf);
        export_fig(sprintf('%s/rca_full_cond%s_cov.pdf', data_location, sprintf('%.0d',use_conds)),'-pdf','-transparent',rca_h);

        close all;
        % create axxRCA
        full_w = [rca_full(1).W, w_comparison];
        axx_rca_full = axxRCAmake(path_names, full_w, use_conds, opt.axx_type); % create struct
        
        % oddball RCA, first two frequencies 1F1 and 2F1
        % since this is just a subset of the previous RCA, set force_source to false
        rca_odd = rcaSweep(path_names, use_bins, use_freqs(1:end/2), use_conds, use_trials, n_reg, n_comp, opt.data_type, compare_chan,[], rc_plotstyle, false);
        rca_h = grabCovFig(gcf);
        export_fig(sprintf('%s/rca_odd_cond%s_cov.pdf', data_location, sprintf('%.0d',use_conds)),'-pdf','-transparent',rca_h);
               
        close all;
        odd_w = [rca_odd(1).W, w_comparison];

        % carrier RCA, last two frequencies 1F2 and 2F2
        rca_carr = rcaSweep(path_names, use_bins, use_freqs(end/2+1:end), use_conds, use_trials, n_reg, n_comp, opt.data_type, compare_chan,[], rc_plotstyle, false);
        rca_h = grabCovFig(gcf);
        export_fig(sprintf('%s/rca_carr_cond%s_cov.pdf', data_location, sprintf('%.0d',use_conds)),'-pdf','-transparent',rca_h);
        
        carr_w = [rca_carr(1).W, w_comparison];
        
        axx_rca_odd = axxRCAmake(path_names, odd_w, use_conds, opt.axx_type); % create struct
        axx_rca_carr = axxRCAmake(path_names, carr_w, use_conds, opt.axx_type); % create struct

        close all;
        save(save_name,'rca_full', 'rca_odd', 'rca_carr', 'axx_rca_*', '-append')
    end
    
    %% NOW SPLIT BASED ON BEHAVIOR
    if opt.beh_split 
        % load trial idx, and convert it to a more manageable variable
        load(sprintf('%s/trialIdx.mat',data_location));
        if ~iscell(TrialIdx)
            TrialIdx = permute(TrialIdx,[1,3,2]);
            tSize = size(TrialIdx);
            trialIdx = cell(size(rca_full(1).data));
            if any(size(trialIdx) ~= tSize(2:end))
                error('mismatch between index and data');
            else
            end
            for q = 1 : prod(tSize(2:end)) 
                trialIdx{q} = TrialIdx(:,q);
            end
            clear TrialIdx;
        else
        end
    
        [rca_full, rca_full_incorrect] = splitRCA(rca_full, trialIdx);
        [rca_odd, rca_odd_incorrect] = splitRCA(rca_odd, trialIdx);
        [rca_carr, rca_carr_incorrect] = splitRCA(rca_carr, trialIdx);
    end
    
    % beh_split label
    if opt.beh_split
        trial_type = 'cor_trials';
    else
        trial_type = 'all_trials';
    end
    
    %% COMPUTE VALUES FOR PLOTTING

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
                % add the aggregate data and stats to the struct
                curRCA(c,f) = aggregateData(curRCA(c,f),keepConditions,errorType,opt.trial_error,doNR);
            end
            % AxxRCA
            n_tps = size(curAxxRCA(c).Projected{1,1},1);
            n_rca = size(curAxxRCA(c).Projected{1,1},2);
            n_subs = size(curAxxRCA(c).Projected,2);
            n_cond = length(use_conds);
            curAxxRCA(c).ProjMat = reshape(cell2mat(curAxxRCA(c).Projected),n_tps,n_cond,n_rca,n_subs);
            curAxxRCA(c).avg = nanmean(curAxxRCA(c).ProjMat,4);
            curAxxRCA(c).std = nanstd(curAxxRCA(c).ProjMat,0,4);
            curAxxRCA(c).sem = curAxxRCA(c).std./sqrt(n_subs);
            for r = 1:size(curAxxRCA(c).ProjMat,3)
                for t=1:size(contrast_order,1) % number of tTests
                    [curAxxRCA(c).realT(:,r,t), curAxxRCA(c).realP(:,r,t), ~, curAxxRCA(c).corrT(:,r,t)] = ...
                        ttest_permute_sstats(squeeze(curAxxRCA(c).ProjMat(:,contrast_order(t,1),r,:) - curAxxRCA(c).ProjMat(:,contrast_order(t,2),r,:)), 5000, 'mass');
                end
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
    % get the data and run the between group tests
    for f = 1:max(use_freqs)
        for rc_type = 1:2
            if rc_type == 1
                [rca_data_real,rca_data_imag] = getRealImag(rca_full(f).data);
                [comp_data_real,comp_data_imag] = getRealImag(rca_full(f).comparisonData);
                rca_proj = squeeze(rca_full(f).subjects.proj_amp_signal);
            else
                if f > length(rca_odd)
                    [rca_data_real,rca_data_imag] = getRealImag(rca_carr(f-length(rca_odd)).data);
                    [comp_data_real,comp_data_imag] = getRealImag(rca_carr(f-length(rca_odd)).comparisonData);
                    rca_proj = squeeze(rca_split(f-length(rca_odd)).subjects.proj_amp_signal);
                else
                    [rca_data_real,rca_data_imag] = getRealImag(rca_odd(f).data);
                    [comp_data_real,comp_data_imag] = getRealImag(rca_odd(f).comparisonData);
                    rca_proj = squeeze(rca_split(f).subjects.proj_amp_signal);
                end
            end
            % concatenate comparison and rca data
            rca_data_real = cellfun(@(x,y) cat(2,x,y), rca_data_real, comp_data_real, 'uni', false);
            rca_data_imag = cellfun(@(x,y) cat(2,x,y), rca_data_imag, comp_data_imag, 'uni', false);
            
            for r = 1:6
                x_data = cell2mat(cellfun(@(x) squeeze(nanmean(x(1,r,:),2)), rca_data_real, 'uni', false));
                y_data = cell2mat(cellfun(@(x) squeeze(nanmean(x(1,r,:),2)), rca_data_imag, 'uni', false));
                xy_data = cat(2,permute(x_data,[2,3,1]),permute(y_data,[2,3,1]));
                % compute vector mean amplitude, to get sign
                temp_amp = sqrt(nanmean(x_data,2).^2+nanmean(y_data,2).^2);
                
                for t=1:size(contrast_order,1) % number of tTests
                    results = tSquaredFourierCoefs(squeeze(xy_data(:,:,contrast_order(t,:))));
                    between_t2sig(f, r, rc_type, t) = results.H;
                    between_t2p(f, r, rc_type, t) = results.pVal;
                    between_t2Stat(f, r, rc_type, t) = results.tSqrd;
                    between_sign(f, r, rc_type, t) = sign(temp_amp(contrast_order(t,1))-temp_amp(contrast_order(t,2)));
                    % do the same thing for projected amplitudes
                    [between_tstatSig(f, r, rc_type, t), between_tstatP(f, r, rc_type, t), ci, sts] = ... 
                        ttest(squeeze(rca_proj(r,:,contrast_order(t,1))), squeeze(rca_proj(r,:,contrast_order(t,2))), 'alpha', 0.05, 'dim', 2, 'tail','right');
                    between_tstatVal(f, r, rc_type, t) = sts.tstat;
                    
                end
            end
        end
    end
    
    
    [ rca_rel_expl, rca_var_expl ,pca_var_expl ] = rcaExplained(rca_full(1),2);

    %% NEW PLOTTING
    close all
    axx_ticks = [0 500 1000 1500 2000];
    
    % significance colors
    p_colormap = jmaColors('pval');
    p_colormap(end,:) = [1 1 1];
    
    lWidth = 1;
    fSize = 10;
    gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

    yFigs = 2;

    if opt.plot_split == 0    
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
    bar_width = 0.3;
    
    for f = 1:length(carriers) % carriers (i.e. 6 & 8)
        % Axx xVals
        nTps = size(axx_rca_full(1).ProjMat,1);
        axx_xvals = linspace(0, 1000/carrier_freq * axx_reps , nTps+1);
        axx_xvals = axx_xvals(2:end);
        stim_onset = axx_xvals(floor(linspace(1,length(axx_xvals),axx_reps + 1)));
        stim_onset = stim_onset(1:end-1);
        if opt.plot_split
            cur_order = 1:length(cond_names);
        else
            cur_order = (1:length(cond_names))+(f-1)*length(cond_names);
        end
        figure;
        for r = 1:yFigs 
            if r < 6            
                if opt.plot_split == 0
                    curRCA = rca_full;
                    egiH(r,1) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs);
                    hold on
                    rcaColorBar = [min(curRCA(1).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    [~, colorH] = mrC.plotOnEgi(curRCA(1).A(:,r), rcaColorBar, r==yFigs);
                    rcaType = 'full';
                else
                    curRCA = rca_split;
                    egiH(r,1) = subplot(yFigs,xFigs,3+(r-1)*xFigs); %New Candidate
                    hold on
                    rcaColorBar = [min(curRCA(1).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    [~, colorH(1)] = mrC.plotOnEgi(curRCA(1).A(:,r),rcaColorBar, r==yFigs);
                    hold off
                    egiH(r,2) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs); %I think this one is correct
                    hold on
                    rcaColorBar = [min(curRCA(5).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme]; 
                    [~, colorH(2)] = mrC.plotOnEgi(curRCA(5).A(:,r),rcaColorBar, r==yFigs);
                    rcaType = 'split';
                end
                if r == yFigs
                    for z = 1:length(colorH)
                        set(colorH(z),'units','centimeters','location','north','fontsize',12);
                        color_pos = get(colorH(z),'position');
                        color_pos(1) = color_pos(1) + color_pos(3)*2.1; 
                        color_pos(2) = color_pos(2) - color_pos(4)*10.2;
                        color_pos(4) = color_pos(4) * 0.7;
                        color_pos(3) = color_pos(3) * 1.5;
                        set(colorH(z),'position', color_pos, 'AxisLocation','out'); 
                    end
                else
                end
                hold off
            else
            end

            % axx plot
            for z = 1:2
                if opt.plot_split == 0
                    % full axx
                    subplot(yFigs,xFigs,xFigs+(r-1)*xFigs);
                    subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs, xFigs+(r-1)*xFigs]);
                    cur_axx = axx_rca_full;
                    title_str = 'single cycle waveform (all)';
                else
                    if z == 1
                        subplot(yFigs,xFigs,[1+(r-1)*xFigs, 2+(r-1)*xFigs]);
                        cur_axx = axx_rca_odd;
                        title_str = 'single cycle waveform (odd)';
                    else
                        subplot(yFigs,xFigs,[(xFigs-1)+(r-1)*xFigs,xFigs+(r-1)*xFigs]);
                        cur_axx = axx_rca_carr;
                        title_str = 'single cycle waveform (carrier)';
                    end
                end
                dataToPlot = squeeze(cur_axx.avg(:,cur_order,r));
                errorToPLot = squeeze(cur_axx.sem(:,cur_order,r));
                axx_ymin = floor(min((min(cur_axx.avg(:,cur_order,r)-cur_axx.sem(:,cur_order,r)))/2))*2;
                axx_ymax = ceil(max((max(cur_axx.avg(:,cur_order,r)+cur_axx.sem(:,cur_order,r)))/2))*2;
                sig_pos(1) = axx_ymax;
                sig_pos(2) = axx_ymax-( axx_ymax-axx_ymin) *.1; 
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
                yLine = repmat(get(gca,'YLim'),axx_reps,1)';
                line(repmat(stim_onset,2,1),yLine,'Color','black');
                if r == 1
                    title(title_str,'fontweight','normal');
                elseif r == 6
                    xlabel('time (ms)');
                else
                end
                
                % plot corrected t-values
                num_tests = length(contrast_order)/length(carriers);
                test_idx = (1:num_tests)+(f-1)*num_tests;
                regionIdx = bwlabel(cur_axx.corrT(:,r, test_idx(1)));
                for m=1:max(regionIdx)
                    tmp = regionprops(regionIdx == m,'centroid');
                    idx = round(tmp.Centroid(2));
                    hTxt = text(axx_xvals(idx),sig_pos(1),'*','fontsize',18,'fontname','Helvetica','horizontalalignment','center','verticalalignment','cap');
                end

                % plot uncorrected t-values
                cur_p = repmat( cur_axx.realP(:,r, test_idx(1))',20,1 );
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
            hold off
            for t = 1:2
                subplot(yFigs,xFigs,(xFigs-4)+(r-1)*xFigs+(t-1));
                hold on
                freq_set = max(use_freqs)/2;
                cond_set = length(cur_order);
                curIdx = (use_freqs(1:end/2))+(t-1)*freq_set;
                bar_spacing = (0:bar_width:(bar_width*(cond_set-1))) - mean(0:bar_width:(bar_width*(cond_set-1)));
                x_vals = repmat((1:freq_set),cond_set,1) + repmat(bar_spacing,freq_set,1)';
                for c=1:length(cur_order)
                    if opt.use_projected
                        amp_vals = arrayfun(@(x) curRCA(x).mean.amp_signal(1, 1, r, cur_order(c)), curIdx);
                        amp_temp = arrayfun(@(x) squeeze(nanmean(curRCA(x).subjects.proj_amp_signal(1, 1, r,:, cur_order(c)), 4)), curIdx);
                        df = arrayfun(@(x) sum(~isnan(curRCA(x).subjects.proj_amp_signal(1, 1, r,:, cur_order(c)))), curIdx);
                        if sum(amp_vals - amp_temp) > 1e-10
                            error('average based on projected values is different');
                        else
                            amp_vals = amp_temp;
                        end
                        err_ub = arrayfun(@(x) nanstd(curRCA(x).subjects.proj_amp_signal(1, 1, r, :, cur_order(c)), 0, 4), curIdx)./sqrt(df);
                        err_lb = arrayfun(@(x) nanstd(curRCA(x).subjects.proj_amp_signal(1, 1, r, :, cur_order(c)), 0, 4), curIdx)./sqrt(df);
                        within_sig = arrayfun(@(x) curRCA(x).stats.t_p(1, 1, r, cur_order(c))<0.05, curIdx);
                        
                        if c == 1
                            % DO LME
                            anova_ready = cell2mat(arrayfun(@(x) curRCA(x).subjects.proj_amp_signal(1, 1, r,:,cur_order), curIdx, 'uni', false));
                            anova_ready = anova_ready(:);
                            if do_exp == 2 
                                anova_names = {'dist3','dist0','dist1'}; 
                            else
                                anova_names = {'dist3','dist0'};
                            end
                            cond_idx = repmat(1:length(anova_names), size(anova_ready,1)/length(anova_names),1);
                            cond_lbl = arrayfun(@(x) anova_names(x), cond_idx(:));
                            if t == 1
                                title_str = 'oddball';
                                harm_names = {'1F1','2F1','3F1','4F1'};
                                %r_table = sprintf('%s/%s_oddball_rc%d_freq%d_%s_projected_%s.csv', fig_location, opt.data_type, r, f, rcaType, trial_type);
                                r_table = sprintf('%s/%s_oddball_rc%d_carr%d_%s_projected_%s', fig_location, opt.data_type, r, carriers(f), rcaType, trial_type);  
                            else
                                title_str = 'carrier';
                                harm_names = {'1F2','2F2','3F2','4F2'};
                                r_table = sprintf('%s/%s_carrier_rc%d_carr%d_%s_projected_%s', fig_location, opt.data_type, r, carriers(f), rcaType, trial_type);  
                            end
                            if opt.train_data > 0
                                r_table = sprintf('%s_%s.csv',r_table,train_stim);
                            else
                                r_table = sprintf('%s.csv',r_table);
                            end
                            anova_subs = size(anova_ready,1)/length(anova_names)/length(harm_names);
                            harm_idx = repmat(1:length(harm_names), anova_subs, 1);
                            harm_idx = repmat(harm_idx(:),length(anova_names), 1);
                            harm_lbl = arrayfun(@(x) harm_names(x), harm_idx);
                            
                            subj_lbl = arrayfun(@(x) ...
                                num2str(x, 's%.02d'), repmat(1:anova_subs,1,length(anova_names)*length(harm_names)), 'uni', false);
                            tbl = table(cond_lbl, subj_lbl', harm_lbl, anova_ready,'VariableNames', {'condition','subject','harmonic','data'});
                            tbl.subject = categorical(tbl.subject);
                            tbl.condition = categorical(tbl.condition);
                            tbl.harmonic = categorical(tbl.harmonic);
                            writetable(tbl, r_table, 'Delimiter',',','QuoteStrings',true);
                        else
                        end  
                    else
                        amp_vals = arrayfun(@(x) curRCA(x).mean.amp_signal(1, 1, r, cur_order(c)), curIdx);
                        err_ub = arrayfun(@(x) curRCA(x).stats.amp_up_err(1, 1, r, cur_order(c)), curIdx);
                        err_lb = arrayfun(@(x) curRCA(x).stats.amp_lo_err(1, 1, r, cur_order(c)), curIdx);
                        within_sig = arrayfun(@(x) curRCA(x).stats.t2_p(1, 1, r, cur_order(c))<0.05, curIdx);
                    end
                    
                    if c == 1
                        y_max = ceil(max(cell2mat(...
                            arrayfun(@(x) max(curRCA(x).mean.amp_signal(1, 1, r,:)+curRCA(x).stats.amp_up_err(1, 1, r,:)), curIdx,'uni',false))));
                    else
                    end
                    errorbar(x_vals(c,:), amp_vals, err_lb, err_ub, '.k', 'LineWidth',lWidth, 'marker','none');
                    
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
                if opt.use_projected
                    cur_p = squeeze(between_tstatP(curIdx, r, opt.plot_split+1, test_idx));
                    between_idx = (cur_p < 0.05);
                    between_idx = between_idx + (cur_p < 0.005);
                else
                    cur_p = squeeze(between_t2p(curIdx,r,opt.plot_split+1,test_idx));
                    cur_sign = squeeze(between_sign(curIdx,r,opt.plot_split+1,test_idx) == 1);
                    between_idx = (cur_p < 0.05) & (cur_sign == 1);
                    between_idx = between_idx + ((cur_p < 0.005) & (cur_sign == 1));
                end
                if do_exp == 2
                    for c = 1:3
                        if any(between_idx(:,c))
                            sig_x = [-.2, .2, 0];
                            sig_y = [1.05,1.05,.95];
                            sig_labels = {'3','1','N'};
                            freq = find(between_idx(:,c));
                            arrayfun(@(x) text(x+sig_x(c), y_max*sig_y(c),sig_labels{c},'fontsize',10,'HorizontalAlignment','center'), freq,'uni',false);
                        end
                    end
                else
                    if any(between_idx)
                        arrayfun(@(x) text(x,y_max*.95,'*','fontsize',20,'HorizontalAlignment','center'), find(between_idx > 0),'uni',false);
                        arrayfun(@(x) text(x,y_max*.85,'*','fontsize',20,'HorizontalAlignment','center'), find(between_idx == 2),'uni',false);
                    else
                    end
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
                elseif r==yFigs
                    if t==2
                        lH = legend(amp_h, cond_names);
                        legend boxoff
                        l_pos = get(lH,'position');
                        l_pos(1) = l_pos(1) + .10;
                        l_pos(2) = l_pos(2) + .13;
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
        for r = 1:yFigs
            for z = 1:size(egiH,2)
                if opt.plot_split
                    addVal = 0.8;
                    shiftLeft = 0.02;
                else
                    addVal = 1.5;
                    shiftLeft = 0.01;
                end
                newPos = get(egiH(r,z),'position');
                newPos(1) = newPos(1)-(newPos(3)*addVal/2) - shiftLeft;
                newPos(2) = newPos(2)-(newPos(4)*addVal/2)+(newPos(4)*addVal/30)*(r);
                newPos(3:4) = newPos(3:4)*(1+addVal);
                set(egiH(r,z),'position',newPos);
            end
        end
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = figHeight;
        figPos(3) = figWidth;
        set(gcf,'pos',figPos);

        if opt.use_projected    
            fig_name = sprintf('%s/%s_carr%d_%s_projected_%s', fig_location, opt.data_type, carriers(f), rcaType, trial_type);
        else
            fig_name = sprintf('%s/%s_carr%d_%s_%s', fig_location, opt.data_type, carriers(f), rcaType, trial_type);
        end
        
        if opt.train_data > 0
            fig_name = sprintf('%s_%s',train_stim);
        else
        end

        export_fig(sprintf('%s.png', fig_name),'-png','-opengl','-m5','-transparent',gcf);
    end
%% Test plot just one condition
% close all
% use_freqs= 1:8;
% t = 1; % 1 = odd; 2 = carrier 
% numFreqs = max(use_freqs)/2;
% curIdx = (use_freqs(1:end/2))+(t-1)*numFreqs;
% f = 1; % carrier either 1 = 6 or 2 = 8
% rc = 2; %rca component number
% new_orders = [[5 4 6]; [2 3 1]]; %6v6,6v5,6v9; 8v8,8v9,8v5
% cur_order = cur_order;
% hold on
% for c = 1:3
%     cur_order = cur_order;
%     curRange = ampVals(curIdx,rc,cur_order,1,1) + errUB(curIdx,rc,cur_order,1,1);
%     ampH(c)=plot(1:numFreqs,ampVals(curIdx,rc,cur_order(c),1,opt.plot_split+1),'-','LineWidth',lWidth,'Color',subColors(c,:));
%     errorbar(1:numFreqs,ampVals(curIdx,rc,cur_order(c),1,1),errLB(curIdx,rc,cur_order(c),1,1),errUB(curIdx,rc,cur_order(c),1,1),'Color',subColors(c,:),'LineWidth',lWidth);
%     yMax = ceil(max(curRange(:)));
%     zeroSig = tVs0Pval(curIdx,rc,cur_order(c),1,1)<0.05;
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
%     %curSig = curSig+(tPval(curIdx,r,1,opt.plot_split+1,:,f)<0.005);
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
