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
    %                     [0](vector)/1(amplitude only)/2(zSNR)
    %
    %   behSplit        - only use correct trials
    %                     true/[false]
    
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
            'chanToCompare'	, 75, ...
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
    top_folder = '/Volumes/Denali_DATA1/kohler/EEG_EXP/DATA/numeroOddball';
    do_exp = 2;

    %% IDENTIFY DATA LOCATION
    data_location = sprintf('%s/Experiment%.0d',top_folder,do_exp);    
    fig_location = sprintf('%s/Experiment%.0d/figures',top_folder,do_exp);
    if ~exist(fig_location,'dir')
        mkdir(fig_location);
    else
    end
    folder_names=subfolders(sprintf('%s/*20*',data_location),1);
    if ~isempty(opt.foldersToUse)
        folder_names = folder_names(opt.foldersToUse);
    else
    end
    % and string for saving the data
    save_str = datestr(clock,26);
    save_str(strfind(save_str,'/')) ='';
    save_name = sprintf('%s/rcaData_%s.mat',data_location, save_str); % include the date as a string;

    %% SETUP INPUTS
    use_bins = 0; % use average bin, the zeroeth one
    use_freqs = 1:8; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
   
    train_data = opt.trainData;
    if train_data == 0
        use_conds = 1:6;
        %Reorder to make plotting easier
        carriers = [6, 8];
        new_orders = [[5 4 6]; [2 3 1]]; %6v6,6v5,6v9; 8v8,8v9,8v5
    elseif train_data == 1
        use_conds = 4:6;
        new_orders = [2 1 3];
        carriers = 6;
        train_stim = 'train6';
    elseif train_data == 2
        use_conds = 1:3;
        new_orders = [2 3 1];
        carriers = 8;
        train_stim = 'train8';
    end
    
    use_trials = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    n_reg = 10; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    n_comp = 5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    n_freq = length(use_freqs);
    n_cond = length(use_conds);
    data_type = 'RLS'; % can also be 'DFT' if you have DFT exports
    rc_plotstyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

    %% LAUNCH RCA OR LOAD DATA

    if ~opt.launchRCA % if users says analysis is not to be launched
        temp_name = subfiles([data_location,'/rcaData_*.mat'],1);
        temp_name = temp_name{end}; % grab newest file and loose .mat suffix;
        load(temp_name);
    else
        % generate "filter" for comparison data
        w_comparison = zeros(128,1); w_comparison(opt.chanToCompare) = 1; 
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
        rca_full = rcaSweep(path_names, use_bins, use_freqs, use_conds, use_trials, n_reg, n_comp, data_type, opt.chanToCompare, [], rc_plotstyle, opt.forceSourceData); % TRUE

        rca_h = grabCovFig(gcf);
        export_fig(sprintf('%s/full_rca_cond%s_cov.pdf',fig_location,sprintf('%.0d',use_conds)),'-pdf','-transparent',rca_h);
        close all;
        % create axxRCA
        full_w = [rca_full(1).W, w_comparison];
        axx_rca_full = axxRCAmake(path_names, full_w, use_conds, opt.axxType); % create struct
        
        % oddball RCA, first two frequencies 1F1 and 2F1
        % since this is just a subset of the previous RCA, set forceSourceData to false
        rca_odd = rcaSweep(path_names, use_bins, use_freqs(1:end/2), use_conds, use_trials, n_reg, n_comp, data_type, opt.chanToCompare,[], rc_plotstyle, false);
        rca_h = grabCovFig(gcf);
        export_fig(sprintf('%s/oddRCA_cond%s_cov.pdf', fig_location, sprintf('%.0d',use_conds)),'-pdf','-transparent',rca_h);
        close all;
        odd_w = [rca_odd(1).W, w_comparison];

        % carrier RCA, last two frequencies 1F2 and 2F2
        rca_carr = rcaSweep(path_names, use_bins, use_freqs(end/2+1:end), use_conds, use_trials, n_reg, n_comp, data_type, opt.chanToCompare,[], rc_plotstyle, false);
        rca_h = grabCovFig(gcf);
        export_fig(sprintf('%s/carrierRCA_cond%s_cov.pdf', fig_location, sprintf('%.0d',use_conds)),'-pdf','-transparent',rca_h);
        
        carr_w = [rca_carr(1).W, w_comparison];
        
        axx_rca_odd = axxRCAmake(path_names, odd_w, use_conds, opt.axxType); % create struct
        axx_rca_carr = axxRCAmake(path_names, carr_w, use_conds, opt.axxType); % create struct

        close all;
        save(save_name,'rca_full', 'rca_odd', 'rca_carr', 'axx_rca_*', '-append')
    end
    
    %% NOW SPLIT BASED ON BEHAVIOR
    if opt.behSplit 
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
    
    %% COMPUTE VALUES FOR PLOTTING

    keepConditions = true;
    errorType = 'SEM';
    doNR = false(1,6,6); % 1 freq, 5 RCs, 6 conditions
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
                rcStruct = aggregateData(curRCA(c,f),keepConditions,errorType,opt.trialError,doNR);
                % RC
                curRCA(c,f).stats.Amp = squeeze(rcStruct.ampBins(1,:,:,:));
                curRCA(c,f).stats.SubjectAmp = squeeze(rcStruct.subjectAmp(1,:,:,:,:));
                curRCA(c,f).stats.ErrLB = squeeze(rcStruct.ampErrBins(1,:,:,:,1));
                curRCA(c,f).stats.ErrUB = squeeze(rcStruct.ampErrBins(1,:,:,:,2));
                curRCA(c,f).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins(1,:,:,:));
                curRCA(c,f).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise(1,:,:,:,:));
                % within group t-values
                curRCA(c,f).stats.tSqrdP = squeeze(rcStruct.tSqrdP(1,:,:,:));
                curRCA(c,f).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig(1,:,:,:));
                curRCA(c,f).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal(1,:,:,:));
            end
            % AxxRCA
            n_tps = size(curAxxRCA(c).Projected{1,1},1);
            n_rca = size(curAxxRCA(c).Projected{1,1},2);
            n_subs = size(curAxxRCA(c).Projected,2);
            curAxxRCA(c).ProjMat = reshape(cell2mat(curAxxRCA.Projected),n_tps,n_cond,n_rca,n_subs);
            curAxxRCA(c).avg = nanmean(curAxxRCA(c).ProjMat,4);
            curAxxRCA(c).std = nanstd(curAxxRCA(c).ProjMat,0,4);
            curAxxRCA(c).sem = curAxxRCA(c).std./sqrt(n_subs);
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
    % think how to compute the appropriate pairs to do the
    % ttests and how they are stored
    contrast_order = [[4 5];[6 5];[6 4];[3 2];[1 2];[1 3]];
    store_order = [[1 1];[2 1];[3 1];[1 2];[2 2];[3 2]];
    if train_data == 1
        contrast_order = [[1 2];[3 2];[3 1]];
    elseif train_data == 2
        contrast_order = [[3 2];[1 2];[1 3]];
    else
    end
    % get the data and run the between group tests
    for f = 1:max(use_freqs)
        for rcType = 1:2
            if rcType == 1
                [rca_data_real,rca_data_imag] = getRealImag(rca_full(f).data);
                [comp_data_real,comp_data_imag] = getRealImag(rca_full(f).comparisonData);
            else
                if f > length(rca_odd)
                    [rca_data_real,rca_data_imag] = getRealImag(rca_carr(f-length(rca_odd)).data);
                    [comp_data_real,comp_data_imag] = getRealImag(rca_carr(f-length(rca_odd)).comparisonData);
                else
                    [rca_data_real,rca_data_imag] = getRealImag(rca_odd(f).data);
                    [comp_data_real,comp_data_imag] = getRealImag(rca_odd(f).comparisonData);
                end
            end
            % concatenate comparison and rca data
            rca_data_real = cellfun(@(x,y) cat(2,x,y), rca_data_real, comp_data_real, 'uni', false);
            rca_data_imag = cellfun(@(x,y) cat(2,x,y), rca_data_imag, comp_data_imag, 'uni', false);
            
            for r = 1:6
                x_data = cell2mat(cellfun(@(x) squeeze(nanmean(x(1,r,:),3)), rca_data_real, 'uni', false));
                y_data = cell2mat(cellfun(@(x) squeeze(nanmean(x(1,r,:),3)), rca_data_imag, 'uni', false));
                xy_data = cat(2,permute(x_data,[2,3,1]),permute(y_data,[2,3,1]));
                % compute vector mean amplitude, to get sign
                temp_amp = sqrt(nanmean(x_data,2).^2+nanmean(y_data,2).^2);
                for t=1:size(contrast_order,1) % number of tTests
                    results = tSquaredFourierCoefs(squeeze(xy_data(:,:,contrast_order(t,:))));
                    between_t2sig(f,r,rcType,store_order(t,1),store_order(t,2)) = results.H;
                    between_t2p(f,r,rcType,store_order(t,1),store_order(t,2)) = results.pVal;
                    between_t2Stat(f,r,rcType,store_order(t,1),store_order(t,2)) = results.tSqrd;
                    between_sign(f,r,rcType,store_order(t,1),store_order(t,2)) = sign(temp_amp(contrast_order(t,1))-temp_amp(contrast_order(t,2)));
%                     % amplitude paired ttest
%                     [h,p,ci,ampResults] = ttest(ampData(:,contrast_order(t,1)),ampData(:,contrast_order(t,2)));
%                     AmptSig(f,r,fPair,rcType,store_order(t,1),store_order(t,2)) = h;
%                     AmptPval(f,r,fPair,rcType,store_order(t,1),store_order(t,2)) = p;
%                     AmptStat(f,r,fPair,rcType,store_order(t,1),store_order(t,2)) = ampResults.tstat;
%                     %zSNR paired ttest
%                     [h,p,ci,zSNRResults] = ttest(zSNRData(:,contrast_order(t,1)),zSNRData(:,contrast_order(t,2)));
%                     zSNRtSig(f,r,fPair,rcType,store_order(t,1),store_order(t,2)) = h;
%                     zSNRtPval(f,r,fPair,rcType,store_order(t,1),store_order(t,2)) = p;
%                     zSNRtStat(f,r,fPair,rcType,store_order(t,1),store_order(t,2)) = zSNRResults.tstat;
                end
            end
        end
    end


    %% NEW PLOTTING
    close all
    plotSNR = false;
    x_vals = rca_full(1).settings.freqLabels; % frequencies are x-values
    % Which dimension corresponds to what in axxProjMat
    tPtsIdx = 1;
    condsIdx = 2;
    RCAIdx = 3;
    
    % axxProjMat.avg = axxProjMat(:,newOrder,:);
    % axxProjMat.std = axxProjMat(:,newOrder,:);
    % axxProjMat.sem = axxProjMat(:,newOrder,:);
    % axxProjMat.ProjMat = axxProjMat.ProjMat(:,newOrder,:,:);
    axxTicks = {[0 500 1000 1500 2000], [0 500 1000 1500 2000]};
    bgColors = [1,1,1; 0,0,0];
    subColors =  repmat(distinguishable_colors(4,bgColors),2,1);
    subColors = subColors(1:6,:); %Change this from use_conds for consistency in colors between plots
    lWidth = 1;
    fSize = 10;
    gcaOpts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

    yFigs = n_comp+1;

    if opt.plotSplit == 0    
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
    barH = [0,0];
    bar_width = .3;
    yUnit = 1;
    
    for f = 1:length(carriers) % carriers (i.e. 6 & 8)
        % Axx xVals
        nTps = size(axx_rca_full.ProjMat,tPtsIdx);
        carrierFreq = [3.0 3.0];
        nReps = [6 6];
        axx_xvals = linspace(0,1000/carrierFreq(f)*nReps(f),nTps+1);
        axx_xvals = axx_xvals(2:end);
        stimOnset = axx_xvals(floor(linspace(1,length(axx_xvals),nReps(f)+1)));
        stimOnset = stimOnset(1:end-1);
        figure;
        for r = 1:yFigs 
            if r < 6            
                if opt.plotSplit == 0
                    curRCA = rca_full;
                    egiH(r,1) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs);
                    hold on
                    rcaColorBar = [min(curRCA(1).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA(1).A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'full';
                else
                    curRCA = rca_split;
                    egiH(r,1) = subplot(yFigs,xFigs,3+(r-1)*xFigs); %New Candidate
                    hold on
                    rcaColorBar = [min(curRCA(1).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA(1).A(:,r),rcaColorBar);
                    hold off
                    egiH(r,2) = subplot(yFigs,xFigs,(xFigs-2)+(r-1)*xFigs); %I think this one is correct
                    hold on
                    rcaColorBar = [min(curRCA(5).A(:,r)),max(curRCA(1).A(:,r))];
                    newExtreme = max(abs(rcaColorBar));
                    rcaColorBar = [-newExtreme,newExtreme];
                    mrC.plotOnEgi(curRCA(5).A(:,r),rcaColorBar);
                    hold off
                    rcaType = 'split';
                end
            else
            end

            % axx plot
            for z = 1:2
                if opt.plotSplit == 0
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
                dataToPlot = squeeze(cur_axx.avg(:,new_orders(f,:),r));
                errorToPLot = squeeze(cur_axx.sem(:,new_orders(f,:),r));
                axx_ymin = floor(min((min(cur_axx.avg(:,new_orders(f,:),r)-cur_axx.sem(:,new_orders(f,:),r)))/2))*2;
                axx_ymax = ceil(max((max(cur_axx.avg(:,new_orders(f,:),r)+cur_axx.sem(:,new_orders(f,:),r)))/2))*2;
                if ( axx_ymax-axx_ymin ) > 10
                    axx_yunit = 2;
                else
                    axx_yunit = 1;
                end
                hold on
                for c=1:size(new_orders(f,:),2)
                    plot(axx_xvals,dataToPlot(:,c),'-','LineWidth',lWidth,'Color',bold_colors(c,:));
                    ErrorBars(axx_xvals',dataToPlot(:,c),errorToPLot(:,c),'color',bold_colors(c,:));
                end
                ylim([axx_ymin, axx_ymax]);
                xlim([0,axx_xvals(end)]);
                set(gca,gcaOpts{:},'xtick',axxTicks{f},'xticklabel',cellfun(@(x) num2str(x),num2cell(axxTicks{f}),'uni',false),'ytick',axx_ymin:axx_yunit:axx_ymax);
                yLine = repmat(get(gca,'YLim'),nReps(f),1)';
                line(repmat(stimOnset,2,1),yLine,'Color','black');
                if r == 1
                    title(title_str,'fontweight','normal');
                elseif r == 6
                    xlabel('time (ms)');
                else
                end        
                s = get(gca, 'Position');
                set(gca, 'Position', [s(1)+0.01, s(2), s(3)*0.9, s(4)]);
                if opt.plotSplit == 0
                    continue
                else
                end
            end

            for t = 1:2
                subplot(yFigs,xFigs,(xFigs-4)+(r-1)*xFigs+(t-1));
                hold on
                if opt.plotSplit
                    cur_order = new_orders(f,:);
                else
                    cur_order = new_orders(f,:);
                end
                freq_set = max(use_freqs)/2;
                cond_set = length(cur_order);
                curIdx = (use_freqs(1:end/2))+(t-1)*freq_set;
                x_vals = repmat((1:freq_set),cond_set,1) + repmat(linspace(-bar_width,bar_width,cond_set),freq_set,1)';
                for c=1:length(cur_order)
                    if plotSNR
                        curRange = snrVals(curIdx,:,cur_order,1,opt.plotSplit+1);
                        %valSet = snrVals(curIdx,r,c,f,opt.plotSplit+1)
                      
                    else
                        amp_vals = arrayfun(@(x) curRCA(x).stats.Amp(r,cur_order(c)), curIdx);
                        err_ub = arrayfun(@(x) curRCA(x).stats.ErrUB(r,cur_order(c)), curIdx);
                        err_lb = arrayfun(@(x) curRCA(x).stats.ErrLB(r,cur_order(c)), curIdx);
                        errorbar(x_vals(c,:), amp_vals, err_lb, err_ub, '.k', 'LineWidth',lWidth, 'marker','none');
                        if c == 1
                            y_max = ceil(max(cell2mat(...
                                arrayfun(@(x) max(curRCA(x).stats.Amp(r,cur_order)+curRCA(x).stats.ErrUB(r,cur_order)), curIdx,'uni',false))));
                        else
                        end
                    end
                    within_sig = arrayfun(@(x) curRCA(x).stats.tSqrdP(r,cur_order(c))<0.05, curIdx);
                    
                    amp_h(c) = bar(x_vals(c,:), amp_vals,'BarWidth',bar_width,'edgecolor','none','facecolor',bold_colors(c,:));

                    for z = 1:length(within_sig)
                        if within_sig(z)
                            patch_x = [ x_vals(c,z)-bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)+bar_width/2, x_vals(c,z)-bar_width/2];
                            patch_y = [0, 0, y_max, y_max];
                            pa_h = patch(patch_x, patch_y, bold_colors(c,:),'edgecolor','none','facealpha',.25); 
                            uistack(pa_h,'bottom');
                            %amp_h(c) = bar(x_vals(c,z), amp_vals(z), 'bar_width',bar_width, 'edgecolor', 'none', 'facecolor', weak_colors(c,:));
                        else
                        end
                    end

                    if opt.ampTest == 0
                        cur_p = squeeze(between_t2p(curIdx,r,opt.plotSplit+1,:,f));
                    elseif opt.ampTest == 1
                        cur_p = squeeze(AmptPval(curIdx,r,opt.plotSplit+1,:,f));
                    elseif opt.ampTest == 2
                        cur_p = squeeze(zSNRtPval(curIdx,r,opt.plotSplit+1,:,f));
                    end
                    cur_sign = squeeze(between_sign(curIdx,r,opt.plotSplit+1,:,f) == 1);
                    between_idx = (cur_p(:,c) < 0.05) & (cur_sign(:,c) == 1);
                    if any(between_idx)
                        sig_x = [-.2, 0, .2];
                        sig_y = [1.05,1.05,.95];
                        sig_labels = {'1','3','N'};
                        freq = find(between_idx);
                        arrayfun(@(x) text(x+sig_x(c), y_max*sig_y(c),sig_labels{c},'fontsize',10,'HorizontalAlignment','center'), freq,'uni',false);
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
                elseif r==6
                    if t==2
                        lH = legend(amp_h,{'control','dist 1','dist 3'});
                        legend boxoff
                        l_pos = get(lH,'position');
                        l_pos(1) = l_pos(1) + .10;
                        l_pos(2) = l_pos(2) + .05;
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
        for r = 1:5
            for z = 1:size(egiH,2)
                if opt.plotSplit
                    addVal = 0.72;
                    shiftLeft = 0.02;
                else
                    addVal = 0.8;
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
        end
        if plotSNR        
            fig_name = sprintf('%s/%s_car%d_%s_snr.pdf',fig_location,data_type,carriers(f),rcaType);
        else
            if train_data == 0
                fig_name = sprintf('%s/%s_car%d_%s_%s_%s.pdf',fig_location,data_type,carriers(f),rcaType,testplot,trialLab);
            elseif train_data == 1
                fig_name = sprintf('%s/%s_car%d_%s_%s_%s_%s.pdf',fig_location,data_type,carriers(f),rcaType,train_stim,testplot,trialLab);
            elseif train_data == 2
                fig_name = sprintf('%s/%s_car%d_%s_%s_%s_%s.pdf',fig_location,data_type,carriers(f),rcaType,train_stim,testplot,trialLab);
            else
            end
        end
        export_fig(fig_name,'-pdf','-painters','-transparent',gcf);
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
% cur_order = new_orders(f,:);
% hold on
% for c = 1:3
%     cur_order = new_orders(f,:);
%     curRange = ampVals(curIdx,rc,cur_order,1,1) + errUB(curIdx,rc,cur_order,1,1);
%     ampH(c)=plot(1:numFreqs,ampVals(curIdx,rc,cur_order(c),1,opt.plotSplit+1),'-','LineWidth',lWidth,'Color',subColors(c,:));
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
