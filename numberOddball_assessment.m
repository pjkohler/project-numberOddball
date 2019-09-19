clear all; close all;
load('xDiva_params.mat');

% add paths
addpath(genpath('/Users/kohler/code/git/mrC/tools/colors'))
addpath(genpath('/Users/kohler/code/git/export_fig'))
addpath(genpath('/Users/kohler/code/git/stimulus_assessment'));
addpath(genpath('./socmodel'))
addpath('./minbound_suite')
addpath(genpath('/Users/kohler/code/git/xDiva'))

fig_folder = '/Users/kohler/Google Drive/Dropbox/WRITING/Articles/2019_KohlerNumerositySSVEP/figures';

carr_vals = {'5','6','8','9'};
odd_vals = {'5','6','8','9'};

% experiment 1: carr8_odd5, carr8_odd8
% experiment 2: carr8_odd5, carr8_odd8, carr8_odd9, 
%               carr6_odd5, carr6_odd6, carr6_odd9
% experiment 3: carr5_odd5, carr5_odd8
%               carr6_odd6, carr6_odd9
%               carr8_odd8, carr8_odd5
%               carr9_odd9, carr9_odd6 

model_sets = {'carr5_odd5', 'carr5_odd8', 'carr6_odd5', 'carr6_odd6', 'carr6_odd9', 'carr8_odd5', 'carr8_odd8', 'carr8_odd9', 'carr9_odd6', 'carr9_odd9'};

load_values = false;

% run without prelude
timing{3,2}=0;

if ~load_values
    global img_out; global seq_out;
    img_array = [];
    for o = 1:length(odd_vals)
        parameters{3}{11} = odd_vals{o};
        for c = 1:length(carr_vals)
            parameters{4}{11} = carr_vals{c};
            cur_set = (sprintf('carr%s_odd%s', carr_vals{c}, odd_vals{o}));
            if ~ismember(cur_set ,model_sets)
                continue;
            else
            end
            close all;
            all_images = [];
            for z = 1:20
                pmf_FastOddball_Numerosity_testing('MakeMovie',parameters,timing,videoMode,1);
                % reorder images
                all_images = cat(3, all_images, squeeze(img_out(:,:,:,seq_out(seq_out>0))));
            end

            % make indices
            num_images = size(all_images,3); 
            odd_even_ratio = parameters{3}{1,2}/parameters{4}{1,2};
            val_idx = zeros(1, num_images);
            val_idx(odd_even_ratio:odd_even_ratio:end) = 1;
            diff_idx = zeros(1, num_images-1);
            diff_idx(odd_even_ratio-1:odd_even_ratio:end) = 1;

            odd.(cur_set).images = all_images(:,:, val_idx == 1);
            carr.(cur_set).images = all_images(:,:, val_idx == 0);
            
            % stimulus values
            [stim_vals, e1] = number_assessment(all_images);
            
            % compute euclidean distances
            % z-score across oddball and carrier (4 stimulus params, excluding number)
            z_vals = zscore(stim_vals(:,1:4), 0, 1);
            % compute the Euclidean norm across 4 params, for each sample
            dist_vals = vecnorm(z_vals,2,2);
            % add to stim values
            stim_vals(:,6) = stim_vals(:,5);
            stim_vals(:,5) = dist_vals; 
            stim_diff = diff(stim_vals);
            
            if length(unique(stim_vals(:,6))) > 2
                msg = '\n more than two numerosities in stimulus set \n';
                error(msg);
            else
            end
            
            % response values
            our_fov = 10;
            [resp_vals, e2] = soc_assessment(all_images, our_fov);
            resp_diff = diff(resp_vals);
  
            % extract values and differences for odd and carrier
            odd.(cur_set).stim.val = stim_vals(val_idx==1, :);
            odd.(cur_set).stim.diff = stim_diff(diff_idx==1, :);    
            carr.(cur_set).stim.val = stim_vals(val_idx==0, :);
            carr.(cur_set).stim.diff = stim_diff(diff_idx==0, :);
            odd.(cur_set).resp.val = resp_vals(val_idx==1, :);
            odd.(cur_set).resp.diff = resp_diff(diff_idx==1, :);
            carr.(cur_set).resp.val = resp_vals(val_idx==0, :);
            carr.(cur_set).resp.diff = resp_diff(diff_idx==0, :); 
            
            % examples
            example.(cur_set).input = e1.input;
            example.(cur_set).hull = e1.hull;
            example.(cur_set).resp = e2.resp;
            clear all_images resp_* stim_*;
        end
    end
    save(sprintf('%s/assessment_data.mat', fig_folder), 'odd', 'carr', 'example', '-v7.3');
else
    load(sprintf('%s/assessment_data.mat', fig_folder), 'odd', 'carr', 'example');
end

%% figure params
f_size = 12;
gcaOpts = {'tickdir', 'out', 'fontsize', f_size, 'fontname', 'Helvetica', 'linewidth', 1,'box','off', 'TickLength',[0.05, 0.05]};
cBrewer = load('colorBrewer_new.mat');
odd_color = cBrewer.rgb20(7,:);
carr_color = cBrewer.rgb20(15,:);
odd_string = ['{\bf\color[rgb]{',num2str(odd_color,' %0.2f'), '}oddball}'];
carr_string = ['{\bf\color[rgb]{',num2str(carr_color,' %0.2f'), '}carrier}'];
text_params = {'fontsize', f_size, 'fontname', 'Helvetica','fontweight','normal', 'horizontalalignment','center'};

%% generate example images
close all;
figure;
subplot(1,3,1);
imagesc(example.carr5_odd5.input(:,:,1),[0 254]); axis image tight; colormap(gray); title('input', text_params{:}); axis off;
subplot(1,3,2);
imagesc(example.carr5_odd5.hull,[0 254]); axis image tight; colormap(gray); title('stimulus (with convex hull)', text_params{:}); axis off;
subplot(1,3,3);
imagesc(example.carr5_odd5.resp(:,:,1),[0 1.5]); axis image tight; colormap(gray); title('V1 response', text_params{:}); axis off;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 30; fig_pos(4) = 10;
set(gcf,'pos',fig_pos);
export_fig(sprintf('%s/examples/stimulus_example.pdf',fig_folder),'-transparent', '-nocrop');

close all;
figure;
for c = 1:6
    if c == 6
        imshow(odd.carr8_odd5.images(:,:,1));
    else
        imshow(carr.carr8_odd5.images(:,:,c));
    end
    axis off
    set(gcf,'units','centimeters');
    fig_pos = get(gcf,'pos');
    fig_pos(3) = 10; fig_pos(4) = 10;
    set(gcf,'pos',fig_pos);
    export_fig(sprintf('%s/examples/stimulus%d.png',fig_folder, c),'-png','-opengl','-m5','-transparent',gcf);
end
close all

%% 

%% size and area plotting, Experiment 1
% close all
% figure;
% for z = 1:2
%     subplot(1,2,z);
%     if z == 1
%         cur_set = 'set_8vs5';
%         title_str = '5 odd vs 8 carr';
%     else
%         cur_set = 'set_8vs8';
%         title_str = '8 odd vs 8 carr';
%     end
%     hold on
%     plot(carr.(cur_set).size.diff, carr.(cur_set).area.diff,'bo'); 
%     plot(odd.(cur_set).size.diff, odd.(cur_set).area.diff,'ro'); 
%     x_min = -1; x_max = 1; x_units = .5;
%     y_min = -30; y_max = 30; y_units = 10;
%     xlim([x_min,x_max]);
%     ylim([y_min,y_max]);
%     set(gca,'xtick', x_min:x_units:x_max, 'ytick', y_min:y_units:y_max , gcaOpts{:})
%     if z == 1
%         ylabel('area (% of image size)')
%         xlabel('dot size (% of image size)');
%     else
%     end
%     axis square
%     text(x_min*.95, y_max*.95, title_str, 'fontsize', f_size, 'fontname', 'Helvetica')
% end
% set(gcf,'units','centimeters');
% fig_pos = get(gcf,'pos');
% fig_pos(3) = 20; fig_pos(4) = 10;
% set(gcf,'pos',fig_pos);

%% do correlations
y_label = {'dot size', 'total dot area', 'convex hull', 'mean occupancy', 'euclidean distance', 'numerosity'};
all_carr_stim = []; all_odd_stim = []; all_carr_resp = []; all_odd_resp = [];
model_sets = {'carr5_odd5', 'carr5_odd8', 'carr6_odd5', 'carr6_odd6', 'carr6_odd9', 'carr8_odd5', 'carr8_odd8', 'carr8_odd9', 'carr9_odd6', 'carr9_odd9'};
for m = 1:length(model_sets)
    cur_set = model_sets{m};
    all_carr_stim = cat(1, all_carr_stim, carr.(cur_set).stim.val);
    all_odd_stim = cat(1, all_odd_stim, odd.(cur_set).stim.val);
    all_carr_resp = cat(1, all_carr_resp, carr.(cur_set).resp.val);
    all_odd_resp = cat(1, all_odd_resp, odd.(cur_set).resp.val);
end
figure;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 48; fig_pos(4) = 8;
set(gcf,'pos',fig_pos);
for z = 1:6
    corr_r2(z) = corr([all_carr_stim(:,z); all_odd_stim(:,z)], ...
         [all_carr_resp(:,1); all_odd_resp(:,1)]).^2;
    subplot(1,6,z);
    hold on
    plot(all_carr_stim(:,z), all_carr_resp(:,1),'o','color', carr_color)
    plot(all_odd_stim(:,z), all_odd_resp(:,1),'o','color', odd_color)
    title(sprintf('%s \n (r^{2} = %.2f)',y_label{z}, corr_r2(z)), text_params{:});
    switch z
        case 6
            xlabel('# dots',text_params{:});
        case 5
            xlabel('normalized distance',text_params{:});
        case 4
            xlabel('convex hull/# dots',text_params{:});
        otherwise
            xlabel('% of image pixels',text_params{:});
            if z == 1
                ylabel('SOC V1 response',text_params{:});
            else
            end
    end
    set(gca, 'ylim', [0,.15] , 'clipping','off', gcaOpts{:})
    % set x-ticks
    x_unit = min(diff(get(gca, 'xtick')));
    x_lim = get(gca, 'xlim');
    x_lim(1) = x_lim(1)-mod(x_lim(1),x_unit);
    if mod(x_lim(2),x_unit)
        x_lim(2) = x_lim(2)-mod(x_lim(2),x_unit)+x_unit;
    else
    end
    set(gca, 'xlim', x_lim, 'xtick', x_lim(1):x_unit:x_lim(2));        
    axis square
    hold off
end
export_fig(sprintf('%s/assessment/correlation.pdf',fig_folder),'-transparent');


%% full assessment, all stimuli 
close all;
plot_values = false;
y_label = {'dot size', 'total dot area', 'convex hull', 'mean occupancy', 'euclidean distance', 'numerosity'};
y_lims = {[-1,1,.5],[-8,8,4],[-40,40,20],[-4,4,2],[-6,6,2], [-4,4,2]};

figure
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 50; fig_pos(4) = 25;
set(gcf,'pos',fig_pos);

a_idx = [1,2,3,4,6];

for a = 1:length(a_idx)
    for z = 1:10
        switch z
            case 1
                cur_set = 'carr5_odd5';
                title_str = '5 odd, 5 carr';
                plot_spacing = .1;
            case 2
                cur_set = 'carr5_odd8';
                title_str = '8 odd, 5 carr';
                plot_spacing = .1;
            case 3
                cur_set = 'carr6_odd5';
                title_str = '5 odd, 6 carr';
                plot_spacing = .3;
            case 4
                cur_set = 'carr6_odd6';
                title_str = '6 odd, 6 carr';
                plot_spacing = .1;
            case 5
                cur_set = 'carr6_odd9';
                title_str = '9 odd, 6 carr';
                plot_spacing = .1;
            case 6 
                cur_set = 'carr8_odd9';
                title_str = '9 odd, 8 carr';
                plot_spacing = .3;
            case 7 
                cur_set = 'carr8_odd8';
                title_str = '8 odd, 8 carr';
                plot_spacing = .1;
            case 8
                cur_set = 'carr8_odd5';
                title_str = '5 odd, 8 carr';
                plot_spacing = .1;
            case 9
                cur_set = 'carr9_odd9';
                title_str = '9 odd, 9 carr';
                plot_spacing = .3;
            case 10
                cur_set = 'carr9_odd6';
                title_str = '6 odd, 9 carr';
                plot_spacing = .1;
            otherwise
        end
        subplot(length(a_idx),10,z+(a-1)*10);
        hold on
        if plot_values
            odd_plot = odd.(cur_set).stim.val(:,a_idx(a));
            carr_plot = carr.(cur_set).stim.val(:,a_idx(a));
        else
            odd_plot = odd.(cur_set).stim.diff(:,a_idx(a));
            carr_plot = carr.(cur_set).stim.diff(:,a_idx(a));
        end
        boxplot(odd_plot, 'positions',  1 , 'width', 1, 'color', odd_color, 'symbol', 'o' )
        boxplot(carr_plot, 'positions', 2, 'width', 1, 'color', carr_color, 'symbol', 'o' )
        x_min = 0; x_max = 3;
        if max(abs([odd_plot; carr_plot])) > 10
            y_max = ceil(max(abs([odd_plot; carr_plot]))/20)*20;
            y_units = 10;
        elseif max(abs([odd_plot; carr_plot])) > 1
            y_max = ceil(max(abs([odd_plot; carr_plot]))/2)*2;
            y_units = 1;
        else
            y_max = ceil(max(abs([odd_plot; carr_plot]))/.2)*.2;
            y_units = .1;
        end
        if plot_values
            y_min = 0;
        else
            y_min = -y_max;
        end
        if max([odd_plot; carr_plot]) == y_max
            y_max = y_max + y_units;
        else
        end
        if min([odd_plot; carr_plot]) == y_min
            y_min = y_min - y_units;
        else
        end
        while (y_max - y_min)/y_units >= 5
            y_units = y_units * 2;
        end
        while (y_max - y_min)/y_units <= 3
            y_units = y_units / 2;
        end
        
        set(gca, 'xtick', [], 'xticklabels', {''})            
        xlim([x_min,x_max]);
        ylim([y_min,y_max]);
        set(gca, 'ytick', y_min:y_units:y_max , 'clipping','off', gcaOpts{:})
        if a == 1
            title(title_str, text_params{:});
        else
        end 
        if z == 1
            switch a_idx(a)
                case 6
                    y_label{a_idx(a)} = sprintf('%s \n (# dots)',y_label{a_idx(a)});
                case 5
                     y_label{a_idx(a)} = sprintf('%s \n (normalized)',y_label{a_idx(a)});
                case 4
                    y_label{a_idx(a)} = sprintf('%s \n (convex hull/# dots)',y_label{a_idx(a)});
                otherwise
                    y_label{a_idx(a)} = sprintf('%s \n ',y_label{a_idx(a)});
                    y_label{a_idx(a)} = [y_label{a_idx(a)}, '(% of image pixels)'];
            end
            y_lh = ylabel(y_label{a_idx(a)}, text_params{:}, 'position', [x_min-1.5, y_min+(y_max-y_min)*.5, 0]);
             
        else
            set(gca, 'yticklabels', {''}) 
        end
        if a == 1
            box_pos(:,z,a) = get(gca,'position');
            if z > 1
                box_pos([2,3,4],z,a) = box_pos([2,3,4],1,a);
                box_pos(1,z,a) = box_pos(1,z-1,a)+box_pos(3,z-1,a)*(1+plot_spacing);
            else
            end 
        else
            box_pos(:,z,a) = get(gca,'position');
            % inherit box positions from scatter positions
            box_pos([1,3,4],z,a) = box_pos([1,3,4],z,1);
            if z > 1
                box_pos(2,z,a) = box_pos(2,1,a);
               
            else
            end 
            
        end
        set(gca,'position', box_pos(:,z,a));
        ax_h(z,a) = gca;
        hold off
    end
    g_ymin = min(cell2mat(arrayfun(@(x) min(get(x,'ytick')),ax_h(:,a),'uni',false)));
    g_ymax = max(cell2mat(arrayfun(@(x) max(get(x,'ytick')),ax_h(:,a),'uni',false)));
    g_yunit = max(cell2mat(arrayfun(@(x) min(diff(get(x,'ytick'))),ax_h(:,a),'uni',false)));
    arrayfun(@(x) set(x,'ytick',g_ymin:g_yunit:g_ymax, 'ylim',[g_ymin,g_ymax]), ax_h(:,a),'uni',false);
end
export_fig(sprintf('%s/assessment/full_assessment_all.pdf',fig_folder),'-transparent');

%% response, all stimuli

figure;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 25; fig_pos(4) = 25;
set(gcf,'pos',fig_pos);

y_label = {'dot size', 'total dot area', 'convex hull', 'mean occupancy', 'euclidean distance', 'numerosity'};

for z = 1:10
    switch z
        case 1
            cur_set = 'carr5_odd5';
            title_str = '5 odd, 5 carr';
            sub_pos = 2;
        case 2
            cur_set = 'carr5_odd8';
            title_str = '8 odd, 5 carr';
            sub_pos = 3;
        case 3
            cur_set = 'carr6_odd5';
            title_str = '5 odd, 6 carr';
            sub_pos = 4;
        case 4
            cur_set = 'carr6_odd6';
            title_str = '6 odd, 6 carr';
             sub_pos = 5;
        case 5
            cur_set = 'carr6_odd9';
            title_str = '9 odd, 6 carr';
             sub_pos = 6;
        case 6 
            cur_set = 'carr8_odd9';
            title_str = '9 odd, 8 carr';
            sub_pos = 7;
        case 7 
            cur_set = 'carr8_odd8';
            title_str = '8 odd, 8 carr';
            sub_pos = 8;
        case 8
            cur_set = 'carr8_odd5';
            title_str = '5 odd, 8 carr';
            sub_pos = 9;
        case 9
            cur_set = 'carr9_odd9';
            title_str = '9 odd, 9 carr';
            sub_pos = 11;
        case 10
            cur_set = 'carr9_odd6';
            title_str = '6 odd, 9 carr';
            sub_pos = 12;
        otherwise
    end
    subplot(4,3,sub_pos);
    hold on;
    % adjust values
    mean_values = mean(cat(1, odd.(cur_set).resp.val, carr.(cur_set).resp.val));
    odd_resp = odd.(cur_set).resp.diff./repmat(mean_values,size(odd.(cur_set).resp.diff,1),1).*100;
    carr_resp = carr.(cur_set).resp.diff./repmat(mean_values,size(carr.(cur_set).resp.diff,1),1).*100;
    
    hold on
    boxplot(odd_resp, 'positions',  (1:4)-.15 , 'width', .3, 'color', odd_color, 'symbol', 'o' )
    boxplot(carr_resp, 'positions', (1:4)+.15, 'width', .3, 'color', carr_color, 'symbol', 'o' )
    x_min = 0.5; x_max = 4.5; x_units = 1;
    y_min = -60; y_max = 60; y_units = 20;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', 1:4, 'ytick', y_min:y_units:y_max , 'xticklabels', {'V1','V2','V3','V4'}, 'clipping','off', gcaOpts{:})
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:})
    if sub_pos == 2
        text(x_min+(x_max-x_min)*.5, y_max+(y_max-y_min)*.2, sprintf('change from previous image update (%s vs %s)', odd_string, carr_string), text_params{:});
    elseif sub_pos == 11
        %tH = text(x_min-(x_max-x_min)*.3, y_max+(y_max-y_min)*.4, 'SOC model response (% of mean)', text_params{:});
        %set(tH, 'rotation', 90);
        ylabel([sprintf('%s\n', 'SOC model response'),'(% of mean)'], text_params{:});
    else
    end
    ax_pos(:,z) = get(gca,'position');
    if z > 1
        ax_pos(3:4,z) = ax_pos(3:4,1);
    end
    if ismember(sub_pos, [3,5,6,8,9,12])
        ax_pos(3:4,z) = ax_pos(3:4,z-1);
        ax_pos(2,z) = ax_pos(2,z-1);
        ax_pos(1,z) = ax_pos(1,z-1)+ax_pos(3,z-1)*1.3;
    else
    end
    if ismember(sub_pos, [5,8,11])
        ax_pos(1,z) = ax_pos(1,1);
    elseif ismember(sub_pos, [6,9,12])
        ax_pos(1,z) = ax_pos(1,2);
    elseif ismember(sub_pos, [4,7])
        ax_pos(1,z) = ax_pos(1,1)-ax_pos(3,1)*1.3;
    end
    set(gca,'position',ax_pos(:,z));
end
export_fig(sprintf('%s/assessment/response_all.pdf',fig_folder),'-transparent');
    
%% size and area plotting, Experiment 2
figure;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 30; fig_pos(4) = 40;
set(gcf,'pos',fig_pos);
for z = 1:6
    switch z
        case 1
            cur_set = 'carr8_odd5';
            title_str = '5 odd, 8 carr';
        case 2
            cur_set = 'carr8_odd8';
            title_str = '8 odd, 8 carr';
        case 3
            cur_set = 'carr8_odd9';
            title_str = '9 odd, 8 carr';
        case 4
            cur_set = 'carr6_odd5';
            title_str = '5 odd, 6 carr';
        case 5
            cur_set = 'carr6_odd6';
            title_str = '6 odd, 6 carr';
        case 6
            cur_set = 'carr6_odd9';
            title_str = '9 odd, 6 carr';
        otherwise
    end
    
    % plot size and area
    if z > 3
        subplot(4,3,z+3);
    else
        subplot(4,3,z);
    end
    
    hold on
    plot(carr.(cur_set).stim.diff(:,1), carr.(cur_set).stim.diff(:,3),'o', 'color', carr_color); 
    plot(odd.(cur_set).stim.diff(:,1), odd.(cur_set).stim.diff(:,3),'o', 'color', odd_color); 
    x_min = -1; x_max = 1; x_units = .5;
    y_min = -30; y_max = 30; y_units = 10;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', x_min:x_units:x_max, 'ytick', y_min:y_units:y_max , gcaOpts{:})
    if z ~= 1 && z ~= 4
        set(gca, 'yticklabels', {''});
    else
    end
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:})
    if z == 1
        text(x_min-.3*(x_max-x_min), y_max+.05*(y_max-y_min), 'A', text_params{:}, 'fontsize', f_size*2);
    elseif z == 2
        text(x_min+(x_max-x_min)*.5, y_max*1.3, sprintf('change from previous image update (%s vs %s)', odd_string, carr_string),text_params{:});
    elseif z == 4
        y_label_h = ylabel('area (% of image pixels)',text_params{:});
        set(y_label_h, 'Units', 'Normalized', 'Position', [-0.25, 0.5, 0]);
        xlabel('dot size (% of image pixels)',text_params{:});
        text(x_min-.3*(x_max-x_min), y_max+.05*(y_max-y_min), 'C', text_params{:}, 'fontsize', f_size*2);
    else
    end
    scatter_pos(:,z) = get(gca,'position');
    scatter_pos(1,z) = scatter_pos(1,z)-scatter_pos(3,z)*.3*(z-1);
    scatter_pos(3:4,z) = scatter_pos(3:4,z)*1.2;
    scatter_pos(2,z) = scatter_pos(2,z)-scatter_pos(4,z)*.2*(z>3);
    
    % use top positions for bottom row
    if z > 3
        scatter_pos(1,z) = scatter_pos(1,z-3);
    else
    end
    set(gca,'position', scatter_pos(:,z));
    axis square
    
    if z > 3
        subplot(4,3,z+6);
    else
        subplot(4,3,z+3);
    end
    
    hold on
    for p = 1:2
        if p == 1
            cur_color = odd_color;
            cur_dist = odd.(cur_set).stim.diff(:,5);
        else
            cur_color = carr_color;
            cur_dist = carr.(cur_set).stim.diff(:,5);
        end
        b_h = boxplot(gca, cur_dist, 'positions', p ,'symbol', 'o',  'width', .75, 'color', cur_color, 'boxstyle' ,'outline');
        set(findobj(b_h,'tag','Upper Whisker'),'linestyle', '-', 'color', cur_color);
        set(findobj(b_h,'tag','Lower Whisker'),'linestyle', '-', 'color', cur_color);
        set(findobj(b_h, 'tag', 'Outliers'), 'visible','off');
        tag_outl = findobj(b_h, 'tag', 'Outliers');
        num_outl = length(tag_outl.XData);
        if num_outl > -1
            perc_outl = num_outl/length(cur_dist)*100;
            text(p, -.4, [sprintf('%.1f', perc_outl), '%'], text_params{:}, 'fontsize', 10);
            text(p, -.8, 'outliers', text_params{:}, 'fontsize', 10);
            
        else
        end
    end
    hold on
    x_min = 0.5; x_max = 2.5; x_units = 1;
    y_min = -1; y_max = 4; y_units = 1;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', 1:2, 'ytick', y_min:y_units:y_max , 'xticklabels', {'odd','carrier'}, 'clipping','off', gcaOpts{:})
    if z ~= 1 && z ~= 4
        set(gca, 'yticklabels', {''});
    else
    end
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:})
    axis square
    if z == 1
        text(x_min-.3*(x_max-x_min), y_max+.05*(y_max-y_min), 'B', text_params{:}, 'fontsize', f_size*2);
    elseif z == 4
        y_label_h = ylabel('normalized distance', text_params{:});
        set(y_label_h, 'Units', 'Normalized', 'Position', [-0.25, 0.5, 0]);
        text(x_min-.3*(x_max-x_min), y_max+.05*(y_max-y_min), 'D', text_params{:}, 'fontsize', f_size*2);
    else
    end
    axis square
    box_pos(:,z) = get(gca,'position');
    % inherit box positions from scatter positions
    box_pos([1,3,4],z) = scatter_pos([1,3,4],z);
    box_pos(2,z) = box_pos(2,z)-box_pos(4,z)*(.15 + .30*((z>3)));
    set(gca,'position', box_pos(:,z));
end
export_fig(sprintf('%s/assessment/size_area_exp2.pdf',fig_folder),'-transparent');

%% response plotting, Experiment 2
figure;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 30; fig_pos(4) = 20;
set(gcf,'pos',fig_pos);
for z = 1:6
    switch z
        case 1
            cur_set = 'carr8_odd5';
            title_str = '5 odd, 8 carr';
        case 2
            cur_set = 'carr8_odd8';
            title_str = '8 odd, 8 carr';
        case 3
            cur_set = 'carr8_odd9';
            title_str = '9 odd, 8 carr';
        case 4
            cur_set = 'carr6_odd5';
            title_str = '5 odd, 6 carr';
        case 5
            cur_set = 'carr6_odd6';
            title_str = '6 odd, 6 carr';
        case 6
            cur_set = 'carr6_odd9';
            title_str = '9 odd, 6 carr';
        otherwise
    end
    subplot(2,3,z);
    % adjust values
    mean_values = mean(cat(1, odd.(cur_set).resp.val, carr.(cur_set).resp.val));
    odd_resp = odd.(cur_set).resp.diff./repmat(mean_values,size(odd.(cur_set).resp.diff,1),1).*100;
    carr_resp = carr.(cur_set).resp.diff./repmat(mean_values,size(carr.(cur_set).resp.diff,1),1).*100;
    
    hold on
    boxplot(odd_resp, 'positions',  (1:4)-.15 , 'width', .3, 'color', odd_color, 'symbol', 'o' )
    boxplot(carr_resp, 'positions', (1:4)+.15, 'width', .3, 'color', carr_color, 'symbol', 'o' )
    x_min = 0.5; x_max = 4.5; x_units = 1;
    y_min = -100; y_max = 100; y_units = 50;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', 1:4, 'ytick', y_min:y_units:y_max , 'xticklabels', {'V1','V2','V3','V4'}, 'clipping','off', gcaOpts{:})
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:})
    axis square
    if z == 2
        text(x_min+(x_max-x_min)*.5, y_max*1.3, sprintf('change from previous image update (%s vs %s)', odd_string, carr_string), text_params{:});
    elseif z == 4
        ylabel('SOC model response (% of mean)', text_params{:});
    else
    end
    
    ax_pos(:,z) = get(gca,'position');
    if z == 1
        ax_pos(3:4,z) = ax_pos(3:4,z)*1.1;
    else
        ax_pos(3:4,z) = ax_pos(3:4,1);
    end
    if z > 3    
        ax_pos(2,z) = ax_pos(2,z)+ax_pos(4,z)*.1;
        ax_pos(1,z) = ax_pos(1,z-3);
    else
    end
    set(gca,'position',ax_pos(:,z));
    set(gca,'position',ax_pos(:,z));
end
export_fig(sprintf('%s/assessment/response_exp2.pdf',fig_folder),'-transparent');

%% now do plotting, Experiment 3
figure;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 40; fig_pos(4) = 40;
set(gcf,'pos',fig_pos);
for z = 1:8
    switch z
        case 1
            cur_set = 'carr5_odd8';
            title_str = '8 odd, 5 carr';
        case 2
            cur_set = 'carr5_odd5';
            title_str = '5 odd, 5 carr';
        case 3
            cur_set = 'carr6_odd9';
            title_str = '9 odd, 6 carr';
        case 4
            cur_set = 'carr6_odd6';
            title_str = '6 odd, 6 carr';
        case 5
            cur_set = 'carr8_odd5';
            title_str = '5 odd, 8 carr';
        case 6
            cur_set = 'carr8_odd8';
            title_str = '8 odd, 8 carr';
        case 7
            cur_set = 'carr9_odd6';
            title_str = '6 odd, 9 carr';
        case 8
            cur_set = 'carr9_odd9';
            title_str = '9 odd, 9 carr';
        otherwise
    end
    % plot size and area
    if z > 4
        subplot(4,4,z+4);
    else
        subplot(4,4,z);
    end
    hold on
    plot(carr.(cur_set).stim.diff(:,1), carr.(cur_set).stim.diff(:,3),'o', 'color', carr_color); 
    plot(odd.(cur_set).stim.diff(:,1), odd.(cur_set).stim.diff(:,3),'o', 'color', odd_color); 
    
    x_min = -1; x_max = 1; x_units = .5;
    y_min = -30; y_max = 30; y_units = 10;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', x_min:x_units:x_max, 'ytick', y_min:y_units:y_max , 'clipping', 'on', gcaOpts{:})
    if z ~= 1 && z ~= 5
        set(gca, 'yticklabels', {''});
    else
    end
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:});
    if z == 1
        text(x_min-.2*(x_max-x_min), y_max+.05*(y_max-y_min), 'A', text_params{:}, 'fontsize', f_size*2);
    elseif z == 3
        text(x_min-(x_max-x_min)*.15, y_max*1.3, sprintf('change from previous image update (%s vs %s)', odd_string, carr_string), text_params{:});
    elseif z == 5
        y_label_h = ylabel('area (% of image pixels)', text_params{:});
        set(y_label_h, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
        xlabel('dot size (% of image pixels)', text_params{:});
        text(x_min-.2*(x_max-x_min), y_max+.05*(y_max-y_min), 'C', text_params{:}, 'fontsize', f_size*2);
    else
    end
    axis square
    scatter_pos(:,z) = get(gca,'position');
    scatter_pos(1,z) = scatter_pos(1,z)-scatter_pos(3,z)*.3*(z-1);
    scatter_pos(3:4,z) = scatter_pos(3:4,z)*1.2;
    scatter_pos(2,z) = scatter_pos(2,z)-scatter_pos(4,z)*.2*(z>4);
    
    % use top positions for bottom row
    if z > 4
        scatter_pos(1,z) = scatter_pos(1,z-4);
    else
    end
    set(gca,'position', scatter_pos(:,z));
    hold off
    
    if z > 4
        subplot(4,4,z+8);
    else
        subplot(4,4,z+4);
    end
    
    hold on
    for p = 1:2
        if p == 1
            cur_color = odd_color;
            cur_dist = odd.(cur_set).stim.diff(:,5);
        else
            cur_color = carr_color;
            cur_dist = carr.(cur_set).stim.diff(:,5);
        end
        b_h = boxplot(gca, cur_dist, 'positions', p ,'symbol', 'o',  'width', .75, 'color', cur_color, 'boxstyle' ,'outline');
        set(findobj(b_h,'tag','Upper Whisker'),'linestyle', '-', 'color', cur_color);
        set(findobj(b_h,'tag','Lower Whisker'),'linestyle', '-', 'color', cur_color);
        set(findobj(b_h, 'tag', 'Outliers'), 'visible','off');
        tag_outl = findobj(b_h, 'tag', 'Outliers');
        num_outl = length(tag_outl.XData);
        if num_outl > 0
            perc_outl = num_outl/length(cur_dist)*100;
            text(p, -.4, [sprintf('%.1f', perc_outl), '%'], text_params{:}, 'fontsize', 10);
            text(p, -.8, 'outliers', text_params{:}, 'fontsize', 10);
            
        else
        end
    end
    
    x_min = 0.5; x_max = 2.5; x_units = 1;
    y_min = -1; y_max = 4; y_units = 1;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', 1:2, 'xticklabels', {'odd','carrier'}, 'ytick', y_min:y_units:y_max , 'clipping', 'on', gcaOpts{:})
    if z ~= 1 && z ~= 5
        set(gca, 'yticklabels', {''});
    else
    end
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:})
    if z == 1
        text(x_min-.25*(x_max-x_min), y_max+.05*(y_max-y_min), 'B', text_params{:}, 'fontsize', f_size*2);
    elseif z == 5
        y_label_h = ylabel('normalized distance', text_params{:});
        set(y_label_h, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
        text(x_min-.25*(x_max-x_min), y_max+.05*(y_max-y_min), 'D', text_params{:}, 'fontsize', f_size*2);
    else
    end
    axis square
    box_pos(:,z) = get(gca,'position');
    % inherit box positions from scatter positions
    box_pos([1,3,4],z) = scatter_pos([1,3,4],z);
    box_pos(2,z) = box_pos(2,z)-box_pos(4,z)*(.15 + .3*((z>4)));
    set(gca,'position', box_pos(:,z));
    hold off
end
export_fig(sprintf('%s/assessment/size_area_exp3.pdf',fig_folder),'-transparent');

%% now do plotting, Experiment 3
figure;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 40; fig_pos(4) = 20;
set(gcf,'pos',fig_pos);
for z = 1:8
    switch z
        case 1
            cur_set = 'carr5_odd8';
            title_str = '8 odd, 5 carr';
        case 2
            cur_set = 'carr5_odd5';
            title_str = '5 odd, 5 carr';
        case 3
            cur_set = 'carr6_odd9';
            title_str = '9 odd, 6 carr';
        case 4
            cur_set = 'carr6_odd6';
            title_str = '6 odd, 6 carr';
        case 5
            cur_set = 'carr8_odd5';
            title_str = '5 odd, 8 carr';
        case 6
            cur_set = 'carr8_odd8';
            title_str = '8 odd, 8 carr';
        case 7
            cur_set = 'carr9_odd6';
            title_str = '6 odd, 9 carr';
        case 8
            cur_set = 'carr9_odd9';
            title_str = '9 odd, 9 carr';
        otherwise
    end
    subplot(2,4,z);
    % adjust values
    mean_values = mean(cat(1, odd.(cur_set).resp.val, carr.(cur_set).resp.val));
    odd_resp = odd.(cur_set).resp.diff./repmat(mean_values,size(odd.(cur_set).resp.diff,1),1).*100;
    carr_resp = carr.(cur_set).resp.diff./repmat(mean_values,size(carr.(cur_set).resp.diff,1),1).*100;
    
    hold on
    b_h = boxplot(odd_resp, 'positions',  (1:4)-.15 , 'widths', 0.3, 'color', odd_color, 'symbol', 'o', 'boxstyle' ,'outline');
    set(findobj(b_h,'tag','Upper Whisker'),'linestyle', '-', 'color', odd_color);
    set(findobj(b_h,'tag','Lower Whisker'),'linestyle', '-', 'color', odd_color);
    b_h = boxplot(carr_resp, 'positions', (1:4)+.15, 'widths', 0.3, 'color', carr_color, 'symbol', 'o', 'boxstyle' ,'outline');
    set(findobj(b_h,'tag','Upper Whisker'),'linestyle', '-', 'color', carr_color);
    set(findobj(b_h,'tag','Lower Whisker'),'linestyle', '-', 'color', carr_color);
    x_min = 0.5; x_max = 4.5; x_units = 1;
    y_min = -100; y_max = 100; y_units = 50;
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    set(gca, 'xtick', 1:4, 'ytick', y_min:y_units:y_max , 'xticklabels', {'V1','V2','V3','V4'}, 'clipping','off', gcaOpts{:})
    text(mean(get(gca, 'xtick')), y_max*.95, title_str, text_params{:})
    if z == 3
        text(x_min-(x_max-x_min)*.15, y_max*1.3, sprintf('change from previous image update (%s vs %s)', odd_string, carr_string), text_params{:});
    elseif z == 5
        ylabel('SOC model response (% of mean)', text_params{:});
    else
    end
    axis square
    ax_pos(:,z) = get(gca,'position');
    if z == 1
        ax_pos(3:4,z) = ax_pos(3:4,z)*1.1;
    else
        ax_pos(3:4,z) = ax_pos(3:4,1);
    end
    if z > 4    
        ax_pos(2,z) = ax_pos(2,z)+ax_pos(4,z)*.1;
        ax_pos(1,z) = ax_pos(1,z-4);
    else
    end
    set(gca,'position',ax_pos(:,z));
    hold off
end
export_fig(sprintf('%s/assessment/response_exp3.pdf',fig_folder),'-transparent');


%%
% odd_h = findobj(gca,'color','r');
% arrayfun(@(x) uistack(x,'top'), odd_h);
% subplot(2,3,2);
% plot(s_vals_carr,d_vals_carr,'bo'); hold on
% plot(s_vals_odd, d_vals_odd,'ro'); hold on
% xlabel('size')
% ylabel('density')
% x_min = 0; x_max = .02;
% y_min = 0; y_max = .6;
% xlim([x_min,x_max]);
% ylim([y_min,y_max]);
% set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
% odd_h = findobj(gca,'color','r');
% arrayfun(@(x) uistack(x,'top'), odd_h);
% title(sprintf('values ref: %d odd: %d',numeroRef,numeroOdd));
% subplot(2,3,3);
% plot(a_vals_carr,d_vals_carr,'bo'); hold on
% plot(a_vals_odd, d_vals_odd,'ro'); hold on
% xlabel('area')
% ylabel('density')
% x_min = 0; x_max = .6;
% y_min = 0; y_max = .6;
% xlim([x_min,x_max]);
% ylim([y_min,y_max]);
% set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
% odd_h = findobj(gca,'color','r');
% arrayfun(@(x) uistack(x,'top'), odd_h);
% % plot differences
% subplot(2,3,4);
% plot(s_diff_carr,a_diff_carr,'bo'); hold on
% plot(s_diff_odd,a_diff_odd,'ro'); hold on
% xlabel('size')
% ylabel('area')
% x_min = -.01; x_max = .01;
% y_min = -.3; y_max = .3;
% xlim([x_min,x_max]);
% ylim([y_min,y_max]);
% set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
% odd_h = findobj(gca,'color','r');
% arrayfun(@(x) uistack(x,'top'), odd_h);
% subplot(2,3,5);
% plot(s_diff_carr,d_diff_carr,'bo'); hold on
% plot(s_diff_odd, d_diff_odd,'ro'); hold on
% xlabel('size')
% ylabel('density')
% x_min = -.01; x_max = .01;
% y_min = -.3; y_max = .3;
% xlim([x_min,x_max]);
% ylim([y_min,y_max]);
% set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
% odd_h = findobj(gca,'color','r');
% arrayfun(@(x) uistack(x,'top'), odd_h);
% title(sprintf('difference ref: %d odd: %d',numeroRef,numeroOdd));
% subplot(2,3,6);
% plot(a_diff_carr,d_diff_carr,'bo'); hold on
% plot(a_diff_odd, d_diff_odd,'ro'); hold on
% xlabel('area')
% ylabel('density')
% x_min = -.3; x_max = .3;
% y_min = -.3; y_max = .3;
% xlim([x_min,x_max]);
% ylim([y_min,y_max]);
% set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
% odd_h = findobj(gca,'color','r');
% arrayfun(@(x) uistack(x,'top'), odd_h);
%             
                
            
            
%             gcf;
%             subplot(2,3,1);
%             plot(s_vals_carr,a_vals_carr,'bo'); hold on
%             plot(s_vals_odd,a_vals_odd,'ro'); hold on
%             xlabel('size')
%             ylabel('area')
%             x_min = 0; x_max = .02;
%             y_min = 0; y_max = .6;
%             xlim([x_min,x_max]);
%             ylim([y_min,y_max]);
%             set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
%             odd_h = findobj(gca,'color','r');
%             arrayfun(@(x) uistack(x,'top'), odd_h);
%             subplot(2,3,2);
%             plot(s_vals_carr,d_vals_carr,'bo'); hold on
%             plot(s_vals_odd, d_vals_odd,'ro'); hold on
%             xlabel('size')
%             ylabel('density')
%             x_min = 0; x_max = .02;
%             y_min = 0; y_max = .6;
%             xlim([x_min,x_max]);
%             ylim([y_min,y_max]);
%             set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
%             odd_h = findobj(gca,'color','r');
%             arrayfun(@(x) uistack(x,'top'), odd_h);
%             title(sprintf('values ref: %d odd: %d',numeroRef,numeroOdd));
%             subplot(2,3,3);
%             plot(a_vals_carr,d_vals_carr,'bo'); hold on
%             plot(a_vals_odd, d_vals_odd,'ro'); hold on
%             xlabel('area')
%             ylabel('density')
%             x_min = 0; x_max = .6;
%             y_min = 0; y_max = .6;
%             xlim([x_min,x_max]);
%             ylim([y_min,y_max]);
%             set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
%             odd_h = findobj(gca,'color','r');
%             arrayfun(@(x) uistack(x,'top'), odd_h);
%             % plot differences
%             subplot(2,3,4);
%             plot(s_diff_carr,a_diff_carr,'bo'); hold on
%             plot(s_diff_odd,a_diff_odd,'ro'); hold on
%             xlabel('size')
%             ylabel('area')
%             x_min = -.01; x_max = .01;
%             y_min = -.3; y_max = .3;
%             xlim([x_min,x_max]);
%             ylim([y_min,y_max]);
%             set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
%             odd_h = findobj(gca,'color','r');
%             arrayfun(@(x) uistack(x,'top'), odd_h);
%             subplot(2,3,5);
%             plot(s_diff_carr,d_diff_carr,'bo'); hold on
%             plot(s_diff_odd, d_diff_odd,'ro'); hold on
%             xlabel('size')
%             ylabel('density')
%             x_min = -.01; x_max = .01;
%             y_min = -.3; y_max = .3;
%             xlim([x_min,x_max]);
%             ylim([y_min,y_max]);
%             set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
%             odd_h = findobj(gca,'color','r');
%             arrayfun(@(x) uistack(x,'top'), odd_h);
%             title(sprintf('difference ref: %d odd: %d',numeroRef,numeroOdd));
%             subplot(2,3,6);
%             plot(a_diff_carr,d_diff_carr,'bo'); hold on
%             plot(a_diff_odd, d_diff_odd,'ro'); hold on
%             xlabel('area')
%             ylabel('density')
%             x_min = -.3; x_max = .3;
%             y_min = -.3; y_max = .3;
%             xlim([x_min,x_max]);
%             ylim([y_min,y_max]);
%             set(gca,'xtick', linspace(x_min,x_max,5), 'ytick',linspace(y_min,y_max,5), 'tickdir', 'out')
%             odd_h = findobj(gca,'color','r');
%             arrayfun(@(x) uistack(x,'top'), odd_h);
% 
%             set(gcf,'units','centimeters');
%             fig_pos = get(gcf,'pos');
%             fig_pos(3) = 40; fig_pos(4) = 20;
%             set(gcf,'pos',fig_pos);
%             fileName = [];
%             for i=1:length(imgSeq)
%                 if imgSeq(i) ~= 0
%                     if isempty(fileName)
%                         fileName = sprintf('~/Desktop/numeroDemo_%s',datestr(now,'yyyymmddHHMMSS'));
%                         imwrite(img(:,:,1,imgSeq(i)),[fileName,'.gif'],'gif','LoopCount',Inf,'DelayTime',1/nFrameCycleRef);
%                     else
%                         imwrite(img(:,:,1,imgSeq(i)),[fileName,'.gif'],'WriteMode','append','DelayTime',1/nFrameCycleRef);
%                     end
%                 end
%             end
%         
%         export_fig(sprintf('ref%s_odd%s', carr_vals{c}, odd_vals{o}),'-transparent','-pdf');
%         %savefig(sprintf('ref%s_odd%s.fig', carr_vals{c}, odd_vals{o}))
%         close gcf;