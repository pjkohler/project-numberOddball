clear all; close all;
load('xDiva_params.mat');

% add paths
addpath(genpath('/Users/kohler/code/git/mrC/tools/colors'))
addpath(genpath('/Users/kohler/code/git/export_fig'))
addpath(genpath('./socmodel'))
addpath('./minbound_suite')
addpath(genpath('/Users/kohler/code/git/xDiva'))

fig_folder = '/Users/kohler/Dropbox/WRITING/Articles/2019_KohlerNumerositySSVEP/figures';

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

load_values = true;

if ~load_values
    % load bpfilter
    load('stimuli.mat', 'bpfilter');
    global img_out; global seq_out;
    img_array = [];
    for o = 1:length(odd_vals)
        parameters{3}{11} = odd_vals{o};
        for c = 1:length(carr_vals)
            parameters{4}{11} = carr_vals{c};
            close all;
            cur_set = (sprintf('carr%s_odd%s', carr_vals{c}, odd_vals{o}));
            odd.(cur_set).size.val = []; odd.(cur_set).size.diff = [];
            odd.(cur_set).density.val = []; odd.(cur_set).density.diff = [];
            odd.(cur_set).area.val = []; odd.(cur_set).area.diff = [];
            carr.(cur_set).size.val = []; carr.(cur_set).size.diff = [];
            carr.(cur_set).density.val = []; carr.(cur_set).density.diff = [];
            carr.(cur_set).area.val = []; carr.(cur_set).area.diff = [];
            all_images = [];
            for z = 1:10
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

            % loop over images
            s_vals = []; a_vals = [];
            for i = 1:size(all_images,3)
                temp_img = all_images(:,:,i);
                mask = zeros(size(temp_img));
                mask(temp_img ~= mode(temp_img(:))) = 1;
                c_hull = bwconvhull(mask);
                strct = regionprops(mask,'PixelList'); % coordinates of convex hull
                test = strct.PixelList;
                [bound_c,bound_r] = minboundcircle(test(:,1),test(:,2));                
                dots = bwlabel(mask,4);
                s_vals(i) = numel(mask);
                for q = 1:(length(unique(dots))-1)
                    temp_size = length(find(dots == q)); % just grab smallest dot (dot should be same size, but can overlap)
                    % normalize by array size
                    temp_size = temp_size./numel(mask) * 100;
                    if s_vals(i) > temp_size
                        s_vals(i) = temp_size;
                    else
                    end
                end
                s_vals(i) = temp_size(1);
                a_vals(i) = length(find(c_hull == 1));
                % use 'native' estimate:
                % a_vals(i) = ceil(bound_r^2*pi);
                % normalize by array size 
                a_vals(i) = a_vals(i)./numel(mask) * 100;
            end

            a_diff = diff(a_vals);
            s_diff = diff(s_vals);

            num_vals = zeros(size(s_vals));
            num_vals(val_idx == 1) = odd_vals{o};
            num_vals(val_idx == 0) = carr_vals{c};
            d_vals = (s_vals.*num_vals)./a_vals;
            d_diff = diff(d_vals);

            % area
            odd.(cur_set).area.val = a_vals(val_idx==1);
            odd.(cur_set).area.diff = a_diff(diff_idx==1);
            carr.(cur_set).area.val = a_vals(val_idx==0);
            carr.(cur_set).area.diff = a_diff(diff_idx==0);
            % size
            odd.(cur_set).size.val = s_vals(val_idx==1);
            odd.(cur_set).size.diff = s_diff(diff_idx==1);
            carr.(cur_set).size.val = s_vals(val_idx==0);
            carr.(cur_set).size.diff = s_diff(diff_idx==0);
            % density
            odd.(cur_set).density.val = d_vals(val_idx==1);
            odd.(cur_set).density.diff = d_diff(diff_idx==1);
            carr.(cur_set).density.val = d_vals(val_idx==0);
            carr.(cur_set).density.diff = d_diff(diff_idx==0);

            if ismember(cur_set, model_sets)
                sd_param = [0.9308, 1.0738, 1.4671, 2.1242];
                n_param = [0.1814, 0.1285, 0.1195, 0.1152]; 
                c_param = [0.9276, 0.9928, 0.9941, 0.9472];
                cache = [];
                
                % downsample image to match Kendrick's images
                kay_res = 256;
                kay_fov = 12.7;
                our_fov = 10;
                %new_size = kay_res/kay_fov*our_fov; 
                % actual size 201.5748, use 200;
                new_size = 200;
                new_images = imresize(all_images, [new_size, new_size]);
                bg_val = mode(new_images(:));
                filt_images = padarray(new_images, [new_size/5, new_size/5], bg_val, 'both'); 
                filt_images = arrayfun(@(x) conv2(filt_images(:,:,x), bpfilter), 1:size(filt_images, 3),'uni',false);
                crop_val =  (size(filt_images{1},1)-new_size)/2;
                filt_images = cellfun(@(x) x(crop_val+1:crop_val+new_size,crop_val+1:crop_val+new_size), filt_images, 'uni', false);
                filt_images = cell2mat(reshape(filt_images, [1, 1, size(filt_images,2)]));
                for v = 1:4    
                    [temp_resp,cache] = socmodel(filt_images,new_size,[],1.2,{new_size/4 -1 1 8 2 .01 2 0}, ...
                                        1,.5, sd_param(v), 1/sd_param(v), n_param(v), c_param(v), cache);                      
                    all_resp(:,v) = squeeze(mean(mean(temp_resp,1),2));
                    if strcmp('carr6_odd6',cur_set)
                        if v == 1
                            example_image = new_images(:,:,1);
                        else
                        end
                        example_resp(:,:,v) = temp_resp(:,:,1);
                    else
                    end
                        
                end
                odd.(cur_set).resp.val = all_resp(val_idx==1, :);
                carr.(cur_set).resp.val = all_resp(val_idx==0, :);
                diff_resp = diff(all_resp);
                odd.(cur_set).resp.diff = diff_resp(diff_idx==1, :);
                carr.(cur_set).resp.diff = diff_resp(diff_idx==0, :);   
            else
            end
            clear all_images all_resp;
        end
    end
    save('assessment_data.mat','odd','carr','example_*', '-v7.3');
else
    load('assessment_data.mat','odd','carr', 'example_*');
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
imagesc(example_resp(:,:,1),[0 1.5]); axis image tight; colormap(gray); title('V1 response', text_params{:}); axis off;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 20; fig_pos(4) = 20;
set(gcf,'pos',fig_pos);
export_fig(sprintf('%s/examples/response_example.pdf',fig_folder),'-transparent', '-nocrop');

mask = zeros(size(example_image));
mask(example_image ~= mode(example_image(:))) = 1;
c_hull = bwconvhull(mask);            
dots = bwlabel(mask,4);
example_out = example_image;
example_out(~c_hull) = 0.5;
figure;
imagesc(example_out,[0 254]); axis image tight; colormap(gray); title('stimulus (with convex hull)', text_params{:}); axis off;
set(gcf,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 20; fig_pos(4) = 20;
set(gcf,'pos',fig_pos);
export_fig(sprintf('%s/examples/stimulus_example.pdf',fig_folder),'-transparent', '-nocrop');

for c = 1:6
    if c == 6
        imshow(odd.carr8_odd5.images(:,:,1));
    else
        imshow(carr.carr8_odd5.images(:,:,c));
    end
    axis off
    set(gcf,'units','centimeters');
    fig_pos = get(gcf,'pos');
    fig_pos(3) = 20; fig_pos(4) = 20;
    set(gcf,'pos',fig_pos);
    export_fig(sprintf('%s/examples/stimulus%d.png',fig_folder, c),'-png','-opengl','-m5','-transparent',gcf);
end

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
    plot(carr.(cur_set).size.diff, carr.(cur_set).area.diff,'o', 'color', carr_color); 
    plot(odd.(cur_set).size.diff, odd.(cur_set).area.diff,'o', 'color', odd_color); 
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
    z_size = zscore([carr.(cur_set).size.diff, odd.(cur_set).size.diff]);
    z_area = zscore([carr.(cur_set).area.diff, odd.(cur_set).area.diff]);
    z_dist = sqrt(z_size.^2 + z_area.^2);
    carr_dist = z_dist(1:length(carr.(cur_set).size.diff));
    odd_dist = z_dist(length(carr.(cur_set).size.diff)+1:end);
    
    hold on
    for p = 1:2
        if p == 1
            cur_color = odd_color;
            cur_dist = odd_dist;
        else
            cur_color = carr_color;
            cur_dist = carr_dist;
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
    plot(carr.(cur_set).size.diff, carr.(cur_set).area.diff,'o', 'color', carr_color); 
    plot(odd.(cur_set).size.diff, odd.(cur_set).area.diff,'o', 'color', odd_color); 
    
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
    z_size = zscore([carr.(cur_set).size.diff, odd.(cur_set).size.diff]);
    z_area = zscore([carr.(cur_set).area.diff, odd.(cur_set).area.diff]);
    z_dist = sqrt(z_size.^2 + z_area.^2);
    carr_dist = z_dist(1:length(carr.(cur_set).size.diff));
    odd_dist = z_dist(length(carr.(cur_set).size.diff)+1:end);
    
    hold on
    for p = 1:2
        if p == 1
            cur_color = odd_color;
            cur_dist = odd_dist;
        else
            cur_color = carr_color;
            cur_dist = carr_dist;
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