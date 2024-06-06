clear; close all;

%% Correlation
%%
load('data_VIP_CC\rev_dataname.mat')

rev_trial_anova_win = [232-142 264-160 242-141 214-141 273-153 260-144 197-141 226-141 254-141 236-143 270-140 186-140];

%% sliding
n_num = num2str(length(beh_file));
win_size = 1.5;
event_thr_set = [0.015 0.025 0.020];
hz_cut = 2;
step_size = 1;
if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
else; win_size_str = ''; end
if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
else; hz_cut_str = ''; end
if step_size~=10; step_str = [' ',num2str(step_size),'step'];
else; step_str = ''; end

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome sliding revonset'...
    win_size_str ' n' n_num hz_cut_str step_str ' 230627\SRC_set.mat'])
load(src_filename_set{1})
x_var_num = size(SRC,2);
step_num=  size(SRC,3);
ses_num= size(src_filename_set,2);
T_value_unify_sliding = cell(1,ses_num);
for ifile = 1:ses_num
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name = ['1234' name];
    name(1:7) = 't_value';
    load(fullfile(pathstr,[name ext]))
    T_value_unify_sliding{ifile} = T_value;
end

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome term sliding revonset'...
    win_size_str ' n' n_num hz_cut_str step_str ' 230627\SRC_set.mat'])
load(src_filename_set{1})
T_value_unify_sliding_rw = cell(1,ses_num);
for ifile = 1:size(src_filename_set,2)
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name = ['1234' name];
    name(1:7) = 't_value';
    load(fullfile(pathstr,[name ext]))
    T_value_unify_sliding_rw{ifile} = T_value_rwd;
end

%% plot & index
rise_stepxmean_r = zeros(1,ses_num);

delete_id = [];

smooth_trial = 40;

%%
for ifile = 1:ses_num
    if isempty(T_value_unify_sliding{ifile}); delete_id = [delete_id ifile]; continue; end
    %% LDI
    load(beh_file{ifile})
    
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;
    
    % lick time
    if length(lickTime)>1
        lick_diff = diff(lickTime, 1);
        short_lick_diff = find(lick_diff(:,1)<=50)+1;
        lickTime_50ms = lickTime;
        lickTime_50ms(short_lick_diff,:)=[];
    elseif length(lickTime)==1
        lickTime_50ms = lickTime;
    end
    
    lick_timings_50ms = cell(nTrial,1);
    delay_timings = stateTime(:,3:4)-stateTime(:,1)+1;
    delay_licknum_50ms = zeros(size(lickNum));
    for j = 1:nTrial
        jtrial_ind = find(lickTime_50ms(:,2)==j);
        lick_timings_50ms{j} = lickTime_50ms(jtrial_ind,1) - stateTime(j,1)+1;
        delay_licknum_50ms(j) = length(find((lick_timings_50ms{j}>=delay_timings(j,1)) & (lick_timings_50ms{j}<delay_timings(j,2))));
    end
    
    outcomeID =  outcomeIdentity(1,:);
    outcomeProb = outcomeProbability(1,:);
    pr_cue = find(outcomeID==3&outcomeProb==75);
    rp_cue = find(outcomeID==2&outcomeProb==75);
    
    pr_trials = odorCue == (pr_cue-1);
    rp_trials = odorCue == (rp_cue-1);
        
    % plot
    rev_step = min_step_num_left+1;
    x_plot = ((win_size_sliding/2):step_size:(step_num*step_size+(win_size_sliding/2)-step_size)); % first trial will be ingored by the next line
    x_plot = x_plot-x_plot(rev_step);
    
    dep_win = 100;
    
    first_window = (idenrev_trial-ceil(win_size_sliding/2)+1):(idenrev_trial+floor(win_size_sliding/2));
           
    temp_window = first_window;
    for istep = 1:420/step_size
        temp_window = temp_window - step_size;
        if temp_window(1)<=0
            temp_window = temp_window + step_size;
            temp_left_step_num = istep-1;
            temp_first_trial = temp_window(1);
            break
        end
    end
    temp_first_trial = temp_first_trial + step_size*(temp_left_step_num-min_step_num_left);
    
    temp_dep = zeros(1,min_step_num);
    istep = 1;
    temp_win = (1:dep_win)-step_size+(temp_first_trial-1)+(win_size_sliding-dep_win)/2;
    for istep = 1:min_step_num
        temp_win = temp_win + step_size;
        temp_lick = delay_licknum_50ms(temp_win);
        temp_pr = pr_trials(temp_win);
        temp_rp = rp_trials(temp_win);
        
        temp_pr_lick = mean(temp_lick(temp_pr));
        temp_rp_lick = mean(temp_lick(temp_rp));
        
        temp_dep(istep) = (temp_rp_lick-temp_pr_lick)/(temp_rp_lick+temp_pr_lick);
    end
    
    %%
    if size(T_value_unify{iepoch,ifile},1)<3; delete_id = [delete_id ifile]; continue; end
    
    %% sliding
    temp_T_value_sliding = T_value_unify_sliding{ifile};
    data_draw_csrw_sliding = squeeze(temp_T_value_sliding(:,1,:));
    
    temp_T_value_rw_sliding = T_value_unify_sliding_rw{ifile};
    data_draw_rw_sliding = squeeze(temp_T_value_rw_sliding(:,3,:));
    
    %%
    r_set = zeros(1,step_num);
    p_set = zeros(1,step_num);
    
    for istep = 1:step_num
        temp_csrw = data_draw_csrw_sliding(:,istep);
        temp_rw = data_draw_rw_sliding(:,istep);
        [r,p] = corrcoef(temp_csrw,temp_rw);
        r_set(istep) = r(1,2);
        p_set(istep) = p(1,2);
    end
    
    %% index
    r_set = movmean(r_set,smooth_trial);

    tmp_plot_id = true(1,ses_num);
    tmp_plot_id([7 9]) = false;
    
    [~,temp_min_id_xmean] = min(r_set(rev_step:round(mean(rev_trial_anova_win(tmp_plot_id)))-1+rev_step));
    try
        [temp_max_after_x,temp_max_id_x] = max(r_set(rev_trial_anova_win(ifile)-1+rev_step:rev_trial_anova_win(ifile)-1+rev_step+100/step_size));
    catch
        [temp_max_after_x,temp_max_id_x] = max(r_set(rev_trial_anova_win(ifile)-1+rev_step:end));
    end
    rise_stepxmean_r(ifile) = rev_trial_anova_win(ifile)-1-temp_min_id_xmean+temp_max_id_x;
    
end

%%
plot_id = true(1,ses_num);
plot_id(delete_id) = false;

f_beh = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 9]);
temp_plot = rise_stepxmean_r;
nan_id = isnan(temp_plot);

temp_plot = temp_plot(plot_id&~nan_id);
temp_anova = rev_trial_anova_win(plot_id&~nan_id);

[r,p] = corrcoef(temp_plot,temp_anova);

ax = axes;
hold on
ax.FontSize = 16;
label_fontsize = 20;

scatter(temp_plot,temp_anova,'filled','sizedata',100,'markerfacecolor',[0.6350 0.0780 0.1840],'markeredgecolor',[0.6350 0.0780 0.1840])
ylim([0 150])
yticks(0:50:200)
xlim([0 200])

ylabel('Trials to criterion','fontsize',label_fontsize)
xlabel('RPE persistence','fontsize',label_fontsize)
if (min(temp_plot) * max(temp_plot)) < 0
    yline1 = line([0 0],[0 ceil(max(temp_anova)/100)*100]);
    yline1.Color = [0.5 0.5 0.5]; yline1.LineStyle = '--';
end
fit1 = polyfit(temp_plot,temp_anova,1);
plot([min(temp_plot) max(temp_plot)],polyval(fit1,[min(temp_plot) max(temp_plot)]),'Color','k','linewidth',1.5)

str1 = {['{\itr} = ',num2str(r(1,2),'%.3f')],['{\itp} = ',num2str(p(1,2),'%.3f')]};
annotation('textbox',[0.6 0.2 0.3 0.2],'String',str1,'FitBoxToText','on','edgecolor','none','fontsize',14)
if p(1,2)<0.05; star_str = '*'; if p(1,2)<0.01; star_str = '**'; if p(1,2)<0.001; star_str = '***'; end; end; else; star_str = ''; end
annotation('textbox',[0.84 0.16 0.2 0.2],'string',star_str,'FitBoxToText','on','edgecolor','none','fontsize',30)
