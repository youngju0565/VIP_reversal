clear; close all;

%% Trace with LDI
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
    name(1:7) = 't_vlaue';
    load(fullfile(pathstr,[name ext]))
    T_value_unify_sliding{ifile} = T_value;
end

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome term sliding revonset'...
    win_size_str ' n' n_num hz_cut_str step_str ' 230627\SRC_set.mat'])
load(src_filename_set{1})
T_value_unify_sliding_rw = cell(1,ses_num);
T_value_unify_sliding_pn = cell(1,ses_num);
for ifile = 1:size(src_filename_set,2)
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name = ['1234' name];
    name(1:7) = 't_vlaue';
    load(fullfile(pathstr,[name ext]))
    T_value_unify_sliding_rw{ifile} = T_value_rwd;
    T_value_unify_sliding_pn{ifile} = T_value_pun;
end

rev_step = min_step_num_left+1;
x_plot = ((win_size_sliding/2):step_size:(step_num*step_size+(win_size_sliding/2)-step_size)); % first trial will be ingored by the next line
x_plot = x_plot-x_plot(rev_step);

%%
r_set_set = cell(1,ses_num);
dep_set = cell(1,ses_num);

r_max_x_sat = nan(1,ses_num);
LDI_min_sat2 = nan(1,ses_num);

delete_id = [];

smooth_trial = 40;

dep_win = 100; 

%%
for ifile = 1:ses_num
    if isempty(T_value_unify_sliding{ifile}); delete_id = [delete_id ifile]; continue; end
    if size(T_value_unify_sliding{ifile},1)<3; delete_id = [delete_id ifile]; continue; end
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
    beh_file_name = beh_file{ifile};
    beh_file_name_split = strsplit(beh_file_name,'_');
    animal_name = beh_file_name_split{1};
    
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
    temp_win = ((1:dep_win)-step_size+(temp_first_trial-1))+(win_size_sliding-dep_win)/2;
    for istep = 1:min_step_num
        temp_win = temp_win + step_size;
        temp_lick = delay_licknum_50ms(temp_win);
        temp_pr = pr_trials(temp_win);
        temp_rp = rp_trials(temp_win);
        
        temp_pr_lick = mean(temp_lick(temp_pr));
        temp_rp_lick = mean(temp_lick(temp_rp));
        
        temp_dep(istep) = (temp_rp_lick-temp_pr_lick)/(temp_rp_lick+temp_pr_lick);
    end
    
    dep_set{ifile} = temp_dep;
    
    dep_smooth = movmean(temp_dep,smooth_trial);

    %% sliding
    temp_T_value_sliding = T_value_unify_sliding{ifile};
    data_draw_csrw_sliding = squeeze(temp_T_value_sliding(:,1,:));
    
    temp_T_value_rw_sliding = T_value_unify_sliding_rw{ifile};
    data_draw_rw_sliding = squeeze(temp_T_value_rw_sliding(:,3,:));
    
    r_set = zeros(1,step_num);
    p_set = zeros(1,step_num);
    for istep = 1:step_num
        temp_csrw = data_draw_csrw_sliding(:,istep);
        temp_rw = data_draw_rw_sliding(:,istep);
        [r,p] = corrcoef(temp_csrw,temp_rw);
        r_set(istep) = r(1,2);
        p_set(istep) = p(1,2);
    end
    
    r_set_set{ifile} = r_set;
    
    r_set_smooth = movmean(r_set,smooth_trial);
    
    %%
    [temp_min_x2,temp_min_id_x2] = min(r_set_smooth(rev_step:100/step_size-1+rev_step)); % new idea
    try
        [temp_max_after_x,temp_max_id_x] = max(r_set_smooth(rev_trial_anova_win(ifile)-1+rev_step:rev_trial_anova_win(ifile)-1+rev_step+100/step_size));
    catch
        [temp_max_after_x,temp_max_id_x] = max(r_set_smooth(rev_trial_anova_win(ifile)-1+rev_step:end));
    end
    r_sat_point_x = (temp_max_after_x-temp_min_x2)*0.95+temp_min_x2;
    r_sat_id = find(r_set_smooth(rev_step+temp_min_id_x2-1:rev_trial_anova_win(ifile)-1+rev_step+temp_max_id_x-1)>r_sat_point_x,1);
    temp_x = x_plot(rev_step+temp_min_id_x2-1:rev_trial_anova_win(ifile)-1+rev_step+temp_max_id_x-1);
    r_max_x_sat(ifile) = x_plot(x_plot==temp_x(r_sat_id));
    
    % max-min
    temp_ldi_sat = max(temp_dep)+(min(temp_dep)-max(temp_dep))*0.95;
    LDI_min_sat2(ifile) = x_plot(find(temp_dep<temp_ldi_sat,1)); % trial num
end

%%
fill_x = [x_plot fliplr(x_plot)];
plot_id2 = cellfun(@(x) ~isempty(x), r_set_set);

smooth_trial_ai = 40;

f_mean_ai = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 9]);
set(f_mean_ai,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
ax = axes;
set(gca,'xcolor','k','ycolor','k')
hold on
ax.FontSize = 16;

yyaxis left
temp_data = cell2mat(r_set_set');
smooth_r = movmean(temp_data,smooth_trial_ai,2);
temp_mean = mean(smooth_r);
temp_sem = std(smooth_r)/sqrt(size(smooth_r,1));
fill_y = [temp_mean+temp_sem fliplr(temp_mean-temp_sem)];
fill_temp = fill(fill_x,fill_y,'k');
alpha(fill_temp,0.15)
fill_temp.EdgeColor = 'none';
plot(x_plot, temp_mean,'linestyle','-','color','k','marker','.','markerfacecolor','k')
ylim([-0.45 0.45])
yticks(-0.4:0.2:0.4)
ylabel('\itr','fontsize',20)

yyaxis right
temp_data = cell2mat(dep_set');
smooth_l = movmean(temp_data,smooth_trial_ai,2);
temp_mean = mean(smooth_l);
temp_sem = std(smooth_l)/sqrt(size(smooth_l,1));
fill_y = [temp_mean+temp_sem fliplr(temp_mean-temp_sem)];
fill_temp = fill(fill_x,fill_y,'b');
alpha(fill_temp,0.2)
fill_temp.EdgeColor = 'none';
plot(x_plot, temp_mean,'linestyle','-','color','b','marker','.','markerfacecolor','b')
ylim([-1 1])
yticks(-1:0.5:1)
ylabel('LDI','fontsize',20)

plot([min(x_plot) max(x_plot)],[0 0],'--','color',[0.5 0.5 0.5])
plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5])
plot([mean(rev_trial_anova_win(plot_id2)) mean(rev_trial_anova_win(plot_id2))],[-2 2],'--','color',[0.5 0.5 0.5])
xlim([min(x_plot) max([max(x_plot) mean(rev_trial_anova_win(plot_id2))])])
xlabel('Trial since reversal','fontsize',20)

plot(nanmean(r_max_x_sat),1,'color','k','marker','v','markerfacecolor','k','markeredgecolor','k','linestyle','none','markersize',10)
plot(nanmean(LDI_min_sat2),1,'color','b','marker','v','markerfacecolor','b','markeredgecolor','b','linestyle','none','markersize',10)
