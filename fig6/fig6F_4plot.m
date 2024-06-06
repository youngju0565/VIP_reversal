clear; close all;

%%
n_num = num2str(12);
win_size = 1.5;
event_thr_set = [0.015 0.025 0.020];
hz_cut = 2;
win_size_sliding=100;
step_size = 1; % 10
if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
else; win_size_str = ''; end
if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
else; hz_cut_str = ''; end
if win_size_sliding~=100; trial_str = [' ',num2str(win_size_sliding),'trials'];
else; trial_str = ''; end
if step_size~=10; step_str = [' ',num2str(step_size),'step'];
else; step_str = ''; end

%%
fit_name = '3var Ras';
fit_var_num = 3;

%%
anccr_varname_set_nofix = {'PR (CSrw Rw)','PR (CSrw Pn)','PR (CSpn Rw)','PR (CSpn Pn)'};

%%
load('data_VIP_CC\rev_dataname.mat')

dep_win = 100;
smooth_trial = 40;

%%
load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Mat_code\ANCCR\fit VIP ' fit_name '\reversal_set regression outcome ANCCR ' num2str(fit_var_num) 'var sliding revonset' ...
    win_size_str ' n' n_num hz_cut_str trial_str step_str ' 230627\SRC_set.mat'])
load(src_filename_set{1})
x_var_num = 14;
SRC_model = SRC(:,1:x_var_num,:);
step_num=  size(SRC_model,3);

T_value_unify = cell(x_var_num,step_num);
dep_set = cell(1,size(src_filename_set,2));

t_sat = nan(1,size(src_filename_set,2));
LDI_sat = nan(1,size(src_filename_set,2));

rev_step = min_step_num_left+1;
x_plot = ((win_size_sliding/2):step_size:(step_num*step_size+(win_size_sliding/2)-step_size));
x_plot = x_plot-x_plot(rev_step);

rev_trial_anova_win = [232-142 264-160 242-141 214-141 273-153 260-144 197-141 226-141 254-141 236-143 270-140 186-140];
success_step_set = rev_trial_anova_win-1;

success_trial = mean(rev_trial_anova_win);

%%
for ifile = 1:size(src_filename_set,2)
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name = ['1234' name];
    name(1:7) = 't_value';
    load(fullfile(pathstr,[name ext]))
    for istep = 1:step_num
        for ivar = 1:x_var_num
            T_value_unify{ivar,istep}(ifile) = mean(abs(T_value(:,ivar,istep)));
        end
    end
    
    %% LDI
    load(beh_file{ifile})
    
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;
    
    % lick time
    if length(lickTime)>1
        lick_diff = diff(lickTime, 1);
        short_lick_diff = find(lick_diff(:,1)<=50)+1;
        lickTime_50ms = lickTime;
        lickTime_50ms(short_lick_diff,:)=[]; % 50ms 이내 간격이면 뒤의것 제외하기
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
    
    %% saturation
    temp_data_retrovar_sliding = movmean(permute(mean(abs(T_value(:,11,:)),1),[2 3 1]),smooth_trial);
    tmp_plot_id = true(size(success_step_set));
    tmp_plot_id(6)=false;
    
    [temp_max_mean,temp_max_after_id_xmean] = max(temp_data_retrovar_sliding(rev_step:round(mean(success_step_set(tmp_plot_id)))+rev_step));
    try
        [temp_min_x,temp_min_id_x] = min(temp_data_retrovar_sliding(success_step_set(ifile)+rev_step:success_step_set(ifile)+rev_step+100/step_size));
    catch
        [temp_min_x,temp_min_id_x] = min(temp_data_retrovar_sliding(success_step_set(ifile)+rev_step:end));
    end
    t_sat_point_x = (temp_min_x-temp_max_mean)*0.95+temp_max_mean;
    t_sat_id = find(temp_data_retrovar_sliding(rev_step+temp_max_after_id_xmean-1:success_step_set(ifile)+rev_step+temp_min_id_x-1)<t_sat_point_x,1);
    temp_x = x_plot(rev_step+temp_max_after_id_xmean-1:success_step_set(ifile)+rev_step+temp_min_id_x-1);
    t_sat(ifile) = x_plot(x_plot==temp_x(t_sat_id));
    
    % max-min
    temp_ldi_sat = max(temp_dep)+(min(temp_dep)-max(temp_dep))*0.95;
    LDI_sat(ifile) = x_plot(find(temp_dep<temp_ldi_sat,1)); % trial num
end

%%
var_name_set = [{'CS_r_w','CS_p_n','Rw','Pn','Lick','Speed',...
    'CS_r_w(t-1)','CS_p_n(t-1)','Rw(t-1)','Pn(t-1)'}, anccr_varname_set_nofix];

color_set = {[99 103 161]/255 [199 87 107]/255 [30/255 45/255 232/255] [255/255 30/255 70/255] ...
    [0.5 0.5 0.5] [0 0 0] ...
    [99 103 161]/255 [199 87 107]/255 [30/255 45/255 232/255] [255/255 30/255 70/255] ...
    [30 45 232]/255 [255 30 70]/255 [30 45 232]/255 [255 30 70]/255 ...
    [30 45 232]/255 [255 30 70]/255 [30 45 232]/255 [255 30 70]/255 ...
    [30 45 232]/255 [255 30 70]/255 [30 45 232]/255 [255 30 70]/255};

f_var_abs = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);
set(f_var_abs,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

%%
ifig = 11;
temp_data_set = abs(cell2mat(T_value_unify(ifig,:)')');
temp_data_set = movmean(temp_data_set,smooth_trial,2);
data_set_mean = mean(temp_data_set);
data_set_sem = std(temp_data_set)/sqrt(size(temp_data_set,1));

clf
ax = axes;
ax.FontSize = 16;
set(gca,'xcolor','k','ycolor','k')
hold on

yyaxis left
if size(temp_data_set,1)>1
    fill_x = [x_plot fliplr(x_plot)];
    fill_y = [data_set_mean+data_set_sem fliplr(data_set_mean-data_set_sem)];
    filled = fill(fill_x,fill_y,[224 240 219]/255);
    alpha(filled,0.4)
    filled.EdgeColor = 'none';
end
plot(x_plot, data_set_mean,'-','color',[13 137 38]/255,'markersize',8,'linewidth',1.5)
ylim([0.65 1.25])
yticks(0:0.1:2)
plot([min(x_plot) max(x_plot)],[0 0],'--','color',[0.5 0.5 0.5])
plot([x_plot(rev_step) x_plot(rev_step)],[-100 100],'--','color',[0.5 0.5 0.5])
plot([success_trial success_trial],[-100 100],'--','color',[0.5 0.5 0.5])
xlim([min(x_plot) max([success_trial max(x_plot)])])
ylabel('|{\itt}-value|','fontsize',20)

yyaxis right
temp_data = cell2mat(dep_set');
smooth_l = movmean(temp_data,smooth_trial,2);
temp_mean = mean(smooth_l);
temp_sem = std(smooth_l)/sqrt(size(smooth_l,1));
fill_y = [temp_mean+temp_sem fliplr(temp_mean-temp_sem)];
fill_temp = fill(fill_x,fill_y,[217 179 255]/255);
alpha(fill_temp,0.2)
fill_temp.EdgeColor = 'none';
plot(x_plot, temp_mean,'linestyle','-','color',[102 0 204]/255,'marker','.','markerfacecolor',[102 0 204]/255)
ylim([-1 1])
yticks(-1:0.5:1)
ylabel('LDI','fontsize',20)

plot([nanmean(t_sat) nanmean(t_sat)],[-2 0.85],'color',[13 137 38]/255,'marker','o','markerfacecolor',[13 137 38]/255,'markeredgecolor',[13 137 38]/255,'linestyle','--','markersize',5)
plot([nanmean(LDI_sat) nanmean(LDI_sat)],[-2 0.85],'color',[102 0 204]/255,'marker','o','markerfacecolor',[102 0 204]/255,'markeredgecolor',[102 0 204]/255,'linestyle','--','markersize',5)

t_sem = nanstd(t_sat)/sqrt(sum(~isnan(t_sat)));
LDI_sem = nanstd(LDI_sat)/sqrt(sum(~isnan(LDI_sat)));

plot(nanmean(t_sat)+[-t_sem/2 t_sem/2],[1 1]*0.85,'-','color',[13 137 38]/255)
plot(nanmean(LDI_sat)+[-LDI_sem/2 LDI_sem/2],[1 1]*0.85,'-','color',[102 0 204]/255)

xlabel('Trial since reversal','fontsize',20)
title('PR[CSrw\leftarrowRw]','fontsize',20)
