clear; close all;

%% RP reversal: Behavior
%%
load('data_VIP_CC\rev_dataname.mat')

ses_num = length(beh_file);

%%
rev_point_set = nan(1,ses_num);
right_num_set = nan(1,ses_num);
rp_set = cell(1,ses_num);
pr_set = cell(1,ses_num);
neut_set = cell(1,ses_num);

%%
for i = 1:ses_num
    %%
    clearvars -except i ses_num beh_file rp_set pr_set neut_set rev_point_set right_num_set
    load(beh_file{i});
    
    %%
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;
    if outcomeIdentity(1,idenrev_cue(1))==3
        idenrev_cue = flipud(idenrev_cue);
    end
    
    lick_diff = diff(lickTime, 1);
    short_lick_diff = find(lick_diff(:,1)<=50)+1;
    lickTime_50ms = lickTime;
    lickTime_50ms(short_lick_diff,:)=[];
    
    lick_timings_50ms = cell(nTrial,1);
    delay_timings = stateTime(:,3:4)-stateTime(:,1)+1;
    delay_licknum_50ms = zeros(size(lickNum));
    for j = 1:nTrial
        jtrial_ind = find(lickTime_50ms(:,2)==j);
        lick_timings_50ms{j} = lickTime_50ms(jtrial_ind,1) - stateTime(j,1)+1;
        delay_licknum_50ms(j) = length(find((lick_timings_50ms{j}>=delay_timings(j,1)) & (lick_timings_50ms{j}<delay_timings(j,2))));
    end
    
    %%
    trial_win = 50;
    step_size = 25;
    
    start_trial = idenrev_trial-100;
    trial_range = start_trial:length(delay_licknum_50ms);
    
    rp_cue = idenrev_cue(1);
    pr_cue = idenrev_cue(2);
    cue_list = unique(odorCue+1);
    neut_cue = cue_list(~ismember(cue_list,idenrev_cue));
    
    temp_win = (1:trial_win)-step_size;
    istep = 0;
    temp_step_num = floor((length(trial_range)-trial_win)/step_size)+1;
    rp_lick = zeros(1,temp_step_num);
    pr_lick = zeros(1,temp_step_num);
    neut_lick = zeros(1,temp_step_num);
    temp_rev_point = nan;
    for istep = 1:temp_step_num
        temp_win = temp_win+step_size;
        if isnan(temp_rev_point) && (min(trial_range(temp_win))+trial_win/2)>=idenrev_trial
            temp_rev_point = istep;
        end
        temp_lick = delay_licknum_50ms(trial_range(temp_win));
        rp_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(rp_cue-1)));
        pr_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(pr_cue-1)));
        neut_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(neut_cue-1)));
    end
    
    rp_set{i} = rp_lick;
    pr_set{i} = pr_lick;
    neut_set{i} = neut_lick;
    rev_point_set(i) = temp_rev_point;
    right_num_set(i) = temp_step_num-temp_rev_point;
end

%%
x_left = min(rev_point_set);
x_right = min(right_num_set);
rp_set_cut = cell(1,ses_num);
pr_set_cut = cell(1,ses_num);
neut_set_cut = cell(1,ses_num);
for i = 1:ses_num
    if isnan(rev_point_set(i)); continue; end
    rp_set_cut{i} = rp_set{i}((rev_point_set(i)-x_left)+1:rev_point_set(i)+x_right);
    pr_set_cut{i} = pr_set{i}((rev_point_set(i)-x_left)+1:rev_point_set(i)+x_right);
    neut_set_cut{i} = neut_set{i}((rev_point_set(i)-x_left)+1:rev_point_set(i)+x_right);
end

%%
close all
x_plot = ((1-x_left:x_right))*step_size;

linewidth_set = [1.5 3];
fontsize_ax = 20;
fontsize_label = 24;
fig_size = [1 1 12 7.5];
r_color = [228 175 175; 139 23 28]/255;
b_color = [165 159 251; 57 89 168]/255;
n_color = [0.8 0.8 0.8; 0.5 0.5 0.5];
ylim_max = 10;
sigmaker_size = 10;

f1 = figure('PaperUnits','Centimeters','PaperPosition',fig_size);
ax1 = axes;
ax1.FontSize = fontsize_ax;
hold on
for i = 1:ses_num
    if isempty(rp_set_cut{i}); continue; end
    plot(x_plot, neut_set_cut{i}, 'color', n_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    plot(x_plot, rp_set_cut{i}, 'color', b_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    plot(x_plot, pr_set_cut{i}, 'color', r_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
end
ylim([0 ylim_max])
yticks(0:5:ylim_max)
xlim([min(x_plot)-step_size max(x_plot)+step_size])
plot([0 0],[0 100],'color',[0.5 0.5 0.5],'linestyle',':')
errorbar(x_plot, mean(cell2mat(neut_set_cut')), std(cell2mat(neut_set_cut'))/sqrt(ses_num-1),...
    'linewidth',linewidth_set(2),'color',n_color(2,:),'capsize',0)
errorbar(x_plot, mean(cell2mat(rp_set_cut')), std(cell2mat(rp_set_cut'))/sqrt(ses_num-1),...
    'linewidth',linewidth_set(2),'color',b_color(2,:),'capsize',0)
errorbar(x_plot, mean(cell2mat(pr_set_cut')), std(cell2mat(pr_set_cut'))/sqrt(ses_num-1),...
    'linewidth',linewidth_set(2),'color',r_color(2,:),'capsize',0)
ylabel('Anticipatory lick','fontsize',fontsize_label)
xlabel('Trial since reversal','fontsize',fontsize_label)

anova_result = zeros(size(x_plot));
data_set = {cell2mat(neut_set_cut') cell2mat(rp_set_cut') cell2mat(pr_set_cut')};
for istep = 1:length(x_plot)
    temp_data = [data_set{1}(:,istep) data_set{2}(:,istep) data_set{3}(:,istep)];
    [tbl,rm] = simple_mixed_anova(temp_data, [], {'Time'},{});
    if tbl.pValue(3) < 0.05
        c=multcompare(rm,'Time','ComparisonType','bonferroni');
        if x_plot(istep)<0
            anova_result(istep) = c.pValue(1)<0.05 & c.pValue(4)<0.05;
        else
            anova_result(istep) = c.pValue(2)<0.05 & c.pValue(4)<0.05;
        end
    end
end
anova_result(anova_result==0)=nan;
plot(x_plot,anova_result*ylim_max*0.95,'linestyle','none','marker','s','markeredgecolor','k','markerfacecolor','k',...
    'markersize',sigmaker_size)
