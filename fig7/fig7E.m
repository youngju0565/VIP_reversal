clear; close all;

%% Optogenetics: anticipatory lick rate
main_path = 'E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\ChRmine';
cd(main_path)

%%
file_list = ls('*.mat');
file_num = size(file_list,1);

if mod(file_num,2)>0; disp('File number error: not even'); return; end

mod_order = repmat([1 2],1,file_num/2)==1;

%%
% parameters and pre-allocations
trial_win = 60;
step_size = 15;

rp_control_set = cell(1,file_num);
rp_laser_set = cell(1,file_num);
pr_control_set = cell(1,file_num);
pr_laser_set = cell(1,file_num);
rev_point_set = zeros(2,file_num); % control - laser
right_num_set = zeros(2,file_num); % control - laser

for ises = 1:file_num
    load(file_list(ises,:))
    first_mod = mod_order(ises);
    ses_name = strsplit(file_list(ises,:),'_');

    %% reversal 
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;

    %% lick time
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

    %%
    trial_cut = [idenrev_trial; nTrial];
    rp_set = cell(1,2);
    pr_set = cell(1,2);
    rev_point = nan(1,2);
    right_num = nan(1,2);

    for irev = 1:2
        %%
        if irev==1
            start_trial = idenrev_trial(1)-100;
        elseif irev==2
            start_trial = idenrev_trial(2)-100;
        end
        trial_range = start_trial:trial_cut(irev+1)-1;

        outcomeID =  outcomeIdentity(idenrev_trial(irev),:);
        pn_cue = find(outcomeID==3);
        rw_cue = find(outcomeID==2);

        [~,rev_cue] = find(diff(outcomeIdentity(idenrev_trial(irev)-1:idenrev_trial(irev),:))); % RP
        rp_cue = pn_cue(ismember(pn_cue,rev_cue));
        pr_cue = rw_cue(ismember(rw_cue,rev_cue));

        temp_win = (1:trial_win)-step_size;
        istep = 0;
        temp_step_num = floor((length(trial_range)-trial_win)/step_size)+1;
        rp_lick = zeros(1,temp_step_num);
        pr_lick = zeros(1,temp_step_num);
        temp_rev_point = nan;
        for istep = 1:temp_step_num
            temp_win = temp_win+step_size;
            if isnan(temp_rev_point) && (min(trial_range(temp_win))+trial_win/2)>=idenrev_trial(irev)
                temp_rev_point = istep;
            end
            temp_lick = delay_licknum_50ms(trial_range(temp_win));
            rp_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(rp_cue-1)));
            pr_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(pr_cue-1)));
        end

        rp_set{irev} = rp_lick;
        pr_set{irev} = pr_lick;
        rev_point(irev) = temp_rev_point;
        right_num(irev) = temp_step_num-temp_rev_point;
    end

    %%
    if first_mod
        rp_control_set{ises} = rp_set{2};
        rp_laser_set{ises} = rp_set{1};
        pr_control_set{ises} = pr_set{2};
        pr_laser_set{ises} = pr_set{1};
        rev_point_set(:,ises) = flipud(rev_point');
        right_num_set(:,ises) = flipud(right_num');
    else
        rp_control_set{ises} = rp_set{1};
        rp_laser_set{ises} = rp_set{2};
        pr_control_set{ises} = pr_set{1};
        pr_laser_set{ises} = pr_set{2};
        rev_point_set(:,ises) = rev_point';
        right_num_set(:,ises) = right_num';
    end
end

%%
% cut
x_left = min(min(rev_point_set));
x_right = min(min(right_num_set));
rp_control_cut = cell(1,file_num);
rp_laser_cut = cell(1,file_num);
pr_control_cut = cell(1,file_num);
pr_laser_cut = cell(1,file_num);
for ises = 1:file_num
    rp_control_cut{ises} = rp_control_set{ises}((rev_point_set(1,ises)-x_left)+1:rev_point_set(1,ises)+x_right);
    pr_control_cut{ises} = pr_control_set{ises}((rev_point_set(1,ises)-x_left)+1:rev_point_set(1,ises)+x_right);
    rp_laser_cut{ises} = rp_laser_set{ises}((rev_point_set(2,ises)-x_left)+1:rev_point_set(2,ises)+x_right);
    pr_laser_cut{ises} = pr_laser_set{ises}((rev_point_set(2,ises)-x_left)+1:rev_point_set(2,ises)+x_right);
end

% animal mean
animal_num = file_num/2;
rp_control_set_cut = cell(1,animal_num);
rp_laser_set_cut = cell(1,animal_num);
pr_control_set_cut = cell(1,animal_num);
pr_laser_set_cut = cell(1,animal_num);
for ianimal = 1:animal_num
    rp_control_set_cut{ianimal} = (rp_control_cut{ianimal*2-1}+rp_control_cut{ianimal*2})/2;
    pr_control_set_cut{ianimal} = (pr_control_cut{ianimal*2-1}+pr_control_cut{ianimal*2})/2;
    rp_laser_set_cut{ianimal} = (rp_laser_cut{ianimal*2-1}+rp_laser_cut{ianimal*2})/2;
    pr_laser_set_cut{ianimal} = (pr_laser_cut{ianimal*2-1}+pr_laser_cut{ianimal*2})/2;
end

%% plot
x_plot = ((1-x_left:x_right))*step_size;

linewidth_set = [1.5 3];
fontsize_ax = 20;
fontsize_label = 24;
r_color = [228 175 175; 139 23 28]/255;
b_color = [165 159 251; 57 89 168]/255;
ylim_max = 7;
sigmaker_size = 10;

f_nolaser = figure('PaperUnits','Centimeters','PaperPosition',[1 1 9.5 9]);
ax_nolaser=axes;
ax_nolaser.FontSize = fontsize_ax;
hold on
title('Laser OFF','fontsize',fontsize_label)
data1 = rp_control_set_cut;
data2 = pr_control_set_cut;
for ises = 1:animal_num
    plot(x_plot, data1{ises}, 'color', b_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    plot(x_plot, data2{ises}, 'color', r_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
end
ylim([0 ylim_max])
yticks([0 3 6])
xlim([min(x_plot)-step_size max(x_plot)+step_size])
plot([0 0],[0 100],'color',[0.5 0.5 0.5],'linestyle',':')
errorbar(x_plot, mean(cell2mat(data1')), std(cell2mat(data1'))/sqrt(animal_num),...
    'linewidth',linewidth_set(2),'color',b_color(2,:),'capsize',0)
errorbar(x_plot, mean(cell2mat(data2')), std(cell2mat(data2'))/sqrt(animal_num),...
    'linewidth',linewidth_set(2),'color',r_color(2,:),'capsize',0)
xlabel('Trial since reversal','fontsize',fontsize_label)
ylabel('Anticipatory lick','fontsize',fontsize_label)
ttest_result = zeros(size(x_plot));
temp1 = cell2mat(data1');
temp2 = cell2mat(data2');
for istep = 1:length(x_plot)
    ttest_result(istep) = ttest(temp1(:,istep),temp2(:,istep));
end
ttest_result(ttest_result==0)=nan;
plot(x_plot,ttest_result*ylim_max*0.95,'linestyle','none','marker','s','markeredgecolor','k','markerfacecolor','k',...
    'markersize',sigmaker_size)
plot([x_plot(find(ttest_result==1&mean(cell2mat(data2'))>mean(cell2mat(data1')),1)) x_plot(find(ttest_result==1&mean(cell2mat(data2'))>mean(cell2mat(data1')),1))],...
    [0 ylim_max*0.95],'--k','linewidth',2)

f_laser = figure('PaperUnits','Centimeters','PaperPosition',[1 1 9.5 9]);
ax_laser=axes;
ax_laser.FontSize = fontsize_ax;
hold on
title('Laser ON','fontsize',fontsize_label)
data1 = rp_laser_set_cut;
data2 = pr_laser_set_cut;
for ises = 1:animal_num
    plot(x_plot, data1{ises}, 'color', b_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    plot(x_plot, data2{ises}, 'color', r_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
end
ylim([0 ylim_max])
yticks([0 3 6])
xlim([min(x_plot)-step_size max(x_plot)+step_size])
plot([0 0],[0 100],'color',[0.5 0.5 0.5],'linestyle',':')
errorbar(x_plot, mean(cell2mat(data1')), std(cell2mat(data1'))/sqrt(animal_num),...
    'linewidth',linewidth_set(2),'color',b_color(2,:),'capsize',0)
errorbar(x_plot, mean(cell2mat(data2')), std(cell2mat(data2'))/sqrt(animal_num),...
    'linewidth',linewidth_set(2),'color',r_color(2,:),'capsize',0)
xlabel('Trial since reversal','fontsize',fontsize_label)
ylabel('Anticipatory lick','fontsize',fontsize_label)
ttest_result = zeros(size(x_plot));
temp1 = cell2mat(data1');
temp2 = cell2mat(data2');
for istep = 1:length(x_plot)
    ttest_result(istep) = ttest(temp1(:,istep),temp2(:,istep));
end
ttest_result(ttest_result==0)=nan;
plot(x_plot,ttest_result*ylim_max*0.95,'linestyle','none','marker','s','markeredgecolor','k','markerfacecolor','k',...
    'markersize',sigmaker_size)
plot([x_plot(find(ttest_result==1&mean(cell2mat(data2'))>mean(cell2mat(data1')),1)) x_plot(find(ttest_result==1&mean(cell2mat(data2'))>mean(cell2mat(data1')),1))],...
    [0 ylim_max*0.95],'--k','linewidth',2)
