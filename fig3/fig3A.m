clear; close all;

%% Example session: lick trace
%%
load('data_VIP_CC\no_rev_dataname.mat')

ses_num = length(beh_file);

%%
i = 3;
load(beh_file{i});

[stateTime, stateTime_original, state_duration, trial_duration, state_time] = stateTime_zerofil(stateTime);

%%
if size(lickTime,1)>1
    lick_diff = diff(lickTime, 1);
    short_lick_diff = find(lick_diff(:,1)<=80)+1;
    lickTime_80ms = lickTime;
    lickTime_80ms(short_lick_diff,:)=[];
elseif size(lickTime,1)==1
    lickTime_80ms = lickTime;
elseif isempty(lickTime)
    fprintf('No lick in %s\n',filename)
    return
end

lick_timings_80ms = cell(nTrial,1);
delay_timings = stateTime(:,3:4)-stateTime(:,1)+1;
delay_licknum_80ms = zeros(size(lickNum));
for j = 1:nTrial
    jtrial_ind = find(lickTime_80ms(:,2)==j);
    lick_timings_80ms{j} = lickTime_80ms(jtrial_ind,1) - stateTime(j,1)+1;
    delay_licknum_80ms(j) = length(find((lick_timings_80ms{j}>=delay_timings(j,1)) & (lick_timings_80ms{j}<delay_timings(j,2))));
end

cylinder_timings = cell(nTrial,1);
delay_cylinder = zeros(size(lickNum));
for itrial = 1:nTrial
    itrial_ind = find(cylinderTime(:,2)==itrial);
    cylinder_timings{itrial} = cylinderTime(itrial_ind,1) - stateTime(itrial,1)+1;
    delay_cylinder(itrial) = length(find((cylinder_timings{itrial}>=delay_timings(itrial,1)) & (cylinder_timings{itrial}<(delay_timings(itrial,1)+1.5*1000))));
end

%%
epoch_cut_trial = 360;

type1 = []; type2 = []; type3 = []; type4 = []; type5 = []; type6 = []; type7 = []; type8 = [];
for k = thr_pass(i):thr_pass(i)+epoch_cut_trial-1
    cue = odorCue(k);
    rewarded = waterReward(k);
    switch cue
        case 0
            if rewarded; type1 = [type1; k];
            else type2 = [type2; k]; end
        case 1
            if rewarded; type3 = [type3; k];
            else type4 = [type4; k]; end
        case 2
            if rewarded; type5 = [type5; k];
            else type6 = [type6; k]; end
        case 3
            if rewarded; type7 = [type7; k];
            else type8 = [type8; k]; end
    end
end
types = {type1, type2, type3, type4, type5, type6, type7, type8};
cues = {sort([type1;type2]), sort([type3;type4]), sort([type5;type6]), sort([type7;type8])};

%%
cmap1 = [30 45 232; 128 128 128; 232 126 58; 255 30 70]/255;

type_mean = zeros(8, trial_duration*10^3/10^2);
type_set_cell = cell(1,8);
for o = 1:8
    if isempty(types{o}); continue; end
    type_trials = types{o};
    type_count = zeros(1, trial_duration*10^3/10^2);
    type_set = zeros(length(type_trials), trial_duration*10^3/10^2);
    for p = 1:length(type_trials)
        p_i = type_trials(p);
        trial_lick_times = lick_timings_80ms{p_i};
        [counts, edges] = histcounts(trial_lick_times, 0:10^2:trial_duration*10^3);
        type_set(p,:) = smoothdata(counts,'movmean',5)*10;
        type_count = type_count + counts;
    end
    type_mean(o,:) = type_count * 10 / size(type_trials,1);
    type_set_cell{o} = type_set;
end

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3]);
ax = axes;
hold on
ax_fontsize = 10;
ax.FontSize = ax_fontsize;
label_fontsize = 10;
set(gca,'xcolor','k','ycolor','k')
legend_vec2={}; legend_idx2 = 0;
fill_x = [edges(2:end) fliplr(edges(2:end))];
for q = 1:8
    if isempty(types{q}); continue; end
    cue_q = ceil(q/2);
    if ~isempty(find([2 4 8]==q,1))
        fill_y = [smoothdata(type_mean(q,:), 'movmean',5)-std(type_set_cell{q})/sqrt(length(types{q})),...
            fliplr(smoothdata(type_mean(q,:), 'movmean',5)+std(type_set_cell{q})/sqrt(length(types{q})))];
        fill_ = fill(fill_x, fill_y, cmap1(cue_q,:));
        fill_.EdgeColor = 'none';
        alpha(fill_, 0.2)
        plot(edges(2:end), smoothdata(type_mean(q,:), 'movmean',5), '--','color', cmap1(cue_q,:), 'linewidth',1)
    else
        legend_idx2 = legend_idx2+1;
        fill_y = [smoothdata(type_mean(q,:), 'movmean',5)-std(type_set_cell{q})/sqrt(length(types{q})),...
            fliplr(smoothdata(type_mean(q,:), 'movmean',5)+std(type_set_cell{q})/sqrt(length(types{q})))];
        fill_ = fill(fill_x, fill_y, cmap1(cue_q,:));
        fill_.EdgeColor = 'none';
        alpha(fill_, 0.2)
        plot(edges(2:end), smoothdata(type_mean(q,:), 'movmean',5), 'color', cmap1(cue_q,:), 'linewidth',1);
    end
end
axisValue = axis;
line([state_time(3)*10^3 state_time(3)*10^3],[0 axisValue(4)],'color', [0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
line([state_time(2)*10^3 state_time(2)*10^3],[0 axisValue(4)],'color', [0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
line([state_time(1)*10^3 state_time(1)*10^3],[0 axisValue(4)],'color', [0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
line([4*10^3 4*10^3],[0 axisValue(4)],'color', [0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
line([8*10^3 8*10^3],[0 axisValue(4)],'color', [0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
ylim([0 axisValue(4)])
xticks ([0.5*10^3 : 10^3 : trial_duration*10^3])
xticklabels([0:1:trial_duration])
xlim([3.5*10^3 9*10^3])
xlim([0*10^3 9*10^3])
ylabel('Lick rate (Hz)','FontSize',label_fontsize)
xlabel('Time from cue onset (s)','FontSize',label_fontsize)
