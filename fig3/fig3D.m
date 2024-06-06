clear; close all;

%% Mean normalized activity
%%
load('data_VIP_CC\no_rev_dataname.mat')

ses_num = length(cell_file);

trace_not_aligned_Rw_R_set = cell(ses_num,1);
trace_not_aligned_Rw_0_set = cell(ses_num,1);
trace_not_aligned_Pn_P_set = cell(ses_num,1);
trace_not_aligned_Pn_0_set = cell(ses_num,1);
trace_not_aligned_0_0_set = cell(ses_num,1);

trial_idx_rec_set = cell(1,ses_num);
neuron_num_set = zeros(1,ses_num);

%%
for i = 1:ses_num
    %% load
    clearvars -except ifile i cell_file beh_file ses_num aligned_file thr_pass neuron_drop rev_trial_anova_win ...
        trace_not_aligned_Rw_R_set trace_not_aligned_Rw_0_set trace_not_aligned_Pn_P_set trace_not_aligned_Pn_0_set trace_not_aligned_0_0_set ...
        trial_idx_rec_set neuron_num_set
    load(cell_file{i},'neuron');
    load(beh_file{i});
    load(aligned_file{i});
        
    %%
    [types, cues, types_rec, no_lick, no_lick1, outcomes] = recorded_trial_types(odorCue, waterReward, outcomeIdentity, stateTime, lickTime, trial_idx_rec, nTrial);
    for itype = 1:8
        eval(sprintf('type%d = types{%d};',itype,itype));
        eval(sprintf('type%d_rec = types_rec{%d};',itype,itype));
    end
    for icue = 1:4
        eval(sprintf('cue%d = cues{%d};',icue,icue));
    end
    
    no_lick = []; no_lick1 = [];
    
    rwded = outcomes{1};
    rwded_rec = outcomes{2};
    unrwded = outcomes{3};
    unrwded_rec = outcomes{4};
    
    puned = outcomes{5};
    puned_rec = outcomes{6};
    unpuned = outcomes{7};
    unpuned_rec = outcomes{8};
    
    epoch_cut_trial = 360;
    after_pass = (trial_idx_rec<thr_pass(i)+epoch_cut_trial)&(trial_idx_rec>=(thr_pass(i)));
    
    %% setting    
    C_raw = neuron.C_raw;

    S = full(neuron.S);
    neuron_num = size(C_raw,1);
    S_rate = zeros(1,neuron_num);
    for icell = 1:neuron_num
        S_rate(icell) = sum(S(icell,:)>0)/(size(S,2)/30);
    end
    event_thr_set = [0.015 0.025 0.020];
    hz_cut = 2;
    neuron_id = find(S_rate>event_thr_set(hz_cut));
    neuron_id(ismember(neuron_id, neuron_drop{i}))=[];
    neuron_num = length(neuron_id);
    neuron_num_set(i) = neuron_num;
    
    if neuron_num == 0; continue; end
    
    state_frame_num_rec = state_frame_num;
    stateTime_rec = stateTime;
    
    odorCue_rec = odorCue(trial_idx_rec);
    waterReward_rec = waterReward(trial_idx_rec);
    outcomeIdentity_rec = zeros(size(odorCue_rec));
    for itrial = 1:length(trial_idx_rec)
        outcomeIdentity_rec(itrial) = outcomeIdentity(trial_idx_rec(itrial),odorCue_rec(itrial)+1);
    end
    Prob = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    Ident = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    Reward = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    Cue = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    temp_itrial1 = 1;
    for itrial1 = trial_idx_rec
        if ismember(itrial1,no_lick)||ismember(itrial1,no_lick1); continue; end
        cue_l = odorCue(itrial1);
        rewarded_l = waterReward(itrial1);
        prob_l = outcomeProbability(itrial1,cue_l+1);
        id_l = outcomeIdentity(itrial1,cue_l+1);
        
        Cue(temp_itrial1) = cue_l+1;
        Ident(temp_itrial1) = id_l; % Outcome identity: 3 pn 2 rw
        Prob(temp_itrial1) = prob_l; % Reward probability
        Reward(temp_itrial1) = rewarded_l; % 0 1
        temp_itrial1 = temp_itrial1+1;
    end
    
    rw_period_set_frame = zeros(length(trial_idx_rec),1);
    rw_period_set_time = zeros(length(trial_idx_rec),1);
    period_diff_time = zeros(length(trial_idx_rec),1);
    j_rec = 0;
    for j = trial_idx_rec
        j_rec=j_rec+1;
        temp_period_diff_frame=0;
        temp_period_diff=0;
        if ismember(j,rwded)
            temp_delay_offset = stateTime(j,4);
            temp_lick_onset = lickTime(find(lickTime>=temp_delay_offset,1));
            temp_period_diff = temp_lick_onset - temp_delay_offset;
            temp_period_diff_frame = floor(temp_period_diff/1000*30);
        end
        period_diff_time(j_rec) = temp_period_diff;
        rw_period_set_frame(j_rec) = state_frame_num(j_rec,4)+temp_period_diff_frame;
        rw_period_set_time(j_rec) = stateTime(j_rec,4)+temp_period_diff;
    end
    
    lick_out_post = zeros(length(trial_idx_rec),1);
    speed_out_post = zeros(length(trial_idx_rec),1);
    for itrial = 1:length(trial_idx_rec)
        out_onset = rw_period_set_time(itrial);
        lick_out_post(itrial) = length(find((lickTime(:,1)>=out_onset) ...
            & (lickTime(:,1)<(out_onset + 1.5*1000))));
        speed_out_post(itrial) = length(find((cylinderTime(:,1)>=out_onset) ...
            & (cylinderTime(:,1)<(out_onset + 1.5*1000))));
    end
    
    %% indeces
    trial_idx_rec_set{i} = trial_idx_rec;
    
    baseline_onset = state_frame_num(:,1);
    baseline_offset = state_frame_num(:,2)-1;
    
    trace_not_aligned_Rw_R = nan(neuron_num,30*13+1);
    trace_not_aligned_Rw_0 = nan(neuron_num,30*13+1);
    trace_not_aligned_Pn_P = nan(neuron_num,30*13+1);
    trace_not_aligned_Pn_0 = nan(neuron_num,30*13+1);
    trace_not_aligned_0_0 = nan(neuron_num,30*13+1);
    
    temp_epoch = after_pass';
    
    for icell = 1:neuron_num
        temp_C_raw_smooth = movmean(C_raw(neuron_id(icell),:),15); % 0.5 s smooth
        temp_C_raw = C_raw(neuron_id(icell),:);

        temp_baseline_smooth = cell(1,length(baseline_onset));
        temp_baseline = cell(1,length(baseline_onset));
        for j3 = 1:length(baseline_onset)
            temp_baseline_smooth{j3} = mean(temp_C_raw_smooth(baseline_onset(j3):baseline_offset(j3)));
            temp_baseline{j3} = mean(temp_C_raw(baseline_onset(j3):baseline_offset(j3)));
        end
        temp_baseline_mat_smooth = cell2mat(temp_baseline_smooth);
        temp_baseline_mat = cell2mat(temp_baseline);
        temp_C_raw_z_smooth = (temp_C_raw_smooth-mean(temp_baseline_mat_smooth))/std(temp_baseline_mat_smooth); % baseline z
        temp_C_raw_z = (temp_C_raw-mean(temp_baseline_mat))/std(temp_baseline_mat); % baseline z

        trial_trace_cue = nan(length(trial_idx_rec),5*30+1);
        trial_trace_out = nan(length(trial_idx_rec),8*30+1);
        for itrial = 1:length(trial_idx_rec)
            trial_trace_cue(itrial,:) = (temp_C_raw_z_smooth(state_frame_num_rec(itrial,2)-0.5*30:state_frame_num_rec(itrial,2)+4.5*30));
            trial_trace_out(itrial,:) = (temp_C_raw_z_smooth(rw_period_set_frame(itrial)-0.5*30:rw_period_set_frame(itrial)+7.5*30));
        end

        % not ialigned
        trace_not_aligned_Rw_R(icell,5*30+1:end) = mean(trial_trace_out(Ident==2&Prob==75&Reward==1&temp_epoch,:));
        trace_not_aligned_Rw_0(icell,5*30+1:end) = mean(trial_trace_out(Ident==2&Prob==75&Reward==0&temp_epoch,:));
        trace_not_aligned_Pn_P(icell,5*30+1:end) = mean(trial_trace_out(Ident==3&Prob==75&Reward==1&temp_epoch,:));
        trace_not_aligned_Pn_0(icell,5*30+1:end) = mean(trial_trace_out(Ident==3&Prob==75&Reward==0&temp_epoch,:));
        trace_not_aligned_0_0(icell,5*30+1:end) = mean(trial_trace_out(Prob==0&temp_epoch,:));

        trace_not_aligned_Rw_R(icell,1:5*30) = mean(trial_trace_cue(Ident==2&Prob==75&temp_epoch,1:5*30));
        trace_not_aligned_Pn_P(icell,1:5*30) = mean(trial_trace_cue(Ident==3&Prob==75&temp_epoch,1:5*30));
        trace_not_aligned_0_0(icell,1:5*30) = mean(trial_trace_cue(Prob==0&temp_epoch,1:5*30));
    end
    
    trace_not_aligned_Rw_R_set{i} = trace_not_aligned_Rw_R;
    trace_not_aligned_Rw_0_set{i} = trace_not_aligned_Rw_0;
    trace_not_aligned_Pn_P_set{i} = trace_not_aligned_Pn_P;
    trace_not_aligned_Pn_0_set{i} = trace_not_aligned_Pn_0;
    trace_not_aligned_0_0_set{i} = trace_not_aligned_0_0;
    
end

%% fig
temp_trace_set = {trace_not_aligned_Rw_R_set, trace_not_aligned_Rw_0_set, trace_not_aligned_Pn_P_set, trace_not_aligned_Pn_0_set, trace_not_aligned_0_0_set};

plot_x1 = -0.5:1/30:4.5;
plot_x2 = -0.5:1/30:7.5;
cmap1 = [30 45 232; 30 45 232; 255 30 70; 255 30 70; 128 128 128]/255;
line_type = {'-','--','-','--','-'};

f_out = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 5]);
s1 = subplot(1,11,1:5);
hold on
fill_x=[plot_x1 fliplr(plot_x1)];
for ifill = 1:5
    temp_data = temp_trace_set{ifill};
    temp_fill = cell2mat(temp_data(:,1));
    temp_fill = temp_fill(:,1:length(plot_x1));
    fill_y=[mean(temp_fill)+std(temp_fill,0,1)/sqrt(size(temp_fill,1)) ...
        fliplr(mean(temp_fill)-std(temp_fill,0,1)/sqrt(size(temp_fill,1)))];
    fill_hist = fill(fill_x,fill_y,cmap1(ifill,:));
    fill_hist.EdgeColor = 'none';
    alpha(fill_hist,0.2)
end
for ifill = 1:5
    temp_data = temp_trace_set{ifill};
    temp_fill = cell2mat(temp_data(:,1));
    temp_fill = temp_fill(:,1:length(plot_x1));
    plot(plot_x1,mean(temp_fill),'color',cmap1(ifill,:),'linestyle',line_type{ifill},'linewidth',1.5)
end
axisVal = axis;
ylim([axisVal(3) axisVal(4)])
plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5]);
plot([1.5 1.5],[-10 10],'--','color',[0.5 0.5 0.5]);
plot([3.5 3.5],[-10 10],'--','color',[0.5 0.5 0.5]);
plot([-0.5 12.5],[0 0],'--','color',[0.5 0.5 0.5]);
xlim([-0.5 4.5-1/30])
xticks(0:4)

s2 = subplot(1,11,6:11);
hold on
fill_x=[plot_x2 fliplr(plot_x2)];
for ifill = 1:5
    temp_data = temp_trace_set{ifill};
    temp_fill = cell2mat(temp_data(:,1));
    temp_fill = temp_fill(:,end-length(plot_x2)+1:end);
    fill_y=[mean(temp_fill)+std(temp_fill,0,1)/sqrt(size(temp_fill,1)) ...
        fliplr(mean(temp_fill)-std(temp_fill,0,1)/sqrt(size(temp_fill,1)))];
    fill_hist = fill(fill_x,fill_y,cmap1(ifill,:));
    fill_hist.EdgeColor = 'none';
    alpha(fill_hist,0.2)
end
for ifill = 1:5
    temp_data = temp_trace_set{ifill};
    temp_fill = cell2mat(temp_data(:,1));
    temp_fill = temp_fill(:,end-length(plot_x2)+1:end);
    plot(plot_x2,mean(temp_fill),'color',cmap1(ifill,:),'linestyle',line_type{ifill},'linewidth',1.5)
end
axisVal = axis;
ylim([axisVal(3) axisVal(4)])
xlim([-0.5 5.5])
plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5]);
plot([2.5 2.5],[-10 10],'--','color',[0.5 0.5 0.5]);
plot([-2.5 12.5],[0 0],'--','color',[0.5 0.5 0.5]);
yticklabels({''})
xticks(0:5)

linkaxes([s1,s2],'y')
ylim([-0.1 0.2])
