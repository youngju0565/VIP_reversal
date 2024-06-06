clear; close all;

%% Cue & Outcome dependent activity
%%
load('data_VIP_CC\no_rev_dataname.mat')

ses_num = length(cell_file);

trial_idx_rec_set = cell(1,ses_num);
neuron_num_set = zeros(1,ses_num);

d1_mean_set = cell(ses_num,1);
out_mean_set = cell(ses_num,1);

%%
for i = 1:ses_num
    %% load
    clearvars -except ifile i cell_file beh_file ses_num aligned_file thr_pass neuron_drop rev_trial_anova_win ...
        trial_idx_rec_set neuron_num_set d1_mean_set out_mean_set
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
    
    temp_epoch = after_pass';
    
    d1_mean = zeros(neuron_num,3);
    out_mean = zeros(neuron_num,4);
    
    for icell = 1:neuron_num
        temp_C_raw = C_raw(neuron_id(icell),:);
        
        temp_baseline = cell(1,length(baseline_onset));
        for j3 = 1:length(baseline_onset)
            temp_baseline{j3} = mean(temp_C_raw(baseline_onset(j3):baseline_offset(j3)));
        end
        temp_baseline_mat = cell2mat(temp_baseline);
        temp_C_raw_z = (temp_C_raw-mean(temp_baseline_mat))/std(temp_baseline_mat); % baseline z
        
        trial_cut_out = zeros(length(trial_idx_rec),1);
        trial_cut_d1 = zeros(length(trial_idx_rec),1);
        for itrial = 1:length(trial_idx_rec)
            trial_cut_out(itrial) = mean(temp_C_raw_z(rw_period_set_frame(itrial):rw_period_set_frame(itrial)+1.5*30));
            trial_cut_d1(itrial) = mean(temp_C_raw_z(state_frame_num_rec(itrial,3):state_frame_num_rec(itrial,3)+1.5*30));
        end
        
        % mean activity per cues/types
        d1_mean(icell,1) = mean(trial_cut_d1(Ident==2&Prob==75&temp_epoch)); % Rw
        d1_mean(icell,2) = mean(trial_cut_d1(Ident==3&Prob==75&temp_epoch)); % Pn
        d1_mean(icell,3) = mean(trial_cut_d1(Prob==0&temp_epoch)); % Nt
        
        out_mean(icell,1) = mean(trial_cut_out(Ident==2&Prob==75&Reward==1&temp_epoch)); % Rw-Rwded
        out_mean(icell,2) = mean(trial_cut_out(Ident==2&Prob==75&Reward==0&temp_epoch)); % Rw-omit
        out_mean(icell,3) = mean(trial_cut_out(Ident==3&Prob==75&Reward==1&temp_epoch)); % Pn-puned
        out_mean(icell,4) = mean(trial_cut_out(Ident==3&Prob==75&Reward==0&temp_epoch)); % Pn-omit
    end
    
    d1_mean_set{i} = d1_mean;
    out_mean_set{i} = out_mean;
end

%%
temp_cue_data = cell2mat(d1_mean_set); % R-P-0
temp_cue_data = temp_cue_data(:,1)-temp_cue_data(:,2);

temp_out_data = cell2mat(out_mean_set); % RR R0 PP P0
temp_out_data = [temp_out_data(:,1)-temp_out_data(:,2) temp_out_data(:,3)-temp_out_data(:,4)];

temp_data = [temp_cue_data temp_out_data];

cmap1 = [0.3 0.3 0.3; 30/255 45/255 232/255; 255/255 30/255 70/255];

%%
f_bar3 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.5 3.5]);
clf
ax = axes;
set(gca,'xcolor','k','ycolor','k')
hold on
for ibar = 1:3
    s1 = scatter(rand(1,size(temp_data,1))*0.5-0.25+ibar,temp_data(:,ibar));
    s1.MarkerFaceColor = cmap1(ibar,:);
    s1.MarkerEdgeColor = 'none';
    s1.SizeData = 3;
    alpha(s1,0.3)
    b1 = bar(ibar,mean(temp_data(:,ibar)),0.5);
    b1.EdgeColor = cmap1(ibar,:); b1.LineWidth = 2; b1.FaceColor = cmap1(ibar,:);
    e1 = errorbar(ibar,mean(temp_data(:,ibar)),std(temp_data(:,ibar))/sqrt(size(temp_data,1)));
    e1.Color = 'k'; e1.CapSize = 0; e1.LineWidth = 2;
end
xticks(1:3)
xlim([0.5 3.5])
xticklabels({'Cue','Rw','Pn'})

ylim([-1 1])
axisVal = axis;
for ibar = 1:3
    [~,p] = ttest(temp_data(:,ibar));
    if p<0.05
        star_str = '*';
        if p<0.01; star_str = '**'; end
        if p<0.001; star_str = '***'; end
        text(ibar,axisVal(3)+225/250*(axisVal(4)-axisVal(3)),star_str,...
            'fontsize',15,'horizontalalignment','center')
    end
end
