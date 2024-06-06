clear; close all;

%% RP reversal: Cue-dependent (Delay1 period)
%%
load('data_VIP_CC\rev_dataname.mat')

ses_num = length(cell_file);

%%
cue_dep_d1_set = cell(1,ses_num);
cue_dep_trial_set = cell(1,ses_num);

idenrev_trial_set = zeros(1,ses_num);
trial_idx_rec_set = cell(1,ses_num);

%%
for i = 1:ses_num
    %% load
    clearvars -except i cell_file beh_file ses_num aligned_file rev_trial_anova_win active_cell_S neurons neuron_drop thr_pass ...
        cue_dep_d1_set idenrev_trial_set trial_idx_rec_set cue_dep_trial_set
    load(cell_file{i},'neuron');
    load(beh_file{i});
    load(aligned_file{i});
    
    %% recorded trial types
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;
    if outcomeIdentity(1,idenrev_cue(1))==3
        idenrev_cue = flipud(idenrev_cue);
    end
    idenrev_trial_set(i) = idenrev_trial;
    
    [types, cues, types_rec, no_lick, no_lick1, outcomes] = recorded_trial_types(odorCue, waterReward, outcomeIdentity, stateTime, lickTime, trial_idx_rec, nTrial);
    for itype = 1:8
        eval(sprintf('type%d = types{%d};',itype,itype));
        eval(sprintf('type%d_rec = types_rec{%d};',itype,itype));
    end
    for icue = 1:4
        eval(sprintf('cue%d = cues{%d};',icue,icue));
    end
    
    if ~isempty(no_lick)||~isempty(no_lick1)
        [~,no_lick_idx] = ismember(no_lick, trial_idx_rec);
        [~,no_lick1_idx] = ismember(no_lick1, trial_idx_rec);
        
        out_trial_rec = false(size(trial_idx_rec));
        out_trial_rec([no_lick_idx; no_lick1_idx])=true;
    else
        out_trial_rec = false(size(trial_idx_rec));
    end

    rwded = outcomes{1};
    rwded_rec = outcomes{2};
    unrwded = outcomes{3};
    unrwded_rec = outcomes{4};
    
    puned = outcomes{5};
    puned_rec = outcomes{6};
    unpuned = outcomes{7};
    unpuned_rec = outcomes{8};
    
    epoch_cut_trial = 100;
    before_rev = ~out_trial_rec&(trial_idx_rec<idenrev_trial)&(trial_idx_rec>=(idenrev_trial-epoch_cut_trial));
    during_rev = ~out_trial_rec&(trial_idx_rec<idenrev_trial+epoch_cut_trial)&(trial_idx_rec>=(idenrev_trial));
    after_rev = ~out_trial_rec&(trial_idx_rec>=idenrev_trial+150)&(trial_idx_rec<(idenrev_trial+150+epoch_cut_trial));
    
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

    if neuron_num == 0; continue; end
    
    state_frame_num_rec = state_frame_num;
    stateTime_rec = stateTime;
    
    odorCue_rec = odorCue(trial_idx_rec);
    waterReward_rec = waterReward(trial_idx_rec);
    outcomeIdentity_rec = zeros(size(odorCue_rec));
    for itrial = 1:length(trial_idx_rec)
        outcomeIdentity_rec(itrial) = outcomeIdentity(trial_idx_rec(itrial),odorCue_rec(itrial)+1);
    end
    Prob = zeros(length(trial_idx_rec), 1);
    Ident = zeros(length(trial_idx_rec), 1);
    Reward = zeros(length(trial_idx_rec), 1);
    Cue = zeros(length(trial_idx_rec), 1);
    temp_itrial1 = 1;
    for itrial1 = trial_idx_rec
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
    
    %% indeces
    trial_idx_rec_set{i} = trial_idx_rec;
    
    baseline_onset = state_frame_num(:,1);
    baseline_offset = state_frame_num(:,2)-1;

    cue_dep_d1 = zeros(neuron_num,3);
    for icell = 1:neuron_num
        temp_C_raw = C_raw(neuron_id(icell),:);
        temp_baseline = cell(1,length(baseline_onset));
        for j3 = 1:length(baseline_onset)
            temp_baseline{j3} = mean(temp_C_raw(baseline_onset(j3):baseline_offset(j3)));
        end
        temp_baseline_mat = cell2mat(temp_baseline);
        temp_C_raw_z = (temp_C_raw-mean(temp_baseline_mat))/std(temp_baseline_mat);
        
        temp_delay1 = zeros(length(trial_idx_rec),3);
        for itrial = 1:length(trial_idx_rec) % for each recorded trial
            rec_itrial = trial_idx_rec(itrial);
            idelay_onset = state_frame_num(itrial,3);
            temp_delay1(itrial,1) = mean(temp_C_raw_z(idelay_onset:idelay_onset+30*1.5));
            temp_delay1(itrial,2) = mean(temp_C_raw_z(idelay_onset:idelay_onset+30*1.5));
            temp_delay1(itrial,3) = mean(temp_C_raw_z(idelay_onset:idelay_onset+30*1.5));
        end

        cue_dep_d1(icell,1) = mean(temp_delay1(before_rev'&odorCue_rec==(idenrev_cue(1)-1),1))-mean(temp_delay1(before_rev'&odorCue_rec==(idenrev_cue(2)-1),1));
        cue_dep_d1(icell,2) = mean(temp_delay1(during_rev'&odorCue_rec==(idenrev_cue(1)-1),2))-mean(temp_delay1(during_rev'&odorCue_rec==(idenrev_cue(2)-1),2));
        cue_dep_d1(icell,3) = mean(temp_delay1(after_rev'&odorCue_rec==(idenrev_cue(1)-1),3))-mean(temp_delay1(after_rev'&odorCue_rec==(idenrev_cue(2)-1),3));
    end
    cue_dep_d1_set{i} = cue_dep_d1;
    
    delay_mean = zeros(neuron_num,length(trial_idx_rec));
    j2 = 0;
    for j1 = neuron_id
        temp_C_raw = C_raw(j1,:);
        temp_baseline = cell(1,length(baseline_onset));
        for j3 = 1:length(baseline_onset)
            temp_baseline{j3} = mean(temp_C_raw(baseline_onset(j3):baseline_offset(j3)));
        end
        temp_baseline_mat = cell2mat(temp_baseline);
        temp_C_raw_z = (temp_C_raw-mean(temp_baseline_mat))/std(temp_baseline_mat); % baseline z
            
        temp_neuron = zeros(1,length(trial_idx_rec));
        for k = 1:length(trial_idx_rec)
            temp_delay = temp_C_raw_z(state_frame_num(k,3):state_frame_num(k,3)+30*1.5); % no lick delay
            temp_neuron(k) = mean(temp_delay);
        end
        j2 = j2+1;
        delay_mean(j2,:)=temp_neuron;
    end
    
    win_size_trial = 41;
    index_mat = zeros(neuron_num, length(trial_idx_rec)-win_size_trial+1);
    for l=1:(length(trial_idx_rec)-win_size_trial+1)
        temp_bin = delay_mean(:,l:l-1+win_size_trial);
        % index = (cue2-cue1)
        temp_trials = trial_idx_rec(l:l-1+win_size_trial);
        temp_cue1_trials = find(odorCue(temp_trials)==(idenrev_cue(1)-1));
        temp_cue2_trials = find(odorCue(temp_trials)==(idenrev_cue(2)-1));
        temp_cue1 = mean(temp_bin(:,temp_cue1_trials),2);
        temp_cue2 = mean(temp_bin(:,temp_cue2_trials),2);
        temp_index = (temp_cue1-temp_cue2);
        index_mat(:,l) = temp_index;
    end
    cue_dep_trial_set{i} = index_mat;
end

%%
colorscheme = othercolor('BuDRd_18');
close all;

%% cue dependent activity scatter
fontsize_ax = 20; 
fontsize_label = 24;

cue_dep_d1_set_mat = cell2mat(cue_dep_d1_set');

f_scatter = figure('PaperUnits','Centimeters','PaperPosition',[2 2 11 11]);
ax_scatter = axes;
hold on
ax_scatter.FontSize = fontsize_ax;

data_draw = cue_dep_d1_set_mat; x = 1; y = 3;
marker_color = 0.4;
scatter(data_draw(:,x), data_draw(:,y),'filled','markeredgecolor',[1 1 1]*marker_color,'markerfacecolor',[1 1 1]*marker_color)
xlim([-1 1])
ylim([-1 1])
xticks(-1:1:1)
yticks(-1:1:1)
xline1 = line([-2000 2000],[0 0]);
xline1.Color = 'k'; xline1.LineStyle = '--';
yline1 = line([0 0],[-2000 2000]);
yline1.Color = 'k'; yline1.LineStyle = '--';
fit1 = polyfit(data_draw(:,x),data_draw(:,y),1);
plot([min(data_draw(:,x)) max(data_draw(:,x))]*1.2,polyval(fit1,[min(data_draw(:,x)) max(data_draw(:,x))]*1.2),'Color','k','linewidth',1.5)
[r,p] = corrcoef(data_draw(:,x),data_draw(:,y));
str1 = {['{\itr} = ',num2str(r(1,2),'%.3f')],['{\itp} = ',num2str(p(1,2),'%.3f')]};
annotation('textbox',[0.7 0.85 0.15 0.1],'String',str1,'FitBoxToText','on','edgecolor','none','fontsize',14)
annotation('textbox',[0.2 0.85 0.15 0.1],'String','Delay','FitBoxToText','on','edgecolor','none','fontsize',24)
xlabel('Before','fontsize',fontsize_label)
ylabel('After','fontsize',fontsize_label)

%% trial window cue dependent activity
before_rev_set = zeros(1,ses_num);
after_rev_set = zeros(1,ses_num);
for ises = 1:ses_num
    temp_data = cue_dep_trial_set{ises};
    if isempty(temp_data); before_rev_set(ises) = nan; after_rev_set(ises) = nan; continue; end
    temp_step_num = size(temp_data,2);
    rev_trial_idx = find(trial_idx_rec_set{ises}==idenrev_trial_set(ises));
    before_reversal_num = (rev_trial_idx-(win_size_trial-1)/2)-1;
    after_reversal_num = temp_step_num-before_reversal_num;
    before_rev_set(ises) = before_reversal_num;
    after_rev_set(ises) = after_reversal_num;
end
before_rev = min(before_rev_set);
after_rev = min(after_rev_set);

cue_dep_set_cut = cell(ses_num,1);
for ises = 1:ses_num
    temp_data = cue_dep_trial_set{ises};
    if isempty(temp_data); continue; end
    rev_trial_idx = find(trial_idx_rec_set{ises}==idenrev_trial_set(ises));
    zero_idx = (rev_trial_idx-(win_size_trial-1)/2);
    cue_dep_set_cut{ises} = temp_data(:,zero_idx-before_rev:zero_idx+after_rev-1);
end
cue_dep_set_cut_mat = cell2mat(cue_dep_set_cut);
neuron_num_sum = size(cue_dep_set_cut_mat,1);

f_trial=figure('PaperUnits','Centimeters','PaperPosition',[2 2 10.5 11]);
ax_trial = axes;
ax_trial.FontSize = fontsize_ax;
hold on
rev_trial_idx = before_rev+1;
[~,sort_idx] = sort(mean(cue_dep_set_cut_mat(:,rev_trial_idx-round(mean(idenrev_trial_set-thr_pass)):(rev_trial_idx)),2));
for m=1:neuron_num_sum
    m2 = sort_idx(m);
    imagesc((1:size(cue_dep_set_cut_mat(m,:)))-rev_trial_idx+(win_size_trial-1)/2,ones(size(cue_dep_set_cut_mat(m,:)))*m,cue_dep_set_cut_mat(m2,:))
end
colormap(colorscheme)
xlim([1 size(cue_dep_set_cut_mat,2)]-rev_trial_idx+(win_size_trial-1)/2)
xticks(-200:100:400)
xlabel("Trial since reversal",'fontsize',fontsize_label)
ylim([1 neuron_num_sum])
yticks([1 neuron_num_sum])
yticklabels([1 neuron_num_sum])
plot([0 0], [1 neuron_num_sum], 'Color','k','LineWidth',3)
clim_idx = 0.7;
set(gca,'CLim',[-clim_idx clim_idx])
title("Cue preference: delay",'fontsize',fontsize_label)
ylabel("Neurons",'fontsize',fontsize_label)
