clear; close all;

%% Reversal onset response: example neuron
%% 
load('data_VIP_CC\rev_dataname.mat')

ses_num = length(cell_file);

%%
trial_idx_rec_set = cell(1,ses_num);
neuron_num_set = zeros(1,ses_num);
idenrev_trial_set = zeros(1,ses_num);

cue_data_set_1 = cell(1,5);

first_unexp_rwd = zeros(1,ses_num);

%%
i = 10;
% load
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

no_lick = []; no_lick1 = [];

rwded = outcomes{1};
rwded_rec = outcomes{2};
unrwded = outcomes{3};
unrwded_rec = outcomes{4};

puned = outcomes{5};
puned_rec = outcomes{6};
unpuned = outcomes{7};
unpuned_rec = outcomes{8};

epoch_cut_trial = 100;
before_rev = (trial_idx_rec<idenrev_trial)&(trial_idx_rec>=(idenrev_trial-epoch_cut_trial));
during_rev = (trial_idx_rec<idenrev_trial+epoch_cut_trial)&(trial_idx_rec>=(idenrev_trial));
after_rev = (trial_idx_rec>=idenrev_trial+150)&(trial_idx_rec<(idenrev_trial+150+epoch_cut_trial));

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

if neuron_num == 0; return; end

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
        temp_period_diff = temp_lick_onset - temp_delay_offset; % millisec
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

trial_num_before = 4;
trial_num_after = 4;

id_csrwpn = Cue==idenrev_cue(1);
type_id_set = {id_csrwpn};

unexp_rwd = find(trial_idx_rec'>=idenrev_trial & ismember(trial_idx_rec,rwded)',1);
unexp_pun = find(trial_idx_rec'>=idenrev_trial & ismember(trial_idx_rec,puned)',1);
first_unexp_rwd(i) = 2*(unexp_rwd<unexp_pun)-1; % rwd 1 pun -1
unexp_first = min([unexp_pun unexp_rwd]);
unexp_set = [unexp_first unexp_rwd unexp_pun];

for irev = 1
    id_set = cell(length(type_id_set),trial_num_after+1);
    for itype = 1:length(type_id_set)
        temp_type_id = type_id_set{itype};
        temp_id1 = find(trial_idx_rec'<idenrev_trial & temp_type_id,trial_num_before,'last');
        temp_id = find(trial_idx_rec'>unexp_set(irev) & temp_type_id,trial_num_after);
        id_set{itype,1} = temp_id1;
        for iid = 1:trial_num_after
            id_set{itype,iid+1} = temp_id(iid);
        end
    end

    for itype = 1
        for iid = 1:size(id_set,2)
            temp_id = id_set{itype,iid};
            icell = 6;
            % for each neuron
            temp_C_raw_smooth_1 = movmean(C_raw(neuron_id(icell),:),30); % 1.0 s smooth
            
            temp_baseline_smooth_1 = cell(1,length(baseline_onset));
            for j3 = 1:length(baseline_onset)
                temp_baseline_smooth_1{j3} = mean(temp_C_raw_smooth_1(baseline_onset(j3):baseline_offset(j3)));
            end
            temp_baseline_mat_smooth_1 = cell2mat(temp_baseline_smooth_1);
            temp_C_raw_z_smooth_1 = (temp_C_raw_smooth_1-mean(temp_baseline_mat_smooth_1))/std(temp_baseline_mat_smooth_1); % baseline z
            
            trial_trace_1 = nan(length(trial_idx_rec),5*30+1);
            for itrial = 1:length(trial_idx_rec)
                trial_trace_1(itrial,:) = (temp_C_raw_z_smooth_1(state_frame_num_rec(itrial,2)-0.5*30:state_frame_num_rec(itrial,2)+4.5*30));
            end
            
            if iid == 1
                cue_data_set_1_sem = std(trial_trace_1(temp_id,:),[],1)/sqrt(length(temp_id));
            end
            cue_data_set_1{iid} = mean(trial_trace_1(temp_id,:),1);
        end
    end
end


%%
cmap_rw = [186 190 247; 186 190 247; 117 125 240; 30 45 232; 15 23 138; 7 11 69]/255;
cmap_pn = [255 179 193; 255 179 193; 255 102 130; 255 30 70; 153 0 28; 51 0 9]/255;
cmap_nt = [0.8 0.8 0.8; 0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2; 0 0 0];
cmap_cs = [236 179 255; 236 179 255; 210 77 255; 191 0 255; 115 0 153; 56 0 77]/255;

cmap_set = {cmap_pn cmap_rw cmap_nt cmap_cs};
plot_x_cue = -0.5:1/30:4.5;
line_width = 2;

linestyle_set = {':','-','-','-','-','-'};

animal_id_set = {1:length(first_unexp_rwd),find(first_unexp_rwd==1),find(first_unexp_rwd==-1)};

%%
temp_ses_id = animal_id_set{1};

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 5]);

irev = 1;
itype = 1;

temp_cmap = cmap_set{itype};
temp_cmap_before = cmap_set{3-itype};

ianimal = 11;
icell = 6;

ax_fontsize = 10;
label_fontsize = 14;

ax = axes;
ax.FontSize = ax_fontsize;
hold on
fill_x = [plot_x_cue fliplr(plot_x_cue)];
fill_y = [cue_data_set_1{1}+cue_data_set_1_sem fliplr(cue_data_set_1{1}-cue_data_set_1_sem)];
fill_cue = fill(fill_x,fill_y,temp_cmap_before(4,:));
fill_cue.EdgeColor = 'none';
alpha(fill_cue, 0.1);
for iid = 1:size(temp_cmap,1)-1
    temp_mean = cue_data_set_1{iid};
    if iid == 1
        plot(plot_x_cue, temp_mean, 'color', temp_cmap_before(4,:), 'linestyle', linestyle_set{iid}, 'linewidth', line_width);
    else
        plot(plot_x_cue, temp_mean, 'color', temp_cmap(iid,:), 'linestyle', linestyle_set{iid}, 'linewidth', line_width);
    end
end

axisVal = axis;
ylim([axisVal(3) axisVal(4)]); xlim([-0.5 4.5-1/30]);
plot([0 0],[-10 10],'--','color',[0.5 0.5 0.5]); plot([1.5 1.5],[-10 10],'--','color',[0.5 0.5 0.5]); plot([3.5 3.5],[-10 10],'--','color',[0.5 0.5 0.5]);
plot([-0.5 12.5],[0 0],'--','color',[0.5 0.5 0.5]);
xticks(0:4)
xlabel({'Time from cue onset (s)'},'fontsize',label_fontsize)
ylabel('Normalized activity','fontsize',label_fontsize)
set(gca,'xcolor','k','ycolor','k')
yticks(-2:1)
