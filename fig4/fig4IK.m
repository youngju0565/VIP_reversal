clear; close all;

%% RP reversal: Outcome-dependent (Scatter)
%%
load('data_VIP_CC\rev_dataname.mat')

ses_num = length(cell_file);

%%
out_mean_set = cell(ses_num,3); % 4 type x epochs

trial_idx_rec_set = cell(1,ses_num);
neuron_num_set = zeros(1,ses_num);
idenrev_trial_set = zeros(1,ses_num);

%%
for i = 1:ses_num
    %% load
    clearvars -except ifile i cell_file beh_file ses_num aligned_file rev_trial_anova_win neurons neuron_drop thr_pass gpio_file ...
        out_mean_set ...
        trial_idx_rec_set neuron_num_set idenrev_trial_set
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
    
    epoch_cut_trial = 75;
    before_rev = (trial_idx_rec<idenrev_trial)&(trial_idx_rec>=(idenrev_trial-epoch_cut_trial));
    during_rev = (trial_idx_rec<idenrev_trial+epoch_cut_trial)&(trial_idx_rec>=(idenrev_trial));
    after_rev = (trial_idx_rec>=idenrev_trial+100)&(trial_idx_rec<(idenrev_trial+100+epoch_cut_trial));
    
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
        rw_period_set_time(j_rec) = stateTime(j,4)+temp_period_diff;
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
    
    epoch_set = {before_rev, during_rev, after_rev};
    for iepoch = 1:3
        temp_epoch = epoch_set{iepoch}';

        out_mean = zeros(neuron_num,4);
        
        for icell = 1:neuron_num
            temp_C_raw = C_raw(neuron_id(icell),:);
            temp_C_raw_smooth = movmean(C_raw(neuron_id(icell),:),15); % 0.5 s smooth

            temp_baseline = cell(1,length(baseline_onset));
            for j3 = 1:length(baseline_onset)
                temp_baseline{j3} = mean(temp_C_raw(baseline_onset(j3):baseline_offset(j3)));
            end
            temp_baseline_mat = cell2mat(temp_baseline);
            temp_C_raw_z = (temp_C_raw-mean(temp_baseline_mat))/std(temp_baseline_mat); % baseline z
            
            trial_cut_out = zeros(length(trial_idx_rec),1);
            for itrial = 1:length(trial_idx_rec)
                trial_cut_out(itrial) = mean(temp_C_raw_z(rw_period_set_frame(itrial):rw_period_set_frame(itrial)+1.5*30));
            end
            
            % mean activity per cues/types
            out_mean(icell,1) = mean(trial_cut_out(Ident==2&Prob==75&Reward==1&temp_epoch)); % Rw-Rwded
            out_mean(icell,2) = mean(trial_cut_out(Ident==2&Prob==75&Reward==0&temp_epoch)); % Rw-omit
            out_mean(icell,3) = mean(trial_cut_out(Ident==3&Reward==1&temp_epoch)); % Pn-puned
            out_mean(icell,4) = mean(trial_cut_out(Ident==3&Reward==0&temp_epoch)); % Pn-omit
        end
        
        out_mean_set{i,iepoch} = out_mean;
    end
end

%%
scatter_data_set = cell(1,8);

empty_idx = cellfun(@isempty,out_mean_set(:,1));

scatter_data_set{1} = cell2mat(cellfun(@(x) x(:,1)-x(:,2),out_mean_set(~empty_idx,:),'uniformoutput',false)); % Rw dep
scatter_data_set{2} = cell2mat(cellfun(@(x) x(:,3)-x(:,4),out_mean_set(~empty_idx,:),'uniformoutput',false)); % Pn dep

fontsize_ax = 20;
fontsize_label = 24;

%% Reward dependent
iscatter = 1;
temp_scatter_data = scatter_data_set{iscatter};
temp_data = temp_scatter_data(:,[1 3]);

f_scatter = figure('PaperUnits','Centimeters','PaperPosition',[2 2 11 11]);
ax_scatter = axes;
hold on
ax_scatter.FontSize = fontsize_ax;
set(gca,'xcolor','k','ycolor','k')

data_draw = temp_data; x = 1; y = 2;
marker_color = 0.4;
scatter(data_draw(:,x), data_draw(:,y),'filled','markeredgecolor',[1 1 1]*marker_color,'markerfacecolor',[1 1 1]*marker_color)
xlim([-2 2])
ylim([-2 2])
xticks(-2:2)
yticks(-2:2)
xline1 = line([-2000 2000],[0 0]);
xline1.Color = 'k'; xline1.LineStyle = '--';
yline1 = line([0 0],[-2000 2000]);
yline1.Color = 'k'; yline1.LineStyle = '--';
fit1 = polyfit(data_draw(:,x),data_draw(:,y),1);
plot([min(data_draw(:,x)) max(data_draw(:,x))]*1.2,polyval(fit1,[min(data_draw(:,x)) max(data_draw(:,x))]*1.2),'Color','k','linewidth',1.5)
[r,p] = corrcoef(data_draw(:,x),data_draw(:,y));
str1 = {['{\itr} = ',num2str(r(1,2),'%.3f')],['{\itp} = ',num2str(p(1,2),'%.3f')]};
annotation('textbox',[0.7 0.85 0.15 0.1],'String',str1,'FitBoxToText','on','edgecolor','none','fontsize',14)
annotation('textbox',[0.2 0.85 0.15 0.1],'String','Reward','FitBoxToText','on','edgecolor','none','fontsize',24)
xlabel('Before','fontsize',fontsize_label)
ylabel('After','fontsize',fontsize_label)

%% Punishment dependent
iscatter = 2;
temp_scatter_data = scatter_data_set{iscatter};

[~,max_during_pn_id] = max(temp_scatter_data(:,2));
[~,min_after_pn_id] = min(temp_scatter_data(:,3));

temp_data = temp_scatter_data(:,[1 3]);
temp_data([max_during_pn_id min_after_pn_id],:) = [];

f_scatter2 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 11 11]);
ax_scatter = axes;
hold on
ax_scatter.FontSize = fontsize_ax;
set(gca,'xcolor','k','ycolor','k')

data_draw = temp_data; x = 1; y = 2;
marker_color = 0.4;
scatter(data_draw(:,x), data_draw(:,y),'filled','markeredgecolor',[1 1 1]*marker_color,'markerfacecolor',[1 1 1]*marker_color)
xlim([-2 2])
ylim([-2 2])
xticks(-2:2)
yticks(-2:2)
xline1 = line([-2000 2000],[0 0]);
xline1.Color = 'k'; xline1.LineStyle = '--';
yline1 = line([0 0],[-2000 2000]);
yline1.Color = 'k'; yline1.LineStyle = '--';
fit1 = polyfit(data_draw(:,x),data_draw(:,y),1);
plot([min(data_draw(:,x)) max(data_draw(:,x))]*1.2,polyval(fit1,[min(data_draw(:,x)) max(data_draw(:,x))]*1.2),'Color','k','linewidth',1.5)
[r,p] = corrcoef(data_draw(:,x),data_draw(:,y));
str1 = {['{\itr} = ',num2str(r(1,2),'%.3f')],['{\itp} = ',num2str(p(1,2),'%.3f')]};
annotation('textbox',[0.7 0.85 0.15 0.1],'String',str1,'FitBoxToText','on','edgecolor','none','fontsize',14)
annotation('textbox',[0.2 0.85 0.15 0.1],'String','Punishment','FitBoxToText','on','edgecolor','none','fontsize',24)
xlabel('Before','fontsize',fontsize_label)
ylabel('After','fontsize',fontsize_label)
