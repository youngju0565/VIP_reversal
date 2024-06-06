clear; close all;

%% RP reversal: Example neurons
%%
load('data_VIP_CC\rev_dataname.mat')

ses_num = length(cell_file);

colorscheme = othercolor('BuDRd_18');
close all;

animal_set = [4 5];
cell_set = [2 5];

%%
for iset = 1:2
    %% load
    clearvars -except iset cell_file beh_file ses_num aligned_file rev_trial_anova_win active_cell_S neurons neuron_drop thr_pass ...
        colorscheme animal_set cell_set
    i = animal_set(iset);
    load(cell_file{i},'neuron');
    load(beh_file{i});
    load(aligned_file{i});
    
    %% recorded trial types
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;
    if outcomeIdentity(1,idenrev_cue(1))==3
        idenrev_cue = flipud(idenrev_cue);
    end
    
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
    
    %%
    baseline_onset = state_frame_num(:,1);
    baseline_offset = state_frame_num(:,2)-1;
    
    Cue_1 = (Cue==idenrev_cue(1)); % Rwd->Pun
    Cue_2 = (Cue==idenrev_cue(2)); % Pun->Rwd
    
    ax_pos_set = [1.95 6.9; 1.95 3.1; 8.8 6.9; 8.8 3.1];
    ax_size = [5.5 3.2]; % cm
    
    fontsize_ax = 20;
    fontsize_label = 24;
    
    smoothing_bin = 5; % frame
    
    title_set = {'Before reversal','After reversal'};
    epoch_set = {before_rev, after_rev};
    
    r_color = [139 23 28]/255;
    b_color = [57 89 168]/255;
        
    %%
    icell = cell_set(iset);
    f_cell = figure('PaperUnits','Centimeters','PaperPosition',[2 2 14 12.3]);
    
    temp_C_raw = C_raw(neuron_id(icell),:);
    temp_baseline = cell(1,length(baseline_onset));
    for j3 = 1:length(baseline_onset)
        temp_baseline{j3} = mean(temp_C_raw(baseline_onset(j3):baseline_offset(j3)));
    end
    temp_baseline_mat = cell2mat(temp_baseline);
    temp_C_raw_z = (temp_C_raw-mean(temp_baseline_mat))/std(temp_baseline_mat);
    
    clf
    ax_heat_rwpn = gobjects(1,2);
    ax_heat_pnrw = gobjects(1,2);
    ax_trace = gobjects(1,2);
    for iepoch = 1:2
        temp_epoch = epoch_set{iepoch}';
        
        rwpn_data = zeros(sum(temp_epoch&Cue_1),5*30+1);
        pnrw_data = zeros(sum(temp_epoch&Cue_2),5*30+1);
        
        rwpn_idx = find(temp_epoch&Cue_1);
        pnrw_idx = find(temp_epoch&Cue_2);
        for itrial = 1:sum(temp_epoch&Cue_1)
            rwpn_data(itrial,:) = temp_C_raw_z(state_frame_num(rwpn_idx(itrial),2)-0.5*30:state_frame_num(rwpn_idx(itrial),2)+4.5*30);
        end
        for itrial = 1:sum(temp_epoch&Cue_2)
            pnrw_data(itrial,:) = temp_C_raw_z(state_frame_num(pnrw_idx(itrial),2)-0.5*30:state_frame_num(pnrw_idx(itrial),2)+4.5*30);
        end
        
        rwpn_data = movmean(rwpn_data,smoothing_bin,2);
        pnrw_data = movmean(pnrw_data,smoothing_bin,2);
        
        % heat map
        ax_heat_rwpn(iepoch) = axes('units','centimeters','position',[ax_pos_set(2*iepoch-1,1) 0.1+ax_pos_set(2*iepoch-1,2)+ax_size(2)*sum(Cue_2&temp_epoch)/sum((Cue_1&temp_epoch)|(Cue_2&temp_epoch)) ax_size(1) ax_size(2)*sum(Cue_1&temp_epoch)/sum((Cue_1&temp_epoch)|(Cue_2&temp_epoch))],'fontsize',fontsize_ax);
        hold on
        xticklabels('')
        title(title_set{iepoch})
        if iepoch==1; ylabel('Trial','fontsize',fontsize_label); end
        imagesc([-0.5 4.5],[1 size(rwpn_data,1)],rwpn_data)
        colormap(colorscheme)
        color_rwpn = caxis;
        xlim([-0.5 4.5])
        ylim([1-0.5 size(rwpn_data,1)+0.5])
        yticks(size(rwpn_data,1)+0.5)
        yticklabels(sum(temp_epoch&Cue_1)+sum(temp_epoch&Cue_2))
        plot([0 0],[-100 100],'--','color','k')
        plot([1.5 1.5],[-100 100],'--','color','k')
        plot([3 3],[-100 100],'--','color','k')
        plot([3.5 3.5],[-100 100],'--','color','k')
        rectangle('position',[-0.68 0.5 0.15 size(rwpn_data,1)],'clipping','off','edgecolor','none','facecolor',b_color)
        
        ax_heat_pnrw(iepoch) = axes('units','centimeters','position',[ax_pos_set(2*iepoch-1,:) ax_size(1) ax_size(2)*sum(Cue_2&temp_epoch)/sum((Cue_1&temp_epoch)|(Cue_2&temp_epoch))],'fontsize',fontsize_ax);
        hold on
        xticklabels('')
        imagesc([-0.5 4.5],[1 size(pnrw_data,1)],pnrw_data)
        colormap(colorscheme)
        color_pnrw = caxis;
        xlim([-0.5 4.5])
        ylim([1-0.5 size(pnrw_data,1)+0.5])
        yticks(0.5)
        yticklabels(1)
        plot([0 0],[-100 100],'--','color','k')
        plot([1.5 1.5],[-100 100],'--','color','k')
        plot([3 3],[-100 100],'--','color','k')
        plot([3.5 3.5],[-100 100],'--','color','k')
        rectangle('position',[-0.68 0.5 0.15 size(pnrw_data,1)],'clipping','off','edgecolor','none','facecolor',r_color)
        
        switch iset
            case 1; clim_val = [-4 3];
            case 2; clim_val = [-1.5 1.5];
        end
        
        set(ax_heat_rwpn(iepoch),'CLim',clim_val); set(ax_heat_pnrw(iepoch),'CLim',clim_val);
        
        % trace
        ax_trace(iepoch) = axes('units','centimeters','position',[ax_pos_set(2*iepoch,:) ax_size],'fontsize',fontsize_ax);
        hold on
        xlabel({'Time from';'cue onset (s)'},'fontsize',fontsize_label)
        if iepoch==1; ylabel('z-score','fontsize',fontsize_label); end
        plot(-0.5:1/30:4.5,mean(rwpn_data),'color',b_color,'linewidth',2)
        plot(-0.5:1/30:4.5,mean(pnrw_data),'color',r_color,'linewidth',2)
        axisVal = axis;
        ylim([axisVal(3) axisVal(4)])
        yticks([axisVal(3) axisVal(4)])
        
        switch iset
            case 1
                ylim([-0.5 2])
                yticks([0 2])
            case 2
                switch iepoch
                    case 1
                        ylim_val = [-0.5 0.5];
                    case 2
                        ylim_val = [-0.5 1];
                end
                ylim(ylim_val)
                yticks([0 ylim_val(2)])
        end
        fill_x = [(-0.5:1/30:4.5) fliplr(-0.5:1/30:4.5)];
        fill_y_rwpn = [mean(rwpn_data)+std(rwpn_data)/sqrt(size(rwpn_data,1)) fliplr(mean(rwpn_data)-std(rwpn_data)/sqrt(size(rwpn_data,1)))];
        fill_rwpn = fill(fill_x,fill_y_rwpn,b_color);
        fill_rwpn.EdgeColor = 'none';
        alpha(fill_rwpn,0.1)
        fill_y_pnrw = [mean(pnrw_data)+std(pnrw_data)/sqrt(size(pnrw_data,1)) fliplr(mean(pnrw_data)-std(pnrw_data)/sqrt(size(pnrw_data,1)))];
        fill_pnrw = fill(fill_x,fill_y_pnrw,r_color);
        fill_pnrw.EdgeColor = 'none';
        alpha(fill_pnrw,0.1)
        xlim([-0.5 4.5])
        xticks([0 1.5 3 3.5])
        xticklabels({'0','1.5','3',''})
        plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5])
        plot([1.5 1.5],[-100 100],'--','color',[0.5 0.5 0.5])
        plot([3 3],[-100 100],'--','color',[0.5 0.5 0.5])
        plot([3.5 3.5],[-100 100],'--','color',[0.5 0.5 0.5])
    end
end