clear; close all;

%% Reversal onset response
%%
load('data_VIP_CC\rev_dataname.mat')

ses_num = length(cell_file);

%%
trial_idx_rec_set = cell(1,ses_num);
neuron_num_set = zeros(1,ses_num);
idenrev_trial_set = zeros(1,ses_num);

delay_bar_set = cell(7,ses_num,3);

first_unexp_rwd = zeros(1,ses_num);

%%
for i = 1:ses_num
    %% load
    clearvars -except i cell_file beh_file ses_num aligned_file rev_trial_anova_win active_cell_S neurons neuron_drop thr_pass ...
        trial_idx_rec_set neuron_num_set idenrev_trial_set ...
        first_unexp_rwd delay_bar_set
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
    
    trial_num = 20;
    
    id_csrwpn = Cue==idenrev_cue(1);
    id_cspnrw = Cue==idenrev_cue(2);
    id_csrev = id_csrwpn|id_cspnrw;
    type_id_set = {id_csrwpn id_cspnrw};
    
    unexp_rwd = trial_idx_rec(find(trial_idx_rec'>=idenrev_trial & ismember(trial_idx_rec,rwded)',1));
    unexp_pun = trial_idx_rec(find(trial_idx_rec'>=idenrev_trial & ismember(trial_idx_rec,puned)',1));
    first_unexp_rwd(i) = 2*(unexp_rwd<unexp_pun)-1;
    unexp_first = min([unexp_pun unexp_rwd]);
    unexp_set = [unexp_first unexp_rwd unexp_pun];
    
    for ibef = 1:27
        for irev = 1:length(unexp_set)
            id_set = cell(length(type_id_set),1);
            temp_norm_id = cell(1,length(type_id_set));
            for itype = 1:length(type_id_set)
                temp_type_id = type_id_set{itype};
                temp_id_bef = find(trial_idx_rec'<idenrev_trial & temp_type_id,trial_num,'last');
                temp_id_dur = find(trial_idx_rec'>unexp_set(irev) & temp_type_id,trial_num); % 지난 후니까 등호 삭제
                temp_id_af = find(temp_type_id & after_rev',trial_num);
                
                id_set{itype} = [temp_id_bef; temp_id_dur; temp_id_af];
                temp_norm_id{itype} = find(id_csrev & before_rev');
            end
            
            for itype = 1:length(type_id_set)
                temp_id_before_norm = temp_norm_id{itype};
                temp_id = id_set{itype};
                temp_delay_bar = nan(neuron_num,trial_num*3);
                for icell = 1:neuron_num
                    % for each neuron
                    temp_C_raw = C_raw(neuron_id(icell),:);
                    
                    trial_delay = nan(length(trial_idx_rec),1);
                    for itrial = 1:length(trial_idx_rec)
                        trial_delay(itrial) = mean(temp_C_raw(state_frame_num_rec(itrial,3):state_frame_num_rec(itrial,3)+1.5*30));
                    end
                    
                    trial_delay_z = (trial_delay - mean(trial_delay(temp_id_before_norm)))/std(trial_delay(temp_id_before_norm));
                    
                    temp_delay_bar(icell,:) = trial_delay_z(temp_id,:)';
                end
                delay_bar_set{itype,i,irev} = temp_delay_bar;
            end
        end
    end
end

%%
cmap_cs = [159 50 205; 91 29 118]/255;

animal_id_set = {1:length(first_unexp_rwd),find(first_unexp_rwd==1),find(first_unexp_rwd==-1)};

%%
for ianimal = 1:length(animal_id_set)
    f_bar = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 9]);
    
    temp_ses_id = animal_id_set{ianimal};
    irev = 1;
        
    temp_data_set_period = delay_bar_set(:,temp_ses_id,irev);
    
    %%
    ax = axes;
    hold on
    
    ax_fontsize = 20;
    label_fontsize = 24;
    ax.FontSize = ax_fontsize;
    set(gca,'xcolor','k','ycolor','k')
    
    temp_bar_data_set1 = cell2mat(permute(temp_data_set_period(1,:),[2 1]));
    temp_bar_data_set2 = cell2mat(permute(temp_data_set_period(2,:),[2 1]));
    temp_bar_data_set = (temp_bar_data_set1+temp_bar_data_set2)/2;
    
    temp_cmap = cmap_cs;
    
    x_plot_val = [-trial_num:-1 1:trial_num 151:150+trial_num];
    trial_num_fig = 4;
    x_plot = 1:trial_num_fig*2;
    x_plot_label = [-trial_num_fig:-1 1:trial_num_fig];
    x_plot_data = ismember(x_plot_val,x_plot_label);
    temp_bar_data_set = temp_bar_data_set(:,x_plot_data);
    
    if ianimal == 1 && ibef==9
        x_idx = find(x_plot_data);
        x_idx = x_idx(5);
        all_animal_first_set{1} = temp_bar_data_set1(:,x_idx);
        all_animal_first_set{2} = temp_bar_data_set2(:,x_idx);
    end
    if ianimal == 2
        divided_animal_first_set{1} = temp_bar_data_set(:,5);
    end
    if ianimal == 3
        divided_animal_first_set{2} = temp_bar_data_set(:,5);
    end
    
    e1 = errorbar(x_plot,mean(temp_bar_data_set),std(temp_bar_data_set,[],1)/sqrt(size(temp_bar_data_set,1)));
    e1.Color = temp_cmap(1,:); e1.CapSize = 0; e1.LineStyle = 'none';
    e1.Marker = 'o'; e1.LineWidth = 1;
    e1.MarkerFaceColor = temp_cmap(2,:); e1.MarkerEdgeColor = temp_cmap(2,:);
    
    xticks(x_plot)
    xticklabels(x_plot_label)
    
    xlim([0 size(temp_bar_data_set,2)+1])
    ylabel('Normalized activity','fontsize',label_fontsize)
    xlabel('Cue appearance','fontsize',label_fontsize)
    xtickangle(45)
    
    axisVal = axis;
    ylim(axisVal(3:4))
    plot([trial_num_fig trial_num_fig]+0.5,axisVal(3:4),'--','color',[0.5 0.5 0.5])
    plot([0 size(temp_bar_data_set,2)+1],[0 0],'-k')
    
    [h,p] = ttest(temp_bar_data_set);
    h(p<0.01) = 2; h(p<0.001) = 3;
    t = gobjects(1,size(temp_bar_data_set,2));
    for ipoint = 1:size(temp_bar_data_set,2)
        if h(ipoint)>0
            t(ipoint) = text(x_plot(ipoint),axisVal(4),repmat('*',1,h(ipoint)),'fontsize',18,'horizontalalignment','center','verticalalignment','cap','fontname','arial');
        end
    end
end
