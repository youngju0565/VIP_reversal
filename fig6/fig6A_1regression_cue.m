clear; close all; warning off

%% Multiple linear regression
%%
load('data_VIP_CC\rev_dataname.mat')

cmap2 = [0.6350 0.0780 0.1840; 232/255 126/255 58/255; 30/255 45/255 232/255; 255/255 30/255 70/255; 0.5 0.5 0.5; 0 0 0];
ses_num = length(cell_file);
cpd_filename_set = cell(1,ses_num);

%%
for i=1:ses_num
    %% load
    clearvars -except i cell_file beh_file ses_num aligned_file active_cell_S cmap2 thr_pass neuron_drop cpd_filename_set
    load(cell_file{i},'neuron');
    load(beh_file{i});
    load(aligned_file{i});
    
    %% recorded trial types
    trial_idx_rec_original = trial_idx_rec;

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
    
    no_lick = []; no_lick1 = [];
    
    epoch_cut_trial = 100;
    before_rev = (trial_idx_rec<idenrev_trial)&(trial_idx_rec>=(idenrev_trial-epoch_cut_trial));
    during_rev = (trial_idx_rec<idenrev_trial+epoch_cut_trial)&(trial_idx_rec>=(idenrev_trial));
    after_rev = (trial_idx_rec>=idenrev_trial+150)&(trial_idx_rec<(idenrev_trial+150+epoch_cut_trial));
    
    epoch_set = {before_rev, during_rev, after_rev};
    
    %% predictors for GLM (no lick num)
    Prob = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    Ident = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    Reward = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
    Punish = zeros(length(trial_idx_rec)-length(no_lick)-length(no_lick1), 1);
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
        if id_l == 2
            Reward(temp_itrial1) = rewarded_l; % 0 1
        end
        if id_l == 3
            Punish(temp_itrial1) = rewarded_l; 
        end
        temp_itrial1 = temp_itrial1+1;
    end
    
    Cue_Rw = double(Ident==2&Prob==75);
    Cue_Pn = double(Ident==3&Prob==75);
    Cue_Rw_1 = [0; Cue_Rw(1:end-1)];
    Cue_Pn_1 = [0; Cue_Pn(1:end-1)];
    Reward_1 = [0; Reward(1:end-1)];
    Punish_1 = [0; Punish(1:end-1)];
    Lick = zeros(size(Reward));
    Speed = zeros(size(Reward));
    GLM_T = table(Cue_Rw, Cue_Pn, Reward, Punish, Lick, Speed, Cue_Rw_1, Cue_Pn_1, Reward_1, Punish_1);

    %% response variable
    win_size = 1.5; % s

    state_frame_num_rec = state_frame_num;
    
    rw_period_set_frame = zeros(length(trial_idx_rec),1);
    rw_period_set_time = zeros(length(trial_idx_rec),1);
    period_diff_time = zeros(length(trial_idx_rec),1);
    j_rec = 0;
    rwded = outcomes{1};
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
        rw_period_set_time(j_rec) = stateTime(j,4)+temp_period_diff;
    end
    
    lick_out_post = zeros(length(trial_idx_rec),1);
    speed_out_post = zeros(length(trial_idx_rec),1);
    for itrial = 1:length(trial_idx_rec)
        out_onset = rw_period_set_time(itrial);
        lick_out_post(itrial) = length(find((lickTime(:,1)>=out_onset) ...
            & (lickTime(:,1)<(out_onset + win_size*1000))));
        speed_out_post(itrial) = length(find((cylinderTime(:,1)>=out_onset) ...
            & (cylinderTime(:,1)<(out_onset + win_size*1000))));
    end
    
    C_raw = neuron.C_raw;
    S = neuron.S;
    event_thr_set = [0.015 0.025 0.020];
    hz_cut = 2;
    neuron_num = size(C_raw,1);
    S_rate = zeros(1,neuron_num);
    for icell = 1:neuron_num
        S_rate(icell) = sum(S(icell,:)>0)/(size(S,2)/30);
    end
    neuron_id = find(S_rate>event_thr_set(hz_cut));
    neuron_id(ismember(neuron_id, neuron_drop{i}))=[];
    neuron_num = length(neuron_id);
    
    if neuron_num == 0; disp('No cells'); continue; end
    
    glm_y = zeros(neuron_num,length(Reward));
    for icell3 = 1:neuron_num % for every neuron
        temp_C_raw = C_raw(neuron_id(icell3),:);
        
        temp_itrial1 = 1;
        for itrial1 = 1:length(trial_idx_rec) % for each recorded trial
            rec_itrial = trial_idx_rec(itrial1);
            if ismember(rec_itrial, no_lick)||ismember(rec_itrial,no_lick1); continue; end
            start_frame_m2 = rw_period_set_frame(itrial1);
            glm_y(icell3,temp_itrial1) = mean(temp_C_raw((start_frame_m2):(start_frame_m2+win_size*30)));
            
            temp_itrial1 = temp_itrial1+1;
        end
    end
        
    %% Regression setting
    x_var_num = size(GLM_T,2);
    
    %% save settings
    session_name = strsplit(beh_file{i},'_');
    animal_name = session_name{1};
    session_date = session_name{2};
    if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
    else; win_size_str = ''; end
    if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
    else; hz_cut_str = ''; end
    if epoch_cut_trial~=100; trial_str = [' ',num2str(epoch_cut_trial),'trials'];
    else; trial_str = ''; end
    fig_dir = ['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\',...
        animal_name,'\',animal_name,'_',session_date,'\regression epochs outcome' win_size_str hz_cut_str trial_str ' 230627\'];
    mkdir(fig_dir);
    
    %% Regression
    fprintf("fitlm start (CPD): 3 epochs\n");
    reg_start = tic;
    % storage
    CPD = zeros(neuron_num,x_var_num,3);
    P_value = zeros(neuron_num,x_var_num,3);
    SRC = zeros(neuron_num,x_var_num,3);
    T_value = zeros(neuron_num,x_var_num,3);

    % regression: whole
    for iepoch = 1:3
        temp_epoch = epoch_set{iepoch};
        
        for icell3 = 1:neuron_num
            tic; fprintf('fitting cell no %d of %d cells\n',icell3, neuron_num);
            
            tmp_T = GLM_T(temp_epoch,:);
            tmp_T.Lick = lick_out_post(temp_epoch);
            tmp_T.Speed = speed_out_post(temp_epoch);
            tmp_T.Y = glm_y(icell3,temp_epoch)';
            
            tmp_mdl = fitlm(tmp_T);
            
            % FON
            P_value(icell3,:,iepoch) = tmp_mdl.Coefficients.pValue(2:end);
            % SRC
            std_x = std(tmp_T{:,1:end-1})';
            std_y = std(tmp_T.Y);
            SRC(icell3,:,iepoch) = tmp_mdl.Coefficients.Estimate(2:end).*std_x/std_y;
            % t-value
            T_value(icell3,:,iepoch) = tmp_mdl.Coefficients.tStat(2:end);
            
            % CPD
            tmp_SSE = tmp_mdl.SSE;
            for ivar = 1:x_var_num
                tmp_T_var = tmp_T;
                tmp_T_var(:,ivar) = [];
                tmp_mdl_var = fitlm(tmp_T_var);
                var_SSE = tmp_mdl_var.SSE;
                CPD(icell3,ivar,iepoch) = (var_SSE - tmp_SSE)/var_SSE;
            end
            
            toc
        end
    end
    disp("Regression done"); toc(reg_start);

    %% save workspace
    clear f1_cpd f1_cpd_1 fill1a plots plots_1  
    save([fig_dir, '\variables','_win_',num2str(win_size),'.mat'],...
        'tmp_mdl','tmp_mdl_var','GLM_T','lick_out_post','speed_out_post','glm_y','neuron_num','iter_num','win_size')
    save([fig_dir, '\CPD','_win_',num2str(win_size),'_draw_ready.mat'],...
        'CPD','neuron_num','win_size')
    save([fig_dir, '\FON','_win_',num2str(win_size),'_draw_ready.mat'],...
        'P_value','neuron_num','win_size')
    save([fig_dir, '\SRC','_win_',num2str(win_size),'_draw_ready.mat'],...
        'SRC','neuron_num','win_size')
    save([fig_dir, '\t_value','_win_',num2str(win_size),'_draw_ready.mat'],...
        'T_value','neuron_num','win_size')
    
    cpd_filename_set{i} = [fig_dir, '\CPD','_win_',num2str(win_size),'_draw_ready.mat'];
    
end

%%
warning on
unify_dir = ['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome'...
    win_size_str ' n' num2str(ses_num)  hz_cut_str trial_str ' 230627\'];
mkdir(unify_dir)
save([unify_dir 'CPD_set.mat'],'cpd_filename_set');
