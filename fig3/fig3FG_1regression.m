clear; close all;

%% Multiple linear regression
%%
load('data_VIP_CC\no_rev_dataname.mat')

ses_num = length(cell_file);

src_filename_set = cell(1,ses_num);

%%
for i=1:ses_num
    %% load
    warning off
    clearvars -except i cell_file beh_file ses_num aligned_file thr_pass neuron_drop src_filename_set
    load(cell_file{i},'neuron');
    load(beh_file{i});
    load(aligned_file{i});
    
    fprintf('Start %d / %d animal.\n',i,ses_num);
    
    %% recorded trial types
    % use only after thr pass    
    [types, cues, types_rec, no_lick, no_lick1, outcomes] = recorded_trial_types(odorCue, waterReward, outcomeIdentity, stateTime, lickTime, trial_idx_rec, nTrial);
    for itype = 1:8
        eval(sprintf('type%d = types{%d};',itype,itype));
        eval(sprintf('type%d_rec = types_rec{%d};',itype,itype));
    end
    for icue = 1:4
        eval(sprintf('cue%d = cues{%d};',icue,icue));
    end

    no_lick = []; no_lick1 = [];
    
    epoch_cut_trial = 360; 
    after_pass = (trial_idx_rec<thr_pass(i)+epoch_cut_trial)&(trial_idx_rec>=(thr_pass(i)));
        
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
    C_raw = neuron.C_raw;
    S = neuron.S;
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
    
    if neuron_num == 0; disp('No cells'); continue; end
    
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
            & (lickTime(:,1)<(out_onset + win_size*1000))));
        speed_out_post(itrial) = length(find((cylinderTime(:,1)>=out_onset) ...
            & (cylinderTime(:,1)<(out_onset + win_size*1000))));
    end
    
    glm_y = zeros(neuron_num,length(Reward));
    for icell3 = 1:neuron_num
        temp_C_raw = C_raw(neuron_id(icell3),:);
        
        temp_itrial1 = 1;
        for itrial1 = 1:length(trial_idx_rec)
            rec_itrial = trial_idx_rec(itrial1);
            if ismember(rec_itrial, no_lick)||ismember(rec_itrial,no_lick1); continue; end
            start_frame_m2 = rw_period_set_frame(itrial1);
            glm_y(icell3,temp_itrial1) = mean(temp_C_raw((start_frame_m2):(start_frame_m2+win_size*30)));
            
            temp_itrial1 = temp_itrial1+1;
        end
    end
    
    %% Regression setting
    % cue=3 reference category
    x_var_num = size(GLM_T,2);
    
    %% save settings
    session_name = strsplit(beh_file{i},'_');
    animal_name = session_name{1};
    session_date = session_name{2};
    if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
    else; win_size_str = ''; end
    if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
    else; hz_cut_str = ''; end
    fig_dir = ['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\',...
        animal_name,'\',animal_name,'_',session_date,'\regression norev outcome term 231108' win_size_str hz_cut_str '\'];
    mkdir(fig_dir);
    
    %% Regression
    fprintf("fitlm start\n");
    reg_start = tic;
    % storage
    CPD = zeros(neuron_num,x_var_num);
    P_value = zeros(neuron_num,x_var_num);
    SRC = zeros(neuron_num,x_var_num);
    T_value = zeros(neuron_num,x_var_num);
    
    CPD_rwd = zeros(neuron_num,x_var_num);
    P_value_rwd = zeros(neuron_num,x_var_num);
    SRC_rwd = zeros(neuron_num,x_var_num);
    T_value_rwd = zeros(neuron_num,x_var_num);
    
    CPD_pun = zeros(neuron_num,x_var_num);
    P_value_pun = zeros(neuron_num,x_var_num);
    SRC_pun = zeros(neuron_num,x_var_num);
    T_value_pun = zeros(neuron_num,x_var_num);
    
    % regression: whole
    temp_epoch = after_pass;
    for icell3 = 1:neuron_num
        tic; fprintf('fitting cell no %d of %d cells\n',icell3, neuron_num);
        
        tmp_T = GLM_T(temp_epoch,:);
        tmp_T.Lick = lick_out_post(temp_epoch);
        tmp_T.Speed = speed_out_post(temp_epoch);
        tmp_T.Y = glm_y(icell3,temp_epoch)';
        tmp_mdl = fitlm(tmp_T);

        % FON
        P_value(icell3,:) = tmp_mdl.Coefficients.pValue(2:end);
        % SRC
        std_x_rwd = std(tmp_T{:,1:end-1})';
        std_y_rwd = std(tmp_T.Y);
        SRC(icell3,:) = tmp_mdl.Coefficients.Estimate(2:end).*std_x_rwd/std_y_rwd;
        % t-value
        T_value(icell3,:) = tmp_mdl.Coefficients.tStat(2:end);
        % CPD
        tmp_SSE = tmp_mdl.SSE;
        for ivar = 1:x_var_num
            tmp_T_var = tmp_T;
            tmp_T_var(:,ivar) = [];
            tmp_mdl_var = fitlm(tmp_T_var);
            var_SSE = tmp_mdl_var.SSE;
            CPD(icell3,ivar) = (var_SSE - tmp_SSE)/var_SSE;
        end
        
        % reward term: only reward cue trials
        tmp_cue_rw = Cue_Rw(temp_epoch);
        tmp_T_rwd = tmp_T(tmp_cue_rw==1,:);
        tmp_mdl_rwd = fitlm(tmp_T_rwd);
        % FON
        P_value_rwd(icell3,:) = tmp_mdl_rwd.Coefficients.pValue(2:end);
        % SRC
        std_x_rwd = std(tmp_T_rwd{:,1:end-1})';
        std_y_rwd = std(tmp_T_rwd.Y);
        SRC_rwd(icell3,:) = tmp_mdl_rwd.Coefficients.Estimate(2:end).*std_x_rwd/std_y_rwd;
        % t-value
        T_value_rwd(icell3,:) = tmp_mdl_rwd.Coefficients.tStat(2:end);
        % CPD
        tmp_SSE = tmp_mdl_rwd.SSE;
        for ivar = 3
            tmp_T_var = tmp_T_rwd;
            tmp_T_var(:,ivar) = [];
            tmp_mdl_var = fitlm(tmp_T_var);
            var_SSE = tmp_mdl_var.SSE;
            CPD_rwd(icell3,ivar) = (var_SSE - tmp_SSE)/var_SSE;
        end
            
        % punish term: only punish cue trials
        tmp_cue_pn = Cue_Pn(temp_epoch);
        tmp_T_pun = tmp_T(tmp_cue_pn==1,:);
        tmp_mdl_pun = fitlm(tmp_T_pun);
        % FON
        P_value_pun(icell3,:) = tmp_mdl_pun.Coefficients.pValue(2:end);
        % SRC
        std_x_pun = std(tmp_T_pun{:,1:end-1})';
        std_y_pun = std(tmp_T_pun.Y);
        SRC_pun(icell3,:) = tmp_mdl_pun.Coefficients.Estimate(2:end).*std_x_pun/std_y_pun;
        % t-value
        T_value_pun(icell3,:) = tmp_mdl_pun.Coefficients.tStat(2:end);
        % CPD
        tmp_SSE = tmp_mdl_pun.SSE;
        for ivar = 4
            tmp_T_var = tmp_T_pun;
            tmp_T_var(:,ivar) = [];
            tmp_mdl_var = fitlm(tmp_T_var);
            var_SSE = tmp_mdl_var.SSE;
            CPD_pun(icell3,ivar) = (var_SSE - tmp_SSE)/var_SSE;
        end
        
        toc
    end
    disp("Regression done"); toc(reg_start);
    
    %% save worksapce
    save([fig_dir, '\variables','_win_',num2str(win_size),'.mat'],...
        'tmp_mdl_rwd','tmp_mdl_var','tmp_mdl_rwd','tmp_mdl_pun','GLM_T','lick_out_post','speed_out_post','glm_y','neuron_num','iter_num','win_size')
    save([fig_dir, '\CPD','_win_',num2str(win_size),'_draw_ready.mat'],...
        'CPD','CPD_rwd','CPD_pun','neuron_num','win_size')
    save([fig_dir, '\FON','_win_',num2str(win_size),'_draw_ready.mat'],...
        'P_value','P_value_rwd','P_value_pun','neuron_num','win_size')
    save([fig_dir, '\SRC','_win_',num2str(win_size),'_draw_ready.mat'],...
        'SRC','SRC_rwd','SRC_pun','neuron_num','win_size')
    save([fig_dir, '\t_value','_win_',num2str(win_size),'_draw_ready.mat'],...
        'T_value','T_value_rwd','T_value_pun','neuron_num','win_size')
    
    src_filename_set{i} = [fig_dir, '\SRC','_win_',num2str(win_size),'_draw_ready.mat'];
    
end
   
%%
warning on
unify_dir = ['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression norev outcome term 231108'...
    win_size_str ' n' num2str(ses_num)  hz_cut_str '\'];
mkdir(unify_dir)
save([unify_dir 'SRC_set.mat'],'src_filename_set','epoch_cut_trial');
