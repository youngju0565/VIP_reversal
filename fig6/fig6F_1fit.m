clear; close all;

%%
rng(2);

load('data_VIP_CC\rev_dataname.mat')
beh_file_rev = beh_file;

load('data_VIP_CC\rev_d1_dataname.mat')

%%
for ianimal = 1:length(beh_file_rev)
    %% transform to eventlog
    load(beh_file_rev{ianimal})
    [idenrev_trial_fig,~] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial_fig = unique(idenrev_trial_fig)+1;
    [eventlog1,IRI_rw1,IRI_pn1] = behmat2eventlog(stateTime,outcomeIdentity,waterReward,...
        nTrial,odorCue,outcomeProbability,lickTime);
    trial_idx1 = 1:nTrial;
    
    load(beh_file_d1{ianimal})
    % only after rev trials
    [idenrev_trial,~] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;
    [eventlog2,IRI_rw2,IRI_pn2] = behmat2eventlog(stateTime(idenrev_trial:nTrial,:),...
        outcomeIdentity(idenrev_trial:nTrial,:),waterReward(idenrev_trial:nTrial,1),...
        length(idenrev_trial:nTrial),odorCue(idenrev_trial:nTrial,1),...
        outcomeProbability(idenrev_trial:nTrial,:),lickTime);
    trial_idx2 = -length(idenrev_trial:nTrial)+1:0;
    
    eventlog = joinEventlogs(eventlog2,eventlog1);
    trial_idx = [trial_idx2 trial_idx1];
    
    %% lick plot
    % 1. one day before
    load(beh_file_d1{ianimal})
    
    lick_diff = diff(lickTime, 1);
    short_lick_diff = find(lick_diff(:,1)<=50)+1;
    lickTime_50ms = lickTime;
    lickTime_50ms(short_lick_diff,:)=[];
    
    lick_timings_50ms = cell(nTrial,1);
    delay_timings = stateTime(:,3:4)-stateTime(:,1)+1;
    delay_licknum_50ms_1 = zeros(size(lickNum));
    for j = idenrev_trial:nTrial % only after rev trials
        jtrial_ind = find(lickTime_50ms(:,2)==j);
        lick_timings_50ms{j} = lickTime_50ms(jtrial_ind,1) - stateTime(j,1)+1;
        delay_licknum_50ms_1(j) = length(find((lick_timings_50ms{j}>=delay_timings(j,1)) & (lick_timings_50ms{j}<delay_timings(j,2))));
    end
    delay_licknum_50ms_1(1:idenrev_trial-1) = [];
    
    % 2. rev day
    load(beh_file_rev{ianimal})
    
    lick_diff = diff(lickTime, 1);
    short_lick_diff = find(lick_diff(:,1)<=50)+1;
    lickTime_50ms = lickTime;
    lickTime_50ms(short_lick_diff,:)=[];
    
    lick_timings_50ms = cell(nTrial,1);
    delay_timings = stateTime(:,3:4)-stateTime(:,1)+1;
    delay_licknum_50ms_2 = zeros(size(lickNum));
    for j = 1:nTrial
        jtrial_ind = find(lickTime_50ms(:,2)==j);
        lick_timings_50ms{j} = lickTime_50ms(jtrial_ind,1) - stateTime(j,1)+1;
        delay_licknum_50ms_2(j) = length(find((lick_timings_50ms{j}>=delay_timings(j,1)) & (lick_timings_50ms{j}<delay_timings(j,2))));
    end
    
    % 3. concatenate
    delay_licknum = [delay_licknum_50ms_1; delay_licknum_50ms_2];
    
    %% ANCCR    
    % fminsearch options
    maxit = 20000;
    maxeval = 20000;
    op=optimset('fminsearch');
    op.MaxIter=maxit;
    op.MaxFunEvals=maxeval;
    
    T = round(mean([IRI_rw1 IRI_pn1 IRI_rw2 IRI_pn2]));
    input_set = {eventlog,delay_licknum,T};
    iter_num = 100;
    error_set = zeros(1,iter_num);
    xpar_set = cell(1,iter_num);
    for i = 1:iter_num
        tic
        fprintf('%d/%d iteration (%d/%d animal)\n',i,iter_num,ianimal,length(beh_file_rev))
        ipar = rand(1,3);
        [xpar,error] = fminsearch(@fun_ANCCR_3var_Ras, ipar, op, input_set);
        if xpar(2)<10^(-3) % YJ: treat max variable size and reduce time
            xpar(2) = 10^(-3);
        end
        error_set(i) = error; % -r
        xpar_set{i} = xpar;
        toc
    end
    [~,min_id] = min(error_set);
    best_xpar = xpar_set{min_id};
    
    animal_name = strsplit(beh_file_rev{ianimal},'_');
    animal_name = animal_name{1};
    
    fig_lick_dir = ['E:\snl Dropbox\Jee Hyun\Ca imaging\Mat_code\ANCCR\fit VIP 3var Ras\' animal_name '\'];
    [~,msg1,~] = mkdir(fig_lick_dir);
    
    save([fig_lick_dir 'ANCCR_3var.mat'],'xpar_set','error_set','min_id','best_xpar')
end