clear; close all;

%%
rng(2);

load('data_VIP_CC\rev_dataname.mat')
beh_file_rev = beh_file;

load('data_VIP_CC\rev_d1_dataname.mat')

r_set_set = cell(1,length(beh_file_rev));

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
    
    %% load fit results
    animal_name = strsplit(beh_file_rev{ianimal},'_');
    animal_name = animal_name{1};
    
    fig_dir = ['E:\snl Dropbox\Jee Hyun\Ca imaging\Mat_code\ANCCR\fit VIP 3var Ras\' animal_name '\'];
    
    load([fig_dir 'ANCCR_3var.mat'])

    %% figs
    % other pars
    T = round(mean([IRI_rw1 IRI_pn1 IRI_rw2 IRI_pn2]));
    R_pn = best_xpar(1);
    eventlog(eventlog(:,3)<0,3) = R_pn;
    
    samplingperiod = best_xpar(3);
    alpha_anccr = best_xpar(2);
    alpha_r = 0.2;
    w = 0.5;
    k = 1;
    minimumrate = 10^(-3);
    maximumjitter = 0.1;
    beta = [0,0,0,1,1,0,0,0];
    threshold = 0.6;
    Tratio = 1.2;
    omidx = [6 4; 7 5];
    
    [DA,ANCCR,PRC,SRC,NC,Rs,~,Mij] = calculateANCCR(eventlog, T*Tratio, alpha_anccr, k,...
        samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,omidx);
    
    DA_data = cell(1,5);
    DA_data{1} = DA(eventlog(:,1)==1); % CSrw-pn
    DA_data{2} = DA(eventlog(:,1)==2); % CSpn-rw
    DA_data{3} = DA(eventlog(:,1)==3); % CSnt
    DA_data{4} = DA(eventlog(:,1)==4); % Rw
    DA_data{5} = DA(eventlog(:,1)==5); % Pn
    
    ANCCR_data = cell(1,8);
    ANCCR_data{1} = ANCCR(1,4,eventlog(:,1)>3); % CSrw-pn & Rw
    ANCCR_data{2} = ANCCR(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    ANCCR_data{3} = ANCCR(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    ANCCR_data{4} = ANCCR(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    ANCCR_data{5} = ANCCR(1,4,eventlog(:,1)<4); % CSrw-pn & Rw
    ANCCR_data{6} = ANCCR(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    ANCCR_data{7} = ANCCR(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    ANCCR_data{8} = ANCCR(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    ANCCR_data = cellfun(@squeeze,ANCCR_data,'uniformoutput',false);
    
    SRC_data = cell(1,8);
    SRC_data{1} = SRC(1,4,eventlog(:,1)>3); % CSrw-pn & Rw
    SRC_data{2} = SRC(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    SRC_data{3} = SRC(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    SRC_data{4} = SRC(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    SRC_data{5} = SRC(1,4,eventlog(:,1)<4); % CSrw-pn & Rw
    SRC_data{6} = SRC(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    SRC_data{7} = SRC(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    SRC_data{8} = SRC(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    SRC_data = cellfun(@squeeze,SRC_data,'uniformoutput',false);
    
    PRC_data = cell(1,8);
    PRC_data{1} = PRC(1,4,eventlog(:,1)>3); % CSrw-pn & Rw
    PRC_data{2} = PRC(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    PRC_data{3} = PRC(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    PRC_data{4} = PRC(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    PRC_data{5} = PRC(1,4,eventlog(:,1)<4); % CSrw-pn & Rw
    PRC_data{6} = PRC(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    PRC_data{7} = PRC(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    PRC_data{8} = PRC(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    PRC_data = cellfun(@squeeze,PRC_data,'uniformoutput',false);
    
    NC_data = cell(1,8);
    NC_data{1} = NC(1,4,eventlog(:,1)>3); % CSrw-pn & Rw
    NC_data{2} = NC(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    NC_data{3} = NC(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    NC_data{4} = NC(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    NC_data{5} = NC(1,4,eventlog(:,1)<4); % CSrw-pn & Rw
    NC_data{6} = NC(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    NC_data{7} = NC(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    NC_data{8} = NC(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    NC_data = cellfun(@squeeze,NC_data,'uniformoutput',false);
    
    PR_data = cell(1,8);
    PR_data{1} = Mij(1,4,eventlog(:,1)>3); % CSrw-pn & Rw
    PR_data{2} = Mij(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    PR_data{3} = Mij(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    PR_data{4} = Mij(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    PR_data{5} = Mij(1,4,eventlog(:,1)<4); % CSrw-pn & Rw
    PR_data{6} = Mij(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    PR_data{7} = Mij(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    PR_data{8} = Mij(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    PR_data = cellfun(@squeeze,PR_data,'uniformoutput',false);
    
    Rs_data = cell(1,8);
    Rs_data{1} = Rs(1,4,eventlog(:,1)>3); % CSrw-pn & Rw
    Rs_data{2} = Rs(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    Rs_data{3} = Rs(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    Rs_data{4} = Rs(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    Rs_data{5} = Rs(1,4,eventlog(:,1)<4); % CSrw-pn & Rw
    Rs_data{6} = Rs(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    Rs_data{7} = Rs(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    Rs_data{8} = Rs(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    Rs_data = cellfun(@squeeze,Rs_data,'uniformoutput',false);

    SRC_R_data = cell(1,8);
    SRC_R = SRC.*Rs;
    SRC_R_data{1} = SRC_R(1,4,eventlog(:,1)>3); % CSrw-pn & Rw % outcome
    SRC_R_data{2} = SRC_R(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
    SRC_R_data{3} = SRC_R(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
    SRC_R_data{4} = SRC_R(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
    SRC_R_data{5} = SRC_R(1,4,eventlog(:,1)<4); % CSrw-pn & Rw % cue
    SRC_R_data{6} = SRC_R(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
    SRC_R_data{7} = SRC_R(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
    SRC_R_data{8} = SRC_R(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
    SRC_R_data = cellfun(@squeeze,SRC_R_data,'uniformoutput',false);
    
    %% save variables for regression
    revses_idx = trial_idx>0;
    save([fig_dir 'variables_3var.mat'],'ANCCR_data','SRC_data','PRC_data','NC_data','PR_data','SRC_R_data','Rs_data','revses_idx')

end
