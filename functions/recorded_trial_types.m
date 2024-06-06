function [types, cues, types_rec, no_lick, no_lick1, outcomes] = recorded_trial_types(odorCue, waterReward, outcomeIdentity, stateTime, lickTime, trial_idx_rec, nTrial)
    type1 = []; type2 = []; type3 = []; type4 = []; type5 = []; type6 = []; type7 = []; type8 = [];
    no_lick = []; no_lick1 = []; no_lick_flag = false;
    type1_rec = []; type2_rec = []; type3_rec = []; type4_rec = []; type5_rec = []; type6_rec = []; type7_rec = []; type8_rec = [];
    k_rec = 0;
    rwded = []; rwded_rec = []; unrwded = []; unrwded_rec = [];
    puned = []; puned_rec = []; unpuned = []; unpuned_rec = [];
    
    for k = 1:nTrial
        cue_k = odorCue(k);
        rewarded_k = waterReward(k);
        iden_k = outcomeIdentity(k,cue_k+1);
        if no_lick_flag
            trial_onset_k = stateTime(k,1);
            if k<nTrial
                trial_offset_k = stateTime(k+1,1);
            else
                trial_offset_k = stateTime(k,5)+3000;
            end
            if isempty(find(((trial_onset_k<=lickTime(:,1))&(lickTime(:,1)<trial_offset_k)),1))
                no_lick1 = [no_lick1; k];
                continue
            else
                no_lick1 = [no_lick1; k];
                no_lick_flag = false;
                continue
            end
        end
        if rewarded_k && iden_k==2
            p0_onset_k = stateTime(k,4);
            if k<nTrial
                trial_offset_k = stateTime(k+1,1);
            else
                trial_offset_k = stateTime(k,5)+3000;
            end
            if isempty(find(((p0_onset_k<=lickTime(:,1))&(lickTime(:,1)<trial_offset_k)),1))
                no_lick = [no_lick; k];
                no_lick_flag = true;
                continue
            end
        end
    end
    
    no_lick(~ismember(no_lick,trial_idx_rec)) = [];
    no_lick1(~ismember(no_lick1,trial_idx_rec)) = [];
    
    if size(trial_idx_rec,1)~=1
        trial_idx_rec = trial_idx_rec';
    end
    
    for k = trial_idx_rec
        k_rec = k_rec+1;
        cue_k = odorCue(k);
        rewarded_k = waterReward(k);
        iden_k = outcomeIdentity(k,cue_k+1);
        if ismember(k,no_lick)||ismember(k,no_lick1)
            continue
        end
        switch cue_k
            case 0
                if rewarded_k; type1 = [type1; k]; type1_rec = [type1_rec; k_rec];
                else; type2 = [type2; k]; type2_rec = [type2_rec; k_rec]; end
            case 1
                if rewarded_k; type3 = [type3; k]; type3_rec = [type3_rec; k_rec];
                else; type4 = [type4; k]; type4_rec = [type4_rec; k_rec]; end
            case 2
                if rewarded_k; type5 = [type5; k]; type5_rec = [type5_rec; k_rec];
                else; type6 = [type6; k]; type6_rec = [type6_rec; k_rec]; end
            case 3
                if rewarded_k; type7 = [type7; k]; type7_rec = [type7_rec; k_rec];
                else; type8 = [type8; k]; type8_rec = [type8_rec; k_rec]; end
        end
        switch iden_k
            case 2 % reward
                if rewarded_k; rwded = [rwded; k]; rwded_rec = [rwded_rec; k_rec];
                else; unrwded = [unrwded; k]; unrwded_rec = [unrwded_rec; k_rec];
                end
            case 3 % punish
                if rewarded_k; puned = [puned; k]; puned_rec = [puned_rec; k_rec];
                else; unpuned = [unpuned; k]; unpuned_rec = [unpuned_rec; k_rec]; end
        end
    end
    types = {type1, type2, type3, type4, type5, type6, type7, type8};
    cues = {sort([type1;type2]), sort([type3;type4]), sort([type5;type6]), sort([type7;type8])};
    types_rec = {type1_rec, type2_rec, type3_rec, type4_rec, type5_rec, type6_rec, type7_rec, type8_rec};
    outcomes = {rwded, rwded_rec, unrwded, unrwded_rec, puned, puned_rec, unpuned, unpuned_rec};
    end