function [eventlog,IRI_rw,IRI_pn] = behmat2eventlog(stateTime,outcomeIdentity,waterReward,nTrial,odorCue,outcomeProbability,lickTime)
%% transform CC behavior mat file to ANCCR input eventlog form
% calculate exact reward consume time
rw_period_set_time = zeros(nTrial,1);
waterReward_lick = waterReward;
for itrial = 1:nTrial
    outcomeID = outcomeIdentity(itrial,:);
    rw_period_set_time(itrial) = stateTime(itrial,4);
    if outcomeID(odorCue(itrial)+1)==2 && waterReward(itrial)
        temp_lick_onset = lickTime(find(lickTime>=stateTime(itrial,4),1));
        if itrial<nTrial && temp_lick_onset < stateTime(itrial+1,1)
            rw_period_set_time(itrial) = temp_lick_onset;
        end
        if itrial==nTrial && temp_lick_onset < stateTime(itrial,5)+3000
            rw_period_set_time(itrial) = temp_lick_onset;
        end
        if itrial<nTrial && temp_lick_onset >= stateTime(itrial+1,1)
            waterReward_lick(itrial) = 0;
        end
        if itrial==nTrial && temp_lick_onset >= stateTime(itrial,5)+3000
            waterReward_lick(itrial) = 0;
        end
    end
end

rw_period_set_time = (rw_period_set_time-stateTime(1,1)+5000)/1000;
cue_time = (stateTime(:,2)-stateTime(1,1)+5000)/1000;

% goal: CSrw or CSrw-pn 1, CSpn or CSpn-rw 2, CSnt 3
cue_set = unique(odorCue)+1;
cue1_set = find(outcomeProbability(1,:)>0&outcomeIdentity(1,:)==2);
cue1 = cue_set(ismember(cue_set,cue1_set));
cue2_set = find(outcomeProbability(1,:)>0&outcomeIdentity(1,:)==3);
cue2 = cue_set(ismember(cue_set,cue2_set));
cue3_set = find(outcomeProbability(1,:)==0);
cue3 = cue_set(ismember(cue_set,cue3_set));

odorCue_input = zeros(size(odorCue));
odorCue_temp = odorCue+1;
odorCue_input(odorCue_temp==cue1)=1;
odorCue_input(odorCue_temp==cue2)=2;
odorCue_input(odorCue_temp==cue3)=3;

% convert
R_pn=-1; % R_pn을 -1로 해두고, 할 때 마다 바꾸는 방식이 계산이 덜할듯
eventlog = zeros(nTrial*2,3);
for itrial = 1:nTrial
    eventlog(2*itrial-1,1) = odorCue_input(itrial);
    eventlog(2*itrial-1,2) = cue_time(itrial);
    eventlog(2*itrial-1,3) = 0;
    
    eventlog(2*itrial,2) = rw_period_set_time(itrial);
    outcomeID = outcomeIdentity(itrial,:);
    otucomeProb = outcomeProbability(itrial,:);
    if outcomeID(odorCue(itrial)+1)==2
        if waterReward_lick(itrial)
            eventlog(2*itrial,1) = 4;
            eventlog(2*itrial,3) = 1;
        elseif otucomeProb(odorCue(itrial)+1)>0
            eventlog(2*itrial,1) = 6;
            eventlog(2*itrial,3) = 0;
        elseif otucomeProb(odorCue(itrial)+1)==0
            eventlog(2*itrial,1) = 8;
            eventlog(2*itrial,3) = 0;
        end
    elseif outcomeID(odorCue(itrial)+1)==3
        if waterReward_lick(itrial)
            eventlog(2*itrial,1) = 5;
            eventlog(2*itrial,3) = R_pn;
        elseif otucomeProb(odorCue(itrial)+1)>0
            eventlog(2*itrial,1) = 7;
            eventlog(2*itrial,3) = 0;
        elseif otucomeProb(odorCue(itrial)+1)==0
            eventlog(2*itrial,1) = 8;
            eventlog(2*itrial,3) = 0;
        end
    end
end

% calculate mean inter-reward-interval (for T)
IRI_rw = round(mean(diff(eventlog(eventlog(:,1)==4,2))));
IRI_pn = round(mean(diff(eventlog(eventlog(:,1)==5,2))));
