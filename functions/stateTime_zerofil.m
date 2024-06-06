function [stateTime, stateTime_original, state_duration, trial_duration, state_time] = stateTime_zerofil(stateTime)
[zerotime_trial,zerotime_state] = find(stateTime==0);
zerotime_trial1 = unique(zerotime_trial);
zerotime_state1 = unique(zerotime_state);
if isempty(zerotime_trial1)
    state_duration = mean(diff(stateTime,1,2))/1000;
    state_time = cumsum(state_duration);
    disp(state_duration)
    temp_trial_duration = sort(diff(stateTime(:,1)));
    trial_duration = ceil(temp_trial_duration(end-1)/1000);
    disp(trial_duration)
else %~ismember(1,zerotime_state1)
    stateTime_nozero = stateTime;
    stateTime_nozero(zerotime_trial1,:) = [];
    state_duration = mean(diff(stateTime_nozero,1,2))/1000;
    state_time = cumsum(state_duration);
    disp(state_duration)
    ref_state = find(~ismember([1 2 3 4 5],zerotime_state1),1);
    temp_trial_duration = sort(diff(stateTime(:,ref_state)));
    trial_duration = ceil(temp_trial_duration(end-1)/1000);
    disp(trial_duration)
end

stateTime_original = stateTime;
state_duration_round = round(state_duration*1000);
while ismember(0,stateTime)
    [temp_trial,temp_state] = find(stateTime==0,1,'last');
    if temp_state<5
        ref_time = stateTime(temp_trial,temp_state+1);
        if ref_time~=0
            stateTime(temp_trial,temp_state) = ref_time-state_duration_round(temp_state);
        end
    else
        if stateTime(temp_trial,4)~=0
            stateTime(temp_trial,5) = stateTime(temp_trial,4)+state_duration_round(4);
        elseif stateTime(temp_trial,3)~=0
            stateTime(temp_trial,4) = stateTime(temp_trial,3)+state_duration_round(3);
            stateTime(temp_trial,5) = stateTime(temp_trial,4)+state_duration_round(4);
        elseif stateTime(temp_trial,2)~=0
            stateTime(temp_trial,3) = stateTime(temp_trial,2)+state_duration_round(2);
            stateTime(temp_trial,4) = stateTime(temp_trial,3)+state_duration_round(3);
            stateTime(temp_trial,5) = stateTime(temp_trial,4)+state_duration_round(4);
        elseif stateTime(temp_trial,1)~=0
            stateTime(temp_trial,2) = stateTime(temp_trial,1)+state_duration_round(1);
            stateTime(temp_trial,3) = stateTime(temp_trial,2)+state_duration_round(2);
            stateTime(temp_trial,4) = stateTime(temp_trial,3)+state_duration_round(3);
            stateTime(temp_trial,5) = stateTime(temp_trial,4)+state_duration_round(4);
        else
            fprintf('Check %d trial: all zeros\n',temp_trial);
            return
        end
    end
end
end