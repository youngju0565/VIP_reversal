clear; close all;

%% Classical conditioning: example session
load('data_VIP_CC\fig1C_example_DMSO_session.mat')
% load('data_VIP_CC\fig1C_example_CNO_session.mat')

%% reversal
[idenrev_trials,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
idenrev_trials = unique(idenrev_trials)+1;
idenrev_trial = idenrev_trials(1);

%% lick time
if length(lickTime)>1
    lick_diff = diff(lickTime, 1);
    short_lick_diff = find(lick_diff(:,1)<=50)+1;
    lickTime_50ms = lickTime;
    lickTime_50ms(short_lick_diff,:)=[];
elseif length(lickTime)==1
    lickTime_50ms = lickTime;
end

lick_timings_50ms = cell(nTrial,1);
delay_timings = stateTime(:,3:4)-stateTime(:,1)+1;
delay_licknum_50ms = zeros(size(lickNum));
for j = 1:nTrial
    jtrial_ind = find(lickTime_50ms(:,2)==j);
    lick_timings_50ms{j} = lickTime_50ms(jtrial_ind,1) - stateTime(j,1)+1;
    delay_licknum_50ms(j) = length(find((lick_timings_50ms{j}>=delay_timings(j,1)) & (lick_timings_50ms{j}<delay_timings(j,2))));
end

%%
trial_range = 30:idenrev_trial-1;

outcomeID =  outcomeIdentity(1,:);
pn_cue = find(outcomeID==3);
rw_cue = find(outcomeID==2);

mean_data = zeros(4,length(trial_range));
sem_data = zeros(4,length(trial_range));
for itrial = 1:length(trial_range)
    temp_range = trial_range(itrial)-29:trial_range(itrial)+30;
    lickNum = delay_licknum_50ms(temp_range);
    cueData = odorCue(temp_range);
    
    % A: rw-rw, B: rw-pn, C: pn-rw, D: pn-pn
    mean_data(1,itrial) = nanmean(lickNum(cueData==rw_cue(ismember(rw_cue,idenrev_cue))-1));
    mean_data(2,itrial) = nanmean(lickNum(cueData==rw_cue(~ismember(rw_cue,idenrev_cue))-1));
    mean_data(3,itrial) = nanmean(lickNum(cueData==pn_cue(ismember(pn_cue,idenrev_cue))-1));
    mean_data(4,itrial) = nanmean(lickNum(cueData==pn_cue(~ismember(pn_cue,idenrev_cue))-1));
    sem_data(1,itrial) = nanstd(lickNum(cueData==rw_cue(ismember(rw_cue,idenrev_cue))-1))/sqrt(sum(cueData==rw_cue(ismember(rw_cue,idenrev_cue))-1));
    sem_data(2,itrial) = nanstd(lickNum(cueData==rw_cue(~ismember(rw_cue,idenrev_cue))-1))/sqrt(sum(cueData==rw_cue(~ismember(rw_cue,idenrev_cue))-1));
    sem_data(3,itrial) = nanstd(lickNum(cueData==pn_cue(ismember(pn_cue,idenrev_cue))-1))/sqrt(sum(cueData==pn_cue(ismember(pn_cue,idenrev_cue))-1));
    sem_data(4,itrial) = nanstd(lickNum(cueData==pn_cue(~ismember(pn_cue,idenrev_cue))-1))/sqrt(sum(cueData==pn_cue(~ismember(pn_cue,idenrev_cue))-1));
end

%%
color_set = [57 89 168; 165 159 251; 228 175 175; 139 23 28]/255;
fill_x = [trial_range fliplr(trial_range)];

f1=figure('PaperUnits','Centimeters','PaperPosition',[1 1 4 3]);
ax = axes;
ax.FontSize = 8;
hold on
set(gca,'xcolor','k','ycolor','k')
for icue = 1:4
    fill_y = [mean_data(icue,:)+sem_data(icue,:) fliplr(mean_data(icue,:)-sem_data(icue,:))];
    fill_cue = fill(fill_x,fill_y,color_set(icue,:));
    fill_cue.EdgeColor = 'none'; alpha(fill_cue,0.5);
    plot(trial_range,mean_data(icue,:),'color',color_set(icue,:),'linewidth',1.5)
end
xlim([min(trial_range) max(trial_range)])
axisVal = axis;
ylim([0 axisVal(4)+1])
ylabel('Anticipatory lick number')
xlabel('Trial')
