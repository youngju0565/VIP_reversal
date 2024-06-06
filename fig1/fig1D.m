clear; close all;

%% Classical conditioning, reward history effect: example session
Rewarded_lick_set = cell(1,2);
Omitted_lick_set = cell(1,2);

file_list = {'fig1D_exmaple_DMSO_session.mat','fig1D_example_CNO_session.mat'};

%%
for ises = 1:2
    load(['data_VIP_CC\' file_list{ises}])
    
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
    trial_range = idenrev_trial-100:idenrev_trial-1;
    
    outcomeID =  outcomeIdentity(1,:);
    rw_cue = find(outcomeID==2);
    
    CSrwd = ismember(odorCue+1,rw_cue);
    CSrwd_1 = [0; CSrwd(1:end-1)];
    Rwt_1 = [0; waterReward(1:end-1)&CSrwd(1:end-1)];
    
    Rewarded = CSrwd&CSrwd_1&Rwt_1;
    Omitted = CSrwd&CSrwd_1&~Rwt_1;
    
    trial_duration = 5;
    Rewarded_lick = zeros(sum(Rewarded(trial_range)),trial_duration*10^3/10^2);
    Omitted_lick = zeros(sum(Omitted(trial_range)),trial_duration*10^3/10^2);
    itrial_rw = 0;
    itrial_om = 0;
    for itrial = trial_range
        if Rewarded(itrial)
            itrial_rw = itrial_rw+1;
            trial_lick_times = lick_timings_50ms{itrial};
            [counts, edges] = histcounts(trial_lick_times, 0:10^2:trial_duration*10^3);
            Rewarded_lick(itrial_rw,:) = smoothdata(counts,'movmean',5)*10;
        end
        if Omitted(itrial)
            itrial_om = itrial_om+1;
            trial_lick_times = lick_timings_50ms{itrial};
            [counts, edges] = histcounts(trial_lick_times, 0:10^2:trial_duration*10^3);
            Omitted_lick(itrial_om,:) = smoothdata(counts,'movmean',5)*10;
        end
    end
    
    Rewarded_lick_set{ises} = Rewarded_lick;
    Omitted_lick_set{ises} = Omitted_lick;
    
end

%% fig
lick_data_set = {Rewarded_lick_set, Omitted_lick_set};

clear clr
ct = cbrewer('seq','Greys',9);
ct = ct([1,4:9],:);
clr{1}=ct(4,:); % dmso
ct = cbrewer('qual','Paired',12);
clr{2}=ct(6,:); % cno

linestyle_set = {'-','--'};

f1 = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3]);
ax = axes;
hold on
ax_fontsize = 10;
ax.FontSize = ax_fontsize;
label_fontsize = 10;
set(gca,'xcolor','k','ycolor','k')
fill_x = [edges(2:end) fliplr(edges(2:end))];
for idrug = 1:2
    for idata = 1:2
        temp_data = lick_data_set{idata}{idrug};
        fill_y = [mean(temp_data)+std(temp_data)/sqrt(size(temp_data,1)) ...
            fliplr(mean(temp_data)-std(temp_data)/sqrt(size(temp_data,1)))];
        fill_lick = fill(fill_x,fill_y,clr{idrug});
        fill_lick.EdgeColor = 'none';
        alpha(fill_lick, 0.2)
    end
end
for idrug = 1:2
    for idata = 1:2
        temp_data = lick_data_set{idata}{idrug};
        plot(edges(2:end), mean(temp_data), linestyle_set{idata},'color', clr{idrug}, 'linewidth',1)
    end
end
xticks(0.5*10^3 : 10^3 : trial_duration*10^3)
xticklabels(0:1:trial_duration)
xlim([edges(2) 4500])
axisValue = axis;
ylim([0 axisValue(4)])
line([0.5*10^3 0.5*10^3],[0 axisValue(4)],'color',[0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
line([1.5*10^3 1.5*10^3],[0 axisValue(4)],'color',[0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
line([2.5*10^3 2.5*10^3],[0 axisValue(4)],'color',[0.5 0.5 0.5],'linestyle','--','lineWidth', 1)
ylabel('Lick rate (Hz)','FontSize',label_fontsize)
xlabel('Time from cue onset (s)','FontSize',label_fontsize)
