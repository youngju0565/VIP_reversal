clear; close all;

%% Reversal: example session
load('data_VIP_CC\fig2B_example_DMSO_session.mat')

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
trial_range = idenrev_trial-90:nTrial-30;

outcomeID =  outcomeIdentity(1,:);
pn_cue = find(outcomeID==3);
rw_cue = find(outcomeID==2);

outcomeID_rev = outcomeIdentity(idenrev_trial,:);
pn_cue_rev = find(outcomeID_rev==3);
rw_cue_rev = find(outcomeID_rev==2);

mean_data = zeros(4,length(trial_range));
sem_data = zeros(4,length(trial_range));
pass_data = zeros(1,length(trial_range));
success_data = zeros(1,length(trial_range));
for itrial = 1:length(trial_range)
    temp_range = trial_range(itrial)-29:trial_range(itrial)+30;
    lickNum = delay_licknum_50ms(temp_range);
    cueData = odorCue(temp_range);
    
    [pANOVA,~,stats] = anova1(lickNum,cueData,'off');
    [c, ~, ~, gnames] = multcompare(stats,'display','off');
    out = (ismember(c(:,1),pn_cue) & ismember(c(:,2),pn_cue)) |...
        (ismember(c(:,1),rw_cue) & ismember(c(:,2),rw_cue));
    
    if pANOVA<0.05 && sum(c(~out,6)<0.05)==4 &&...
            (nanmean(lickNum(cueData==rw_cue(1)-1)) > nanmean(lickNum(cueData==pn_cue(1)-1))) &&...
            (nanmean(lickNum(cueData==rw_cue(1)-1)) > nanmean(lickNum(cueData==pn_cue(2)-1))) &&...
            (nanmean(lickNum(cueData==rw_cue(2)-1)) > nanmean(lickNum(cueData==pn_cue(1)-1))) &&...
            (nanmean(lickNum(cueData==rw_cue(2)-1)) > nanmean(lickNum(cueData==pn_cue(2)-1)))
        pass_data(itrial) = 1;
    end
    
    out_rev = (ismember(c(:,1),pn_cue_rev) & ismember(c(:,2),pn_cue_rev)) |...
        (ismember(c(:,1),rw_cue_rev) & ismember(c(:,2),rw_cue_rev));
    if pANOVA<0.05 && sum(c(~out_rev,6)<0.05)==4 &&...
            (nanmean(lickNum(cueData==rw_cue_rev(1)-1)) > nanmean(lickNum(cueData==pn_cue_rev(1)-1))) &&...
            (nanmean(lickNum(cueData==rw_cue_rev(1)-1)) > nanmean(lickNum(cueData==pn_cue_rev(2)-1))) &&...
            (nanmean(lickNum(cueData==rw_cue_rev(2)-1)) > nanmean(lickNum(cueData==pn_cue_rev(1)-1))) &&...
            (nanmean(lickNum(cueData==rw_cue_rev(2)-1)) > nanmean(lickNum(cueData==pn_cue_rev(2)-1)))
        success_data(itrial) = 1;
    end
    
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
pass_data(pass_data==0)=nan;
success_data(success_data==0)=nan;

%%
trial_range = trial_range-idenrev_trial;

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
pass_plot=plot(trial_range,pass_data*(axisVal(4)+1));
pass_plot.LineStyle = 'none'; pass_plot.Marker = 's'; pass_plot.MarkerEdgeColor = 'none';
pass_plot.MarkerFaceColor = 'k'; pass_plot.MarkerSize = 4;
success_plot=plot(trial_range,success_data*(axisVal(4)+1));
success_plot.LineStyle = 'none'; success_plot.Marker = 's'; success_plot.MarkerEdgeColor = 'none';
success_plot.MarkerFaceColor = 'k'; success_plot.MarkerSize = 4;
plot([0 0],[0 axisVal(4)+1],'--','color',[0.5 0.5 0.5]) % rev onset
ylabel('Anticipatory lick number')
xlabel('Trial since reversal')
