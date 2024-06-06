clear; close all;

%% Optogenetics: trials to criterion
main_path = 'E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\ChRmine';
cd(main_path)

%%
file_list = ls('*.mat');
file_num = size(file_list,1);

if mod(file_num,2)>0; disp('File number error: not even'); return; end

mod_order = repmat([1 2],1,file_num/2)==1;

%%
trials_to_rev_ses = zeros(file_num,2);

for ises = 1:file_num
    load(file_list(ises,:))
    first_mod = mod_order(ises);
    ses_name = strsplit(file_list(ises,:),'_');

    %% reversal 
    [idenrev_trial,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
    idenrev_trial = unique(idenrev_trial)+1;

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

    %% calculation
    trial_cut = [idenrev_trial; nTrial];
    rev_point = zeros(2,1);

    for irev = 1:2
        %% rev point
        outcomeID =  outcomeIdentity(idenrev_trial(irev),:);
        pn_cue = find(outcomeID==3);
        rw_cue = find(outcomeID==2);

        for itrial = trial_cut(irev):trial_cut(irev+1)
            lickNum = delay_licknum_50ms(itrial-49:itrial); % behavior: itrial-49:itrial
            cueData = odorCue(itrial-49:itrial);

            [pANOVA,~,stats] = anova1(lickNum,cueData,'off');
            [c, ~, ~, gnames] = multcompare(stats,'display','off');
            out = (ismember(c(:,1),pn_cue) & ismember(c(:,2),pn_cue)) |...
                (ismember(c(:,1),rw_cue) & ismember(c(:,2),rw_cue)); % R-R, P-P out

            if pANOVA<0.05 && sum(c(~out,6)<0.05)==4 &&...
                    (nanmean(lickNum(cueData==rw_cue(1)-1)) > nanmean(lickNum(cueData==pn_cue(1)-1))) &&...
                    (nanmean(lickNum(cueData==rw_cue(1)-1)) > nanmean(lickNum(cueData==pn_cue(2)-1))) &&...
                    (nanmean(lickNum(cueData==rw_cue(2)-1)) > nanmean(lickNum(cueData==pn_cue(1)-1))) &&...
                    (nanmean(lickNum(cueData==rw_cue(2)-1)) > nanmean(lickNum(cueData==pn_cue(2)-1)))
                rev_point(irev) = itrial;
                break
            end
        end
    end

    %% trials to rev
    trials_to_rev = rev_point - idenrev_trial;

    %% output: OFF - ON order
    if first_mod
        trials_to_rev = flipud(trials_to_rev);
    end

    trials_to_rev_ses(ises,:) = trials_to_rev';
end

%% fig
f_ai = figure('PaperUnits','Centimeters','PaperPosition',[1 1 9 9]);
ax = axes;
ax.FontSize = 20;
hold on

iresult = zeros(size(trials_to_rev_ses,1)/2,2);
for ises = 1:size(trials_to_rev_ses,1)/2
    iresult(ises,:) = mean(trials_to_rev_ses(2*ises-1:2*ises,:));
end

bar(1,mean(iresult(:,1)),'FaceColor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'linewidth',1.5)
alpha(0.2)
errorbar(1,mean(iresult(:,1)),std(iresult(:,1))/sqrt(length(iresult(:,1))),'color',[0.5 0.5 0.5],'capsize',0,'linewidth',1.5)
bar(2,mean(iresult(:,2)),'FaceColor','r','edgecolor','r','linewidth',1.5)
alpha(0.2)
errorbar(2,mean(iresult(:,2)),std(iresult(:,2))/sqrt(length(iresult(:,2))),'r','capsize',0,'linewidth',1.5)

for ises = 1:size(iresult,1)
    plot([1 2],iresult(ises,:),'-','color',[0.5 0.5 0.5])
end

ylabel('Trials to criterion')
yticks(0:50:200)
ylim([0 200])
xticks([1 2])
xlim([0 3])
xticklabels({'OFF','ON'})

% [~,p_t,~,stats] = ttest(iresult(:,1),iresult(:,2));
% 
% fprintf('p = %.3f, t(%d) = %.3f\n',p_t,stats.df,stats.tstat);
