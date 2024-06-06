clear; close all;

%% Trials to reversal success criterion
path_set = {'Gi','Gq'};

for ichemo = 1:2
    %%
    cd(['E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\' path_set{ichemo}])
    
    animal_name = ls;
    animal_name(1:2,:) = [];
    animal_num = size(animal_name,1);
    
    trials_to_pass = zeros(animal_num,3);
    
    main_path = ['E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\' path_set{ichemo}];
    
    trials_to_rev = zeros(animal_num,3);
    first_unexp_rw = nan(animal_num,3);
    
    %%
    for imouse = 1:animal_num
        cd([main_path '\' animal_name(imouse,:)])
        file_list = ls('*.mat');
        
        %%
        for ises = 1:3
            %%
            load(file_list(ises,:))
            ses_name = strsplit(file_list(ises,:),'_');
            
            %% reversal
            [idenrev_trials,idenrev_cue] = find(diff(outcomeIdentity(1:nTrial,:))); % RP
            idenrev_trials = unique(idenrev_trials)+1;
            idenrev_trial = idenrev_trials(1);
            
            %%
            first_outcome_trial = idenrev_trial+find(waterReward(idenrev_trial:end),1)-1;
            first_unexp_rw(imouse,ises) = outcomeIdentity(first_outcome_trial,odorCue(first_outcome_trial)+1);
            
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
            
            %% rev point
            outcomeID =  outcomeIdentity(idenrev_trial,:);
            pn_cue = find(outcomeID==3);
            rw_cue = find(outcomeID==2);
            rev_point = 0;
            
            for itrial = idenrev_trial:nTrial
                lickNum = delay_licknum_50ms(itrial-49:itrial);
                cueData = odorCue(itrial-49:itrial);
                
                [pANOVA,~,stats] = anova1(lickNum,cueData,'off');
                [c, ~, ~, gnames] = multcompare(stats,'display','off');
                out = (ismember(c(:,1),pn_cue) & ismember(c(:,2),pn_cue)) |...
                    (ismember(c(:,1),rw_cue) & ismember(c(:,2),rw_cue));
                
                if pANOVA<0.05 && sum(c(~out,6)<0.05)==4 &&...
                        (nanmean(lickNum(cueData==rw_cue(1)-1)) > nanmean(lickNum(cueData==pn_cue(1)-1))) &&...
                        (nanmean(lickNum(cueData==rw_cue(1)-1)) > nanmean(lickNum(cueData==pn_cue(2)-1))) &&...
                        (nanmean(lickNum(cueData==rw_cue(2)-1)) > nanmean(lickNum(cueData==pn_cue(1)-1))) &&...
                        (nanmean(lickNum(cueData==rw_cue(2)-1)) > nanmean(lickNum(cueData==pn_cue(2)-1)))
                    rev_point = itrial;
                    break
                end
            end
            if rev_point == 0
                rev_point = nTrial;
            end
            
            trials_to_rev(imouse,ises) = rev_point - idenrev_trial +1;
        end
    end
    
    %%
    ct = cbrewer('qual','Paired',12);
    bar_color_set{2}=ct(6,:); % cno
    
    ct = cbrewer('seq','Greys',9);
    ct = ct([1,4:9],:);
    bar_color_set{1}=ct(4,:); %dmso
    bar_color_set{3}=ct(4,:);
    
    f_ai = figure('PaperUnits','Centimeters','PaperPosition',[1 1 9 9]);
    ax = axes;
    ax.FontSize = 20;
    hold on
    set(gca,'xcolor','k','ycolor','k')
    
    iresult = trials_to_rev;
    
    for ibar = 1:3
        bar(ibar,mean(iresult(:,ibar)),'FaceColor',bar_color_set{ibar},'edgecolor',bar_color_set{ibar},'linewidth',1.5)
        alpha(0.2)
        errorbar(ibar,mean(iresult(:,ibar)),std(iresult(:,ibar))/sqrt(length(iresult(:,ibar))),'color',bar_color_set{ibar},'capsize',20,'linewidth',1.5)
    end
    
    for ises = 1:size(iresult,1)
        plot([1 2 3],iresult(ises,:),'-','color',[0.5 0.5 0.5])
    end
    
    ylabel('Trials to criterion')
    yticks(0:100:500)
    ylim([0 500])
    xticks([1 2 3])
    xlim([0 4])
    xticklabels({'DMSO','CNO','DMSO'})
    xtickangle(45)
    
end
