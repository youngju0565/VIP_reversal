clear; close all;

%% Reawrd history effect
path_set = {'Gi','Gq'};

for ichemo = 1:2
    %%
    cd(['E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\' path_set{ichemo}])
    
    animal_name = ls;
    animal_name(1:2,:) = [];
    animal_num = size(animal_name,1);
    
    trials_to_pass = zeros(animal_num,3);
    
    main_path = ['E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\' path_set{ichemo}];
    
    rwt1_lick_set = cell(1,animal_num);
    
    %%
    for imouse = 1:animal_num
        cd([main_path '\' animal_name(imouse,:)])
        file_list = ls('*.mat');
        
        csrw_rwt1_lick = zeros(3,2);
        
        %%
        for ises = 1:3
            %%
            load(file_list(ises,:))
            ses_name = strsplit(file_list(ises,:),'_');
            
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
            
            Reward = double((odorCue==rw_cue(1)-1|odorCue==rw_cue(2)-1)&waterReward);
            
            Reward_t1 = [0; Reward(1:end-1)];
            CS_t = double((odorCue==rw_cue(1)-1|odorCue==rw_cue(2)-1))*2-1; % 1 CSrwd -1 CSpun
            CS_t1 = [0; CS_t(1:end-1)];
            
            delay_lick_range = delay_licknum_50ms(trial_range);
            
            csrw_rwt1_lick(ises,1) = mean(delay_lick_range(CS_t1(trial_range)==1&CS_t(trial_range)==1&Reward_t1(trial_range)==1));
            csrw_rwt1_lick(ises,2) = mean(delay_lick_range(CS_t1(trial_range)==1&CS_t(trial_range)==1&Reward_t1(trial_range)==0));
            
        end
        %%
        rwt1_lick_set{1,imouse} = csrw_rwt1_lick;
    end
    
    %%
    if ichemo==2
        rwt1_lick_set(:,2)=[];
    end
    
    %%
    clear clr
    ct = cbrewer('qual','Paired',12);
    clr{3}=ct(6,:); % cno
    clr{4}=ct(6,:); % cno
    
    ct = cbrewer('seq','Greys',9);
    ct = ct([1,4:9],:);
    clr{1}=ct(4,:); %dmso
    clr{2}=ct(4,:); %dmso
    
    %%
    temp_data1 = cell2mat(cellfun(@(x) x(1,:),rwt1_lick_set(1,:),'uniformoutput',false)');
    temp_data2 = cell2mat(cellfun(@(x) x(3,:),rwt1_lick_set(1,:),'uniformoutput',false)');
    dmso = (temp_data1 + temp_data2)/2;
    
    cno = cell2mat(cellfun(@(x) x(2,:),rwt1_lick_set(1,:),'uniformoutput',false)');
    
    data_set = {dmso,cno};
    
    %%
    bar_color_set = {'b','w'};
    fig_bar = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 4]);
    ax = axes;
    ax.FontSize = 16;
    hold on
    set(gca,'xcolor','k','ycolor','k')
    xlabel_set = {'DMSO','CNO'};
    bar_width = 0.2;
    xtick_set = 0.7;
    for ibar = 1:2
        temp_data = data_set{ibar};
        bar(ibar*xtick_set-0.15,nanmean(temp_data(:,1)),bar_width,'edgecolor',clr{ibar*2-1},'facecolor',clr{ibar*2-1},'facealpha',0.5)
        bar(ibar*xtick_set+0.15,nanmean(temp_data(:,2)),bar_width,'edgecolor',clr{ibar*2},'facecolor','none','facealpha',0.5)
        errorbar(ibar*xtick_set-0.15,nanmean(temp_data(:,1)),nanstd(temp_data(:,1))/...
            sqrt(size(temp_data(:,1),1)),'Color',clr{ibar*2-1},'LineWidth',0.7);
        errorbar(ibar*xtick_set+0.15,nanmean(temp_data(:,2)),nanstd(temp_data(:,2))/...
            sqrt(size(temp_data(:,2),1)),'Color',clr{ibar*2},'LineWidth',0.7);
        for i = 1:size(temp_data(:,1),1)
            line([ibar*xtick_set-0.15 ibar*xtick_set+0.15], ...
                [temp_data(i,1) temp_data(i,2)], 'Color',ct(5,:),'LineWidth',0.3);
        end
    end
    set(ax,'Box','off','TickDir','out','FontSize',8,'XTick',[xtick_set xtick_set*2],'YTick',[0 5 10],...
        'XTickLabel',xlabel_set,'YLim',[0 10],'XLim',[xtick_set-0.15-0.2 2*xtick_set+0.15+0.2]);
    ylabel('Anticipatory lick num','FontSize',8);
    title(path_set{ichemo},'FontSize',8);
    xtickangle(45);
    
end
