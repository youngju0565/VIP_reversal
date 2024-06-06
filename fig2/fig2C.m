clear; close all;

%% Anticipatory lick rate
path_set = {'Gi','Gq'};

for ichemo = 1:2
    %%
    cd(['E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\' path_set{ichemo}])
    
    animal_name = ls;
    animal_name(1:2,:) = [];
    animal_num = size(animal_name,1);
    
    trials_to_pass = zeros(animal_num,3);
    
    main_path = ['E:\snl Dropbox\Jee Hyun\manuscript_CC\Code_and_Data\data_VIP_CC\' path_set{ichemo}];
    
    rev_point_set = zeros(3,animal_num);
    right_num_set = zeros(3,animal_num);
    rp_dmso1_set = cell(1,animal_num);
    rp_cno_set = cell(1,animal_num);
    rp_dmso2_set = cell(1,animal_num);
    pr_dmso1_set = cell(1,animal_num);
    pr_cno_set = cell(1,animal_num);
    pr_dmso2_set = cell(1,animal_num);
    
    %%
    for imouse = 1:animal_num
        cd([main_path '\' animal_name(imouse,:)])
        file_list = ls('*.mat');
        
        trial_win = 60;
        step_size = 15;
        
        rp_set = cell(1,3);
        pr_set = cell(1,3);
        rev_point = zeros(1,3);
        right_num = zeros(1,3);
        
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
            start_trial = idenrev_trial-100;
            trial_range = start_trial:idenrev_trial+200;
            
            outcomeID =  outcomeIdentity(idenrev_trial,:);
            pn_cue = find(outcomeID==3);
            rw_cue = find(outcomeID==2);
            
            [~,rev_cue] = find(diff(outcomeIdentity(idenrev_trial-1:idenrev_trial,:))); % RP
            rp_cue = pn_cue(ismember(pn_cue,rev_cue));
            pr_cue = rw_cue(ismember(rw_cue,rev_cue));
            
            temp_win = (1:trial_win)-step_size;
            istep = 0;
            temp_step_num = floor((length(trial_range)-trial_win)/step_size)+1;
            rp_lick = zeros(1,temp_step_num);
            pr_lick = zeros(1,temp_step_num);
            temp_rev_point = nan;
            for istep = 1:temp_step_num
                temp_win = temp_win+step_size;
                if isnan(temp_rev_point) && (min(trial_range(temp_win))+trial_win/2)>=idenrev_trial
                    temp_rev_point = istep;
                end
                temp_lick = delay_licknum_50ms(trial_range(temp_win));
                rp_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(rp_cue-1)));
                pr_lick(istep) = mean(temp_lick(odorCue(trial_range(temp_win))==(pr_cue-1)));
            end
            
            rp_set{ises} = rp_lick;
            pr_set{ises} = pr_lick;
            rev_point(ises) = temp_rev_point;
            right_num(ises) = temp_step_num-temp_rev_point;
        end
        
        %%
        rp_dmso1_set{imouse} = rp_set{1};
        rp_cno_set{imouse} = rp_set{2};
        rp_dmso2_set{imouse} = rp_set{3};
        pr_dmso1_set{imouse} = pr_set{1};
        pr_cno_set{imouse} = pr_set{2};
        pr_dmso2_set{imouse} = pr_set{3};
        rev_point_set(:,imouse) = rev_point';
        right_num_set(:,imouse) = right_num';
    end
    
    %% cut
    x_left = min(min(rev_point_set));
    x_right = min(min(right_num_set));
    rp_dmso1_set_cut = cell(1,animal_num);
    rp_cno_set_cut = cell(1,animal_num);
    rp_dmso2_set_cut = cell(1,animal_num);
    pr_dmso1_set_cut = cell(1,animal_num);
    pr_cno_set_cut = cell(1,animal_num);
    pr_dmso2_set_cut = cell(1,animal_num);
    for ises = 1:animal_num
        rp_dmso1_set_cut{ises} = rp_dmso1_set{ises}((rev_point_set(1,ises)-x_left)+1:rev_point_set(1,ises)+x_right);
        rp_cno_set_cut{ises} = rp_cno_set{ises}((rev_point_set(2,ises)-x_left)+1:rev_point_set(2,ises)+x_right);
        rp_dmso2_set_cut{ises} = rp_dmso2_set{ises}((rev_point_set(3,ises)-x_left)+1:rev_point_set(3,ises)+x_right);
        pr_dmso1_set_cut{ises} = pr_dmso1_set{ises}((rev_point_set(1,ises)-x_left)+1:rev_point_set(1,ises)+x_right);
        pr_cno_set_cut{ises} = pr_cno_set{ises}((rev_point_set(2,ises)-x_left)+1:rev_point_set(2,ises)+x_right);
        pr_dmso2_set_cut{ises} = pr_dmso2_set{ises}((rev_point_set(3,ises)-x_left)+1:rev_point_set(3,ises)+x_right);
    end
    
    %% plot
    x_plot = ((1-x_left:x_right))*step_size;
    
    linewidth_set = [1.5 3];
    fontsize_ax = 20;
    fontsize_label = 24;
    fig_size = [1 1 12.65 6.85];
    r_color = [228 175 175; 139 23 28]/255;
    b_color = [165 159 251; 57 89 168]/255;
    ylim_max = 8;
    sigmaker_size = 10;
    
    f1 = figure('PaperUnits','Centimeters','PaperPosition',fig_size);
    ax1 = axes;
    ax1.FontSize = fontsize_ax;
    hold on
    title('DMSO','fontsize',fontsize_label)
    data_1 = rp_dmso1_set_cut;
    data_2 = pr_dmso1_set_cut;
    for ises = 1:animal_num
        plot(x_plot, data_1{ises}, 'color', b_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
        plot(x_plot, data_2{ises}, 'color', r_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    end
    ylim([0 ylim_max])
    yticks(0:4:ylim_max)
    xlim([min(x_plot)-step_size max(x_plot)+step_size])
    plot([0 0],[0 100],'color',[0.5 0.5 0.5],'linestyle',':')
    errorbar(x_plot, mean(cell2mat(data_1')), std(cell2mat(data_1'))/sqrt(animal_num),...
        'linewidth',linewidth_set(2),'color',b_color(2,:),'capsize',0)
    errorbar(x_plot, mean(cell2mat(data_2')), std(cell2mat(data_2'))/sqrt(animal_num),...
        'linewidth',linewidth_set(2),'color',r_color(2,:),'capsize',0)
    ylabel('Anticipatory lick','fontsize',fontsize_label,'color','none')
    xlabel('Trial since reversal','fontsize',fontsize_label,'color','none')
    ttest_result = zeros(size(x_plot));
    temp1 = cell2mat(data_1');
    temp2 = cell2mat(data_2');
    for istep = 1:length(x_plot)
        ttest_result(istep) = ttest(temp1(:,istep),temp2(:,istep));
    end
    ttest_result(ttest_result==0)=nan;
    plot(x_plot,ttest_result*ylim_max*0.95,'linestyle','none','marker','s','markeredgecolor','k','markerfacecolor','k',...
        'markersize',sigmaker_size)
    plot([x_plot(find(ttest_result==1&mean(cell2mat(data_2'))>mean(cell2mat(data_1')),1)) x_plot(find(ttest_result==1&mean(cell2mat(data_2'))>mean(cell2mat(data_1')),1))],...
        [0 ylim_max*0.95],'--k','linewidth',2)
    
    f2 = figure('PaperUnits','Centimeters','PaperPosition',fig_size);
    ax2 = axes;
    ax2.FontSize = fontsize_ax;
    hold on
    title('CNO','fontsize',fontsize_label)
    data_1 = rp_cno_set_cut;
    data_2 = pr_cno_set_cut;
    for ises = 1:animal_num
        plot(x_plot, data_1{ises}, 'color', b_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
        plot(x_plot, data_2{ises}, 'color', r_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    end
    ylim([0 ylim_max])
    yticks(0:4:ylim_max)
    xlim([min(x_plot)-step_size max(x_plot)+step_size])
    plot([0 0],[0 100],'color',[0.5 0.5 0.5],'linestyle',':')
    errorbar(x_plot, mean(cell2mat(data_1')), std(cell2mat(data_1'))/sqrt(animal_num),...
        'linewidth',linewidth_set(2),'color',b_color(2,:),'capsize',0)
    errorbar(x_plot, mean(cell2mat(data_2')), std(cell2mat(data_2'))/sqrt(animal_num),...
        'linewidth',linewidth_set(2),'color',r_color(2,:),'capsize',0)
    ylabel('Anticipatory lick','fontsize',fontsize_label)
    xlabel('Trial since reversal','fontsize',fontsize_label,'color','none')
    ttest_result = zeros(size(x_plot));
    temp1 = cell2mat(data_1');
    temp2 = cell2mat(data_2');
    for istep = 1:length(x_plot)
        ttest_result(istep) = ttest(temp1(:,istep),temp2(:,istep));
    end
    ttest_result(ttest_result==0)=nan;
    plot(x_plot,ttest_result*ylim_max*0.95,'linestyle','none','marker','s','markeredgecolor','k','markerfacecolor','k',...
        'markersize',sigmaker_size)
    try
        plot([x_plot(find(ttest_result==1&mean(cell2mat(data_2'))>mean(cell2mat(data_1')),1)) x_plot(find(ttest_result==1&mean(cell2mat(data_2'))>mean(cell2mat(data_1')),1))],...
            [0 ylim_max*0.95],'--k','linewidth',2)
    end
    
    f3 = figure('PaperUnits','Centimeters','PaperPosition',fig_size);
    ax3 = axes;
    ax3.FontSize = fontsize_ax;
    hold on
    title('DMSO','fontsize',fontsize_label)
    data_1 = rp_dmso2_set_cut;
    data_2 = pr_dmso2_set_cut;
    for ises = 1:animal_num
        plot(x_plot, data_1{ises}, 'color', b_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
        plot(x_plot, data_2{ises}, 'color', r_color(1,:), 'linewidth',linewidth_set(1),'marker','none')
    end
    ylim([0 ylim_max])
    yticks(0:4:ylim_max)
    xlim([min(x_plot)-step_size max(x_plot)+step_size])
    plot([0 0],[0 100],'color',[0.5 0.5 0.5],'linestyle',':')
    errorbar(x_plot, mean(cell2mat(data_1')), std(cell2mat(data_1'))/sqrt(animal_num),...
        'linewidth',linewidth_set(2),'color',b_color(2,:),'capsize',0)
    errorbar(x_plot, mean(cell2mat(data_2')), std(cell2mat(data_2'))/sqrt(animal_num),...
        'linewidth',linewidth_set(2),'color',r_color(2,:),'capsize',0)
    ylabel('Anticipatory lick','fontsize',fontsize_label,'color','none')
    xlabel('Trial since reversal','fontsize',fontsize_label)
    ttest_result = zeros(size(x_plot));
    temp1 = cell2mat(data_1');
    temp2 = cell2mat(data_2');
    for istep = 1:length(x_plot)
        ttest_result(istep) = ttest(temp1(:,istep),temp2(:,istep));
    end
    ttest_result(ttest_result==0)=nan;
    plot(x_plot,ttest_result*ylim_max*0.95,'linestyle','none','marker','s','markeredgecolor','k','markerfacecolor','k',...
        'markersize',sigmaker_size)
    plot([x_plot(find(ttest_result==1&mean(cell2mat(data_2'))>mean(cell2mat(data_1')),1)) x_plot(find(ttest_result==1&mean(cell2mat(data_2'))>mean(cell2mat(data_1')),1))],...
        [0 ylim_max*0.95],'--k','linewidth',2)
    
end