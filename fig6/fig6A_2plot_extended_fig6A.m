clear; close all;

%%
n_num = num2str(12);
win_size = 1.5;
event_thr_set = [0.015 0.025 0.020];
hz_cut = 2;
if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
else; win_size_str = ''; end
if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
else; hz_cut_str = ''; end

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome'...
    win_size_str ' n' n_num hz_cut_str ' 230627\CPD_set.mat'])

load(cpd_filename_set{1})
x_var_num = size(CPD,2);
t_unify = cell(x_var_num,3);
p_unify = cell(x_var_num,3);
for iepoch = 1:3
    cell_count = 1;
    for ifile = 1:size(cpd_filename_set,2)
        if isempty(cpd_filename_set{ifile}); continue; end
        [pathstr,name,ext] = fileparts(cpd_filename_set{ifile});
        name(1:3) = 'FON';
        load(fullfile(pathstr,[name ext]))
        name = ['1234' name];
        name(1:7) = 't_value';
        load(fullfile(pathstr,[name ext]))
        for icell = 1:size(T_value,1)
            for ivar = 1:x_var_num
                p_unify{ivar,iepoch}(cell_count) = P_value(icell,ivar,iepoch);
                t_unify{ivar,iepoch}(cell_count) = T_value(icell,ivar,iepoch);
            end
            cell_count = cell_count+1;
        end
    end
end
data_draw_csrw = cell2mat(t_unify(1,:)')';
data_draw_cspn = cell2mat(t_unify(2,:)')';
p_csrw = cell2mat(p_unify(1,:)')';
p_cspn = cell2mat(p_unify(2,:)')';

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome term'...
    win_size_str ' n' n_num hz_cut_str ' 230627\CPD_set.mat'])
load(cpd_filename_set{1})
x_var_num = size(CPD_pun,2);
t_unify_rwd = cell(x_var_num,3); t_unify_pun = cell(x_var_num,3);
p_unify_rwd = cell(x_var_num,3); p_unify_pun = cell(x_var_num,3);
for iepoch = 1:3
    cell_count = 1;
    for ifile = 1:size(cpd_filename_set,2)
        if isempty(cpd_filename_set{ifile}); continue; end
        [pathstr,name,ext] = fileparts(cpd_filename_set{ifile});
        name(1:3) = 'FON';
        load(fullfile(pathstr,[name ext]))
        name = ['1234' name];
        name(1:7) = 't_value';
        load(fullfile(pathstr,[name ext]))
        for icell = 1:size(T_value_pun,1)
            for ivar = 1:x_var_num
                p_unify_rwd{ivar,iepoch}(cell_count) = P_value_rwd(icell,ivar,iepoch);
                t_unify_rwd{ivar,iepoch}(cell_count) = T_value_rwd(icell,ivar,iepoch);
                p_unify_pun{ivar,iepoch}(cell_count) = P_value_pun(icell,ivar,iepoch);
                t_unify_pun{ivar,iepoch}(cell_count) = T_value_pun(icell,ivar,iepoch);
            end
            cell_count = cell_count+1;
        end
    end
end
data_draw_rw = cell2mat(t_unify_rwd(3,:)')';
data_draw_pn = cell2mat(t_unify_pun(4,:)')';
p_rw = cell2mat(p_unify_rwd(3,:)')';
p_pn = cell2mat(p_unify_pun(4,:)')';

neuron_num = cell_count-1;

%%
axis_fontsize = 16;
label_fontsize = 20;

%% reward term: 3rd, punish term: 4th
title_set = {'Before','During','After'};
ylim_max = 5;

for iepoch = 1:3
    temp_csrw = data_draw_csrw(:,iepoch);
    temp_rw = data_draw_rw(:,iepoch);
    data_draw = [temp_csrw temp_rw];
    
    temp_p_csrw = p_csrw(:,iepoch);
    temp_p_rw = p_rw(:,iepoch);
    
    [f_scatter,a_scatter] = compare_scatter_cuedep_fitline_r_vip(data_draw,1,2,(temp_p_csrw<0.05)&(temp_p_rw<0.05));
    a_scatter.FontSize = axis_fontsize;
    title(title_set{iepoch},'fontsize',label_fontsize)
    xlim([-ylim_max ylim_max])
    ylim([-ylim_max ylim_max])
    xlabel('CS_r_w(t)','fontsize',label_fontsize)
    if iepoch == 1
        ylabel('Rw(t)','fontsize',label_fontsize)
    else
        ylabel('Rw(t)','fontsize',label_fontsize,'color','none')
    end

    %%
    temp_cspn = data_draw_cspn(:,iepoch);
    temp_pn = data_draw_pn(:,iepoch);
    data_draw = [temp_cspn temp_pn];
    
    temp_p_csrw = p_cspn(:,iepoch);
    temp_p_rw = p_pn(:,iepoch);
    
    [f_scatter2,a_scatter] = compare_scatter_cuedep_fitline_r_vip(data_draw,1,2,(temp_p_csrw<0.05)&(temp_p_rw<0.05));
    a_scatter.FontSize = axis_fontsize;
    title(title_set{iepoch},'fontsize',label_fontsize)
    xlim([-ylim_max ylim_max])
    ylim([-ylim_max ylim_max])
    xlabel('CS_p_n(t)','fontsize',label_fontsize)
    if iepoch == 1
        ylabel('Pn(t)','fontsize',label_fontsize)
    else
        ylabel('Pn(t)','fontsize',label_fontsize,'color','none')
    end
    
end
