clear; close all;

%%
n_num = num2str(9);
win_size = 1.5;
event_thr_set = [0.015 0.025 0.020];
hz_cut = 2;
if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
else; win_size_str = ''; end
if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
else; hz_cut_str = ''; end

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression norev outcome term 231108'...
    win_size_str ' n' n_num  hz_cut_str '\SRC_set.mat'])

load(src_filename_set{1})
x_var_num = size(SRC,2);
t_unify = cell(x_var_num,1);
p_unify = cell(x_var_num,1);
t_unify_rwd = cell(x_var_num,1); t_unify_pun = cell(x_var_num,1);
p_unify_rwd = cell(x_var_num,1); p_unify_pun = cell(x_var_num,1);
cell_count = 1;
for ifile = 1:size(src_filename_set,2)
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name(1:3) = 'FON';
    load(fullfile(pathstr,[name ext]))
    name = ['1234' name];
    name(1:7) = 't_value';
    load(fullfile(pathstr,[name ext]))
    for icell = 1:size(T_value,1)
        for ivar = 1:x_var_num
            p_unify{ivar}(cell_count) = P_value(icell,ivar);
            t_unify{ivar}(cell_count) = T_value(icell,ivar);
            p_unify_rwd{ivar}(cell_count) = P_value_rwd(icell,ivar);
            t_unify_rwd{ivar}(cell_count) = T_value_rwd(icell,ivar);
            p_unify_pun{ivar}(cell_count) = P_value_pun(icell,ivar);
            t_unify_pun{ivar}(cell_count) = T_value_pun(icell,ivar);
        end
        cell_count = cell_count+1;
    end
end

data_draw_csrw = cell2mat(t_unify(1,:)')';
data_draw_cspn = cell2mat(t_unify(2,:)')';
p_csrw = cell2mat(p_unify(1,:)')';
p_cspn = cell2mat(p_unify(2,:)')';

data_draw_rw = cell2mat(t_unify_rwd(3,:)')';
data_draw_pn = cell2mat(t_unify_pun(4,:)')';
p_rw = cell2mat(p_unify_rwd(3,:)')';
p_pn = cell2mat(p_unify_pun(4,:)')';

neuron_num = cell_count-1;

%%
axis_fontsize = 16;
label_fontsize = 20;

%% reward term: 3rd, punish term: 4th
data_cue = {data_draw_csrw data_draw_cspn};
data_out = {data_draw_rw data_draw_pn};
p_cue = {p_csrw p_cspn};
p_out = {p_rw p_pn};
title_set = {'Reward','Punishment'};
xlabel_set = {'CS_r_w(t)','CS_p_n(t)'};
ylabel_set = {'Rw(t)','Pn(t)'};
ylim_max = 5;

for iout = 1:2
    temp_cs = data_cue{iout};
    temp_out = data_out{iout};
    data_draw = [temp_cs temp_out];
    
    temp_p_cue = p_cue{iout};
    temp_p_out = p_out{iout};
    
    [f_scatter,a_scatter] = compare_scatter_cuedep_fitline_r_vip(data_draw,1,2,(temp_p_cue<0.05)&(temp_p_out<0.05));
    a_scatter.FontSize = axis_fontsize;
    set(gca,'xcolor','k','ycolor','k')
    [~,p] = ttest(data_draw(:,1),data_draw(:,2));
    title('')
    xlim([-ylim_max ylim_max])
    ylim([-ylim_max ylim_max])
    xlabel(xlabel_set{iout},'fontsize',label_fontsize)
    ylabel(ylabel_set{iout},'fontsize',label_fontsize)
end
