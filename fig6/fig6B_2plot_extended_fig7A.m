clear; close all;

%%
n_num = num2str(12);
win_size = 1.5;
step_size = 10; % 10 1
event_thr_set = [0.015 0.025 0.020];
hz_cut = 2;
win_size_sliding = 100;
if win_size~=1.5; win_size_str = [' ',num2str(win_size*1000),'ms'];
else; win_size_str = ''; end
if hz_cut~=2; hz_cut_str = [' ', num2str(event_thr_set(hz_cut))];
else; hz_cut_str = ''; end
if step_size~=10; step_str = [' ',num2str(step_size),'step'];
else; step_str = ''; end
if win_size_sliding~=100; trial_str = [' ',num2str(win_size_sliding),'trials'];
else; trial_str = ''; end

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome sliding revonset' ...
    win_size_str ' n' n_num hz_cut_str trial_str step_str ' 230627\SRC_set.mat'])
load(src_filename_set{1})
x_var_num = size(SRC,2);
step_num=  size(SRC,3);
T_value_unify = cell(x_var_num,step_num);
cell_count = 1;
for ifile = 1:size(src_filename_set,2)
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name = ['1234' name];
    name(1:7) = 't_value';
    load(fullfile(pathstr,[name ext]))
    for icell = 1:size(T_value,1)
        for istep = 1:step_num
            for ivar = 1:x_var_num
                T_value_unify{ivar,istep}(cell_count) = T_value(icell,ivar,istep);
            end
        end
        cell_count = cell_count+1;
    end
end
data_draw_csrw = cell2mat(T_value_unify(1,:)')';
data_draw_cspn = cell2mat(T_value_unify(2,:)')';

load(['E:\snl Dropbox\Jee Hyun\Ca imaging\Recording_data\figs\3_cue_RP\reversal_set regression outcome term sliding revonset' ...
    win_size_str ' n' n_num hz_cut_str trial_str step_str ' 230627\SRC_set.mat'])
load(src_filename_set{1})
x_var_num = size(SRC_pun,2);
step_num=  size(SRC_pun,3);
T_value_unify_rw = cell(x_var_num,step_num);
cell_count = 1;
for ifile = 1:size(src_filename_set,2)
    if isempty(src_filename_set{ifile}); continue; end
    [pathstr,name,ext] = fileparts(src_filename_set{ifile});
    name = ['1234' name];
    name(1:7) = 't_value';
    load(fullfile(pathstr,[name ext]))
    for icell = 1:size(T_value_rwd,1)
        for istep = 1:step_num
            for ivar = 1:x_var_num
                T_value_unify_rw{ivar,istep}(cell_count) = T_value_rwd(icell,ivar,istep);
            end 
        end
        cell_count = cell_count+1;
    end
end
data_draw_rw = cell2mat(T_value_unify_rw(3,:)')';

neuron_num = cell_count-1;

%%
r_set = zeros(1,step_num);
p_set = zeros(1,step_num);

for istep = 1:step_num
    temp_csrw = data_draw_csrw(:,istep);
    temp_rw = data_draw_rw(:,istep);
    [r,p] = corrcoef(temp_csrw,temp_rw);
    r_set(istep) = r(1,2);
    p_set(istep) = p(1,2);
end

%%
rev_step = min_step_num_left+1;
x_plot = ((win_size_sliding/2):step_size:(step_num*step_size+(win_size_sliding/2)-step_size));
x_plot = x_plot-x_plot(rev_step);

rev_trial_anova_win = [232-142 264-160 242-141 214-141 273-153 260-144 197-141 226-141 254-141 236-143 270-140 186-140];
success_trial = mean(rev_trial_anova_win);

%%
f_r = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 9]);
ax = axes;
hold on
ax.FontSize = 16;
plot(x_plot, r_set,'color','k','marker','o','markersize',8,'linewidth',1.5)
plot(x_plot(p_set<0.05),r_set(p_set<0.05),'color','k','marker','o','markersize',8,'markerfacecolor','k','linestyle','none')
plot([min(x_plot) max(x_plot)],[0 0],'--','color',[0.5 0.5 0.5])
plot([x_plot(rev_step) x_plot(rev_step)],[-1 1],'--','color',[0.5 0.5 0.5])
plot([success_trial success_trial],[-1 1],'--','color',[0.5 0.5 0.5])
xlim([min(x_plot) max([success_trial max(x_plot)])])
ylim([-0.6 0.3])
yticks(-0.6:0.3:0.3)
xlabel('Trial since reversal','fontsize',20)
ylabel('\itr','fontsize',20)
title('CS_r_w(t) b/w Rw(t) {\itt}-value','fontsize',20);
set(gca,'xcolor','k','ycolor','k')

%% Extended fig 7A
f_r_fit = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9.1 9]);
ax = axes;
set(gca,'xcolor','k','ycolor','k')
hold on
ax.FontSize = 16;

fit_input_idx_after = x_plot>0;
fit_input_idx = fit_input_idx_after;
function_form_rise=fittype('Rmin+R/(1+exp(-k*(x-a)))');
[cfit_rise,gof_rise] = fit(x_plot(fit_input_idx)',r_set(fit_input_idx)',function_form_rise,...
    'Lower',[0, -1, 40, 0], 'Upper',[2, 1, Inf, Inf], 'startpoint', [0.5 -0.5 mean(rev_trial_anova_win) 0],'MaxFunEvals',3000,'MaxIter',2000);

s1 = scatter(x_plot(fit_input_idx), r_set(fit_input_idx),10,'markeredgecolor','none','markerfacecolor','none');
p_fit = plot(cfit_rise);
legend off
namearray = {'Color','linewidth'};
p_fit_valuearray = {'r',2}; set(p_fit,namearray,p_fit_valuearray);

plot(x_plot, r_set,'color','k','marker','o','markersize',8,'linewidth',1.5)
plot(x_plot(p_set<0.05),r_set(p_set<0.05),'color','k','marker','o','markerfacecolor','k','linestyle','none','markersize',8)
plot([min(x_plot) max(x_plot)],[0 0],'--','color',[0.5 0.5 0.5])
plot([x_plot(rev_step) x_plot(rev_step)],[-1 1],'--','color',[0.5 0.5 0.5])
plot([success_trial success_trial],[-1 1],'--','color',[0.5 0.5 0.5])
xlim([min(x_plot) max([success_trial max(x_plot)])])
ylim([-0.4 0.3])
axisVal = axis;
plot([cfit_rise.a cfit_rise.a+log(19)/cfit_rise.k],[axisVal(4) axisVal(4)],'markersize',10,'marker','v','markerfacecolor','k','markeredgecolor','k','linestyle','none')
yticks(-0.3:0.3:0.3)
xlabel('Trial since reversal','fontsize',20)
ylabel('\itr','fontsize',20)
t = title('Sigmoid fitting','fontsize',20);
