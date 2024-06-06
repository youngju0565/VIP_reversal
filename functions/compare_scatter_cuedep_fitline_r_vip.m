function [fig_name, axis_name] = compare_scatter_cuedep_fitline_r_vip(data_draw, x,y, idx_RP, edgecolor)
    fig_name = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 9]);
    axis_name = axes;
    hold on
    if nargin < 3
        x=1; y=3; idx_RP=[]; edgecolor='none';
    end
    if nargin < 4
        idx_RP=[]; edgecolor='none';
    end
    if nargin < 5
        edgecolor='none';
    end
    marker_color = 0.4;
    scatter(data_draw(:,x), data_draw(:,y),'filled','markeredgecolor',[1 1 1]*marker_color,'markerfacecolor',[1 1 1]*marker_color)
    if ~isempty(idx_RP)
        scatter(data_draw(idx_RP,x), data_draw(idx_RP,y),'filled','r')
    end
    label_set = {'Before','During','After'};
    xlabel(label_set{x})
    ylabel(label_set{y})
    title('Cue-dep. activity')
    xlim([-1 1])
    ylim([-1 1])
    xline1 = line([-2000 2000],[0 0]);
    xline1.Color = 'k'; xline1.LineStyle = '--';
    yline1 = line([0 0],[-2000 2000]);
    yline1.Color = 'k'; yline1.LineStyle = '--';
    fit1 = polyfit(data_draw(:,x),data_draw(:,y),1);
    plot([min(data_draw(:,x)) max(data_draw(:,x))]*1.5,polyval(fit1,[min(data_draw(:,x)) max(data_draw(:,x))]*1.5),'Color','k','linewidth',1.5)
    [r,p] = corrcoef(data_draw(:,x),data_draw(:,y));
    if isequal(num2str(p(1,2),'%.3f'),'0.000')
        p_str = (num2str(p(1,2),'%.1e'));
        if isequal(p_str(end-1),'0')
            p_str_set = [p_str(1:3) '\times10^-^' p_str(end)];
            start_pos_x = 0.52;
        else
            p_str_set = [p_str(1:3) '\times10^-^' p_str(end-1) '^' p_str(end)];
            start_pos_x = 0.54;
        end
        str1 = {['{\itr} = ',num2str(r(1,2),'%.3f')],['{\itp} = ',p_str_set]};
        annotation('textbox',[0.2 0.3 0.15 0.1],'String',str1,'FitBoxToText','on','edgecolor',edgecolor,'fontsize',14)
        if p(1,2)<0.05; star_str = '*'; if p(1,2)<0.01; star_str = '**'; if p(1,2)<0.001; star_str = '***'; end; end; else; star_str = ''; end
        annotation('textbox',[start_pos_x 0.16 0.2 0.2],'string',star_str,'FitBoxToText','on','edgecolor','none','fontsize',30)%,'color',[200 0 0]/255)
    else        
        str1 = {['{\itr} = ',num2str(r(1,2),'%.3f')],['{\itp} = ',num2str(p(1,2),'%.3f')]};
        annotation('textbox',[0.2 0.3 0.15 0.1],'String',str1,'FitBoxToText','on','edgecolor',edgecolor,'fontsize',14)
        if p(1,2)<0.05; star_str = '*'; if p(1,2)<0.01; star_str = '**'; if p(1,2)<0.001; star_str = '***'; end; end; else; star_str = ''; end
        annotation('textbox',[0.43 0.16 0.2 0.2],'string',star_str,'FitBoxToText','on','edgecolor','none','fontsize',30)%,'color',[200 0 0]/255)
    end
end