% plot contour plot for droplet data

% plot true data over time as a contour plot
cmap = interp1([-100; 0; 100], [0 0 1; 1 1 1; 1 0 0], (-100:100));

[plot_t, plot_x] = meshgrid(t, x);
plot_reshape = @(u) reshape(u, size(y_all,1), size(y_all,2))';

% plot_flow(fk_post_save', plot_reshape, plot_t, plot_x)
plot_flow(y_all-fk_post_save', plot_reshape, plot_t, plot_x,cmap)
hold on
yline(x(pk_gks))
hold off
title('GKS-Selected Droplet Profile Error',Interpreter='latex')

%% plot all graphics as subplots
set(gca,'FontSize', 40)

cmap = interp1([-100; 0; 100], [0 0 1; 1 1 1; 1 0 0], (-100:100));

openfig('truecontour.fig')
ax1=gca;
figure;
tcl=tiledlayout(3,2);
ax1.Parent=tcl;
ax1.Layout.Tile=1;

f2 = openfig('gkscontour.fig');
ax2=f2.Children(2);
% ax2 = gca;
ax2.Parent=tcl;
ax2.Layout.Tile=3;

f3 = openfig('gkscontour_error.fig');
% colormap(cmap);
ax3=f3.Children(2);
% ax2 = gca;
colormap(ax3,cmap);
ax3.Parent=tcl;
ax3.Layout.Tile=4;

f4 = openfig('greedycontour.fig');
ax4=f4.Children(2);
% ax2 = gca;
ax4.Parent=tcl;
ax4.Layout.Tile=5;

f5 = openfig('greedycontour_error.fig');
% colormap(cmap);
ax5=f5.Children(2);
colormap(ax5,cmap);
% ax2 = gca;
ax5.Parent=tcl;
ax5.Layout.Tile=6;

%% plotting function

function [fig] = plot_flow(flow_data,plot_reshape, plot_t, plot_x,cmap)
    
    % Convert input into matrix
    % global plot_reshape plot_x plot_y cmap;
    if isvector(flow_data)
        flow_data = plot_reshape(flow_data);
    end
    % size(plot_t)
    % size(plot_x)
    % size(flow_data)
   
 

    % Plot pcolor of data
    figure('Position',[225 308 847 275]);
    fig = pcolor(plot_t, plot_x, flow_data');
    fig.FaceColor='interp';
    fig.EdgeAlpha=0;

    set(gca, "FontSize", 18)
    
    % Set colorbar and axes for data
    max_val = max(abs(flow_data), [], "all");
    clim([-max_val max_val]);
    colormap(cmap); 
    colorbar;
    % axis([-0.5 3 -0.5 0.5])
    xlabel("$t$", "Interpreter","latex");
    ylabel("$x$", "Interpreter","latex");
end