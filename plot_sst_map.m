function plot_sst_map(fk_recon,sensors)
    load("sst_setup.mat",'mask','temp_bounds','x','sea_points')

    figure()
    axis([0 360 -90 90])
    set(gca,'FontSize',36)
    hold on
    maskedfk = fk_recon';
    maskedfk(logical(1-mask')) = NaN;
    h_2 = imagesc(x(:,2), x(:,1), maskedfk);
    colorbar
    % clim(temp_bounds) % for temperature values
    % colormap(sky) % for variance
    % colormap(cool)
    % set(gca, 'FontSize', 16)
    set(h_2, 'AlphaData', ~isnan(maskedfk))
    plot(sea_points(sensors,2), sea_points(sensors,1), 'r*', 'LineWidth',1,'MarkerSize',10)
    hold off
    xlabel('Latitude','Interpreter','Latex')
    ylabel('Longitude','Interpreter','Latex')
end