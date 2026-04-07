function plot_sst_map2(sst_true, sensors_initial, sensors_final)
    load("sst_setup.mat", 'mask', 'temp_bounds', 'x', 'sea_points')
    figure()
    
    % Set the axes background color to grey. This color will show 
    % through the transparent NaN values, forming the landmasses.
    set(gca, 'Color', [0.752 0.752 0.752]) 
    
    axis([0 360 -90 90])
    set(gca, 'FontSize', 24) 
    hold on
    
    % render background
    masked_data = sst_true'; 
    masked_data(logical(1-mask')) = NaN;
    
    h_2 = imagesc(x(:,2), x(:,1), masked_data);
    colorbar
    
    if exist('temp_bounds', 'var')
        clim(temp_bounds); 
    end
    
    set(h_2, 'AlphaData', ~isnan(masked_data))
    
    % partition sensors
    kept_sensors = intersect(sensors_initial, sensors_final);
    discarded_sensors = setdiff(sensors_initial, sensors_final);
    new_sensors = setdiff(sensors_final, sensors_initial);
    
    % overlay scatter plots
    if ~isempty(kept_sensors)
        plot(sea_points(kept_sensors,2), sea_points(kept_sensors,1), ...
            'co', 'LineWidth', 1.5, 'MarkerSize', 8);
    end
    
    if ~isempty(discarded_sensors)
        plot(sea_points(discarded_sensors,2), sea_points(discarded_sensors,1), ...
            'rx', 'LineWidth', 1.3, 'MarkerSize', 9);
    end
    
    if ~isempty(new_sensors)
        plot(sea_points(new_sensors,2), sea_points(new_sensors,1), ...
            'k*', 'LineWidth', 1.12, 'MarkerSize', 8);
    end
    
    hold off
    
    xlabel('Longitude', 'Interpreter', 'Latex')
    ylabel('Latitude', 'Interpreter', 'Latex')
end