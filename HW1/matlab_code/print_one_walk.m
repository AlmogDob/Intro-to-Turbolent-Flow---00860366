clc; clear; close all;

step_size = 1; % [m]
number_of_steps = 5e3;
number_of_walks  = 1e2;

fig1 = figure('Name', '1', 'Position', [100, 250, 900, 600]);

for walk = 1:number_of_walks
    clf
    hold all
    steps = [];
    steps(1,:) = [0,0];
    for i=1:number_of_steps
        [x,y] = one_step(steps(i, 1), steps(i, 2), step_size);
        steps(end+1,:) = [x,y];
    end
    
    plot(steps(:,1), steps(:,2), 'x', 'Color', 'k', 'HandleVisibility','off')
    % quiver(steps(1:end-1,1), steps(1:end-1,2), diff(steps(:,1)), diff(steps(:,2)), 'MaxHeadSize', 0.5)
    plot(steps(1,1),steps(1,2), '^', 'LineWidth', 3, 'Color', 'r')
    plot(steps(end,1),steps(end,2), 'squar', 'LineWidth', 3, 'Color', 'g')

    title(sprintf('walk: %d', walk))
    legend({'start','end'})
    axis square
    ax = gca;
    MaxX = max(abs(ax.XLim));
    MaxY = max(abs(ax.YLim));
    axis([-MaxX MaxX -MaxY MaxY]);
    grid on
    drawnow

    % pause(0.5)
end




% FUNCTIONS ###############################################################
function [des_x, des_y] = one_step(src_x, src_y, step_size)
    theta = rand()*2*pi;
    des_x = src_x + step_size * cos(theta);
    des_y = src_y + step_size * sin(theta);    
end

