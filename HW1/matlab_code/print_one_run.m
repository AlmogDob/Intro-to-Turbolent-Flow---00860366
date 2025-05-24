clc; clear; close all;

step_size = 1; % [m]
num_of_steps = 1e2;
num_of_walks = 1e4;
num_of_runs  = 1e0;

fig1 = figure('Name', '1', 'Position', [100, 250, 900, 600]);

for run = 1:num_of_runs

    clf
    hold all
    steps = [];

    [end_x_vec, end_y_vec] = one_run(0, 0, step_size, num_of_steps, num_of_walks);

    plot(end_x_vec, end_y_vec, 'x', 'LineWidth', 2, 'Color', 'k')
    plot(0,0, '^', 'LineWidth', 3, 'Color', 'r')
    
    ax = gca;
    MaxX = max(abs(ax.XLim));
    MaxY = max(abs(ax.YLim));
    axis([-MaxX MaxX -MaxY MaxY]);
    
    title(sprintf('run: %d', run))
    xlabel('x [m]','FontSize',14,'Interpreter','latex')
    ylabel('y [m]','FontSize',14,'Interpreter','latex')
    legend({'end of one walk','start'})
    box on
    grid on
    grid minor
    
    drawnow

    pause(0.5)
end












% FUNCTIONS ###############################################################
function [des_x, des_y] = one_step(src_x, src_y, step_size)
    theta = rand()*2*pi;
    des_x = src_x + step_size * cos(theta);
    des_y = src_y + step_size * sin(theta);    
end

function [end_x, end_y] = one_walk(start_x, start_y, step_size, num_of_steps)
    x = start_x;
    y = start_y;
    for i=1:num_of_steps
        [x,y] = one_step(x, y, step_size);
    end
    end_x = x;
    end_y = y;
end

function [end_x_vec, end_y_vec] = one_run(start_x, start_y, step_size, num_of_steps, num_of_walks)
    fprintf('preforming a run ...\n');
    steps = [];
    for i=1:num_of_walks
        [x,y] = one_walk(start_x, start_y, step_size, num_of_steps);
        steps(end+1,:) = [x,y];
    end
    end_x_vec = steps(:,1);
    end_y_vec = steps(:,2);
end

function count = num_points_in_band(r, delta_r, x_vec, y_vec)
    
end