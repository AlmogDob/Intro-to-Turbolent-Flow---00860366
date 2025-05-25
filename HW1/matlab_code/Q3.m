clc; clear; close all;

%% Q.3.a PREFORM AND PRINT A RUN #######################################
step_size = 1; % [m]
num_of_steps = 1e4;
num_of_walks = 1e3;
r  = 50;
dr = 4;

fig1 = figure('Name', '1', 'Position', [100, 250, 900, 600]);
hold all

steps = [];
[end_x_vec, end_y_vec, ~] = one_run(0, 0, step_size, num_of_steps, num_of_walks);

[count_of_band, count_of_outer_circ] = num_points_in_band_and_outer_circ(r, dr, end_x_vec, end_y_vec);

plot(end_x_vec, end_y_vec, 'x', 'LineWidth', 2, 'Color', 'k')
plot(0,0, '^', 'LineWidth', 2, 'Color', 'r')
draw_circ(0, 0, r, 'g', 2)
draw_circ(0, 0, r+dr, 'g', 2)

ax = gca;
MaxX = max(abs(ax.XLim));
MaxY = max(abs(ax.YLim));
axis([-MaxX MaxX -MaxY MaxY]);

xlabel('x [m]','FontSize',14,'Interpreter','latex')
ylabel('y [m]','FontSize',14,'Interpreter','latex')
legend({'end of one walk','start'})
axis equal
box on
grid on
grid minor
% exportgraphics(fig1, 'images/graph1.png','Resolution',300);

%% Q.3.a PREFORM AND PRINT A RUN AND PDF ###############################
step_size = 1; % [m]
num_of_steps     = 1e2;
dr               = 0.5;
num_of_walks_vec = [4e5, 3e5, 2e5, 1e5, 7.5e4, 5e4, 2.5e4, 1e4, 7.5e3, 5e3, 2.5e3, 1e3];
% num_of_walks_vec = [7.5e4, 5e4, 2.5e4, 1e4, 5e3, 1e3];
num_of_radiuses  = 1e2;

end_x_vec_vec        = {};
end_y_vec_vec        = {};
distances_vec_of_vec = {};
rs_vec               = {};
pdf_vec_vec          = {};
for index = 1:length(num_of_walks_vec)
    num_of_walks = num_of_walks_vec(index);
    steps = [];
    [current_end_x_vec, current_end_y_vec, current_distances_vec] = one_run(0, 0, step_size, num_of_steps, num_of_walks);
    end_x_vec_vec{index,1} = current_end_x_vec;
    end_y_vec_vec{index,1} = current_end_y_vec;
    distances_vec_of_vec{index,1} = current_distances_vec;
    rs_vec{index,1} = linspace(0, min(max(abs(current_end_x_vec)),max(abs(current_end_y_vec))), num_of_radiuses);
    pdf_vec_vec{index,1} = calc_pdf(rs_vec{index,1}, dr, current_end_x_vec, current_end_y_vec);
end

fig2 = figure('Name', '2', 'Position', [150, 250, 1500, 600]);
hold all
colors = jet(length(num_of_walks_vec));
% colors = cool(length(num_of_walks_vec))*0.8;
lg  = {};
rms = [];
for index = 1:length(num_of_walks_vec)
    num_of_walks = num_of_walks_vec(index);
    lg{end+1} = sprintf('%g', double(num_of_walks_vec(end-index+1)));
    end_x_vec = end_x_vec_vec{index,1};
    end_y_vec = end_y_vec_vec{index,1};

    subplot(1,3,1) % ###################################################
    hold all
    p = plot(0,0, 'x', 'LineWidth', 1, 'Color', 'k', 'HandleVisibility', 'callback');
    plot(end_x_vec, end_y_vec, 'x', 'LineWidth', 2, 'Color', colors(index,:))
    s = plot(0,0, '^', 'LineWidth', 2, 'Color', 'k');

    ax = gca;
    MaxX = max(abs(ax.XLim));
    MaxY = max(abs(ax.YLim));
    axis([-MaxX MaxX -MaxY MaxY]);
    xlabel('x [m]','FontSize',14,'Interpreter','latex')
    ylabel('y [m]','FontSize',14,'Interpreter','latex')
    title('End of Walks Distribution')
    legend([p, s], {'end of one walk','start'})
    axis equal
    % axis square
    box on
    grid on
    grid minor
    
    subplot(1,3,2) % ###################################################
    hold all
    
    rs                  = rs_vec{length(num_of_walks_vec)-index+1, 1};
    analytical_solution = 2/num_of_steps.*rs.*exp(-rs.^2/num_of_steps);
    pdf_vec             = pdf_vec_vec{length(num_of_walks_vec)-index+1, 1};
    
    plot(rs, pdf_vec, '*', 'LineWidth', 2, 'Color', colors(length(num_of_walks_vec)-index+1,:))
    if index == length(num_of_walks_vec)
        plot(rs, analytical_solution, '-', 'LineWidth', 2, 'Color', 'k')
        lg{end+1} = 'analytical solution';
    end
    
    xlabel('r [m]','FontSize',14,'Interpreter','latex')
    ylabel('pdf [-]','FontSize',14,'Interpreter','latex')
    title('PDF as a Function of r')
    legend(lg)
    box on
    grid on
    grid minor

    subplot(1,3,3) % ###################################################
    
    rms(end+1) = sqrt(sum((analytical_solution-pdf_vec).^2));
    
    semilogx(flip(num_of_walks_vec(end-length(rms)+1:end)), rms, '-', 'LineWidth', 2, 'Color', 'k')
    
    xlabel('\# realization [-]','FontSize',14,'Interpreter','latex')
    ylabel('RMS [-]','FontSize',14,'Interpreter','latex')
    title('RMS as a Function of # of realization')
    box on
    grid on
    grid minor
end
% exportgraphics(fig2, 'images/graph2.png','Resolution',300);
% exportgraphics(fig1, 'images/graph1.png','Resolution',300); exportgraphics(fig2, 'images/graph2.png','Resolution',300);

%% Q.3.b Average Displacement ##########################################
% step_size = 1; % [m]
% num_of_steps     = 1e2;
% dr               = 0.5;
% % num_of_walks_vec = [4e5, 3e5, 2e5, 1e5, 7.5e4, 5e4, 2.5e4, 1e4, 7.5e3, 5e3, 2.5e3, 1e3];
% num_of_walks_vec = [7.55e4, 5e4, 2.5e4, 1e4, 5e3, 1e3];
% num_of_radiuses  = 1e2;
% 
% end_x_vec_vec        = {};
% end_y_vec_vec        = {};
% distances_vec_of_vec = {};
% rs_vec               = {};
% pdf_vec_vec          = {};
% for index = 1:length(num_of_walks_vec)
%     num_of_walks = num_of_walks_vec(index);
%     steps = [];
%     [current_end_x_vec, current_end_y_vec, current_distances_vec] = one_run(0, 0, step_size, num_of_steps, num_of_walks);
%     end_x_vec_vec{index,1} = current_end_x_vec;
%     end_y_vec_vec{index,1} = current_end_y_vec;
%     distances_vec_of_vec{index,1} = current_distances_vec;
%     rs_vec{index,1} = linspace(0, min(max(abs(current_end_x_vec)),max(abs(current_end_y_vec))), num_of_radiuses);
%     pdf_vec_vec{index,1} = calc_pdf(rs_vec{index,1}, dr, current_end_x_vec, current_end_y_vec);
% end

ave_disp = [];
for index = 1:length(num_of_walks_vec)
    ave_disp(index) = trapz(rs_vec{index,:}, rs_vec{index,:}.*pdf_vec_vec{index,:});
end

fig3 = figure('Name', '3', 'Position', [300, 250, 900, 600]);

semilogx(num_of_walks_vec, ave_disp, '-', 'LineWidth', 2, 'Color', 'k')
    
xlabel('\# realization [-]','FontSize',14,'Interpreter','latex')
ylabel('average displacment [m]','FontSize',14,'Interpreter','latex')
title('Average Displacment as a Function of # of realization')
box on
grid on
grid minor
% exportgraphics(fig3, 'images/graph3.png','Resolution',300);
% exportgraphics(fig3, 'images/graph3.png','Resolution',300); exportgraphics(fig1, 'images/graph1.png','Resolution',300); exportgraphics(fig2, 'images/graph2.png','Resolution',300);

%% Q.3.c autocovariance and autocorrelation ############################
step_size = 1; % [m]
number_of_steps = 1e2;

fig4 = figure('Name', '4', 'Position', [450, 250, 1350, 400]);

subplot(1,3,1) % #######################################################

hold all
steps = [];
distances = [];
steps(1,:) = [0,0];
for i=1:number_of_steps
    [x,y] = one_step(steps(i, 1), steps(i, 2), step_size);
    steps(end+1,:) = [x,y];
    distances(end+1) = sqrt(x^2+y^2);
end

plot(steps(:,1), steps(:,2), 'x', 'LineWidth', 0.75, 'Color', 'k', 'HandleVisibility','off')
quiver(steps(1:end-1,1), steps(1:end-1,2), diff(steps(:,1)), diff(steps(:,2)), 'MarkerSize', 0.5, 'Color', 'k', 'MaxHeadSize', 0.05, 'HandleVisibility','off')
plot(steps(1,1),steps(1,2), '^', 'LineWidth', 3, 'Color', 'r')
plot(steps(end,1),steps(end,2), 'squar', 'LineWidth', 3, 'Color', 'g')

ax = gca;
MaxX = max(abs(ax.XLim));
MaxY = max(abs(ax.YLim));
axis([-MaxX MaxX -MaxY MaxY]);

xlabel('x position [m]','FontSize',14,'Interpreter','latex')
ylabel('y position [m]','FontSize',14,'Interpreter','latex')
title('One Random Walk')
legend({'start','end'})
axis square
box on
grid on
grid minor

Q = [];
for i = 1:length(distances)
    A = distances(1:length(distances)-i);
    B = distances(1+i:length(distances));
    C = cov(A, B);
    Q(i) = C(1,2);
end

subplot(1,3,2) % #######################################################
plot(1:length(distances), Q)

xlabel('\# of steps','FontSize',14,'Interpreter','latex')
ylabel('autocovariance [m]','FontSize',14,'Interpreter','latex')
title('Autocovariance as a Function of Steps')
box on
grid on
grid minor

subplot(1,3,3) % #######################################################
plot(1:length(distances), Q/Q(1))

xlabel('\# of steps','FontSize',14,'Interpreter','latex')
ylabel('autocovariance [m]','FontSize',14,'Interpreter','latex')
title('Autocovariance as a Function of Steps')
box on
grid on
grid minor
% exportgraphics(fig4, 'images/graph4.png','Resolution',300);
% exportgraphics(fig4, 'images/graph4.png','Resolution',300); exportgraphics(fig3, 'images/graph3.png','Resolution',300); exportgraphics(fig1, 'images/graph1.png','Resolution',300); exportgraphics(fig2, 'images/graph2.png','Resolution',300);

%% Q.3.c ensemble ######################################################
% step_size = 1; % [m]
% num_of_steps     = 1e2;
% dr               = 0.5;
% % num_of_walks_vec = [4e5, 3e5, 2e5, 1e5, 7.5e4, 5e4, 2.5e4, 1e4, 7.5e3, 5e3, 2.5e3, 1e3];
% num_of_walks_vec = [7.55e4, 5e4, 2.5e4, 1e4, 5e3, 1e3];
% num_of_radiuses  = 1e2;
% 
% end_x_vec_vec        = {};
% end_y_vec_vec        = {};
% distances_vec_of_vec = {};
% rs_vec               = {};
% pdf_vec_vec          = {};
% for index = 1:length(num_of_walks_vec)
%     num_of_walks = num_of_walks_vec(index);
%     steps = [];
%     [current_end_x_vec, current_end_y_vec, current_distances_vec] = one_run(0, 0, step_size, num_of_steps, num_of_walks);
%     end_x_vec_vec{index,1} = current_end_x_vec;
%     end_y_vec_vec{index,1} = current_end_y_vec;
%     distances_vec_of_vec{index,1} = current_distances_vec;
%     rs_vec{index,1} = linspace(0, min(max(abs(current_end_x_vec)),max(abs(current_end_y_vec))), num_of_radiuses);
%     pdf_vec_vec{index,1} = calc_pdf(rs_vec{index,1}, dr, current_end_x_vec, current_end_y_vec);
% end

Q_vec = {};
Q_ensemble = {};
for index = 1:length(num_of_walks_vec)
    fprintf('Walk number %d\n', index);
    distances_vec = distances_vec_of_vec{index,1};
    Q_ensemble{index,1} = zeros(1,length(distances));
    for j = 1:length(distances_vec)
        j
        distances = distances_vec{j,1};
        Q = [];
        for i = 1:length(distances)
            A = distances(1:length(distances)-i);
            B = distances(1+i:length(distances));
            C = cov(A, B);
            Q(i) = C(1,2);
        end
        Q_ensemble{index,1} = Q_ensemble{index,1} + Q;
    end
    Q_ensemble{index,1} = Q_ensemble{index,1} / num_of_walks_vec(index);
end
%%
fig5 = figure('Name', '5', 'Position', [600, 250, 900, 600]);

colors = jet(length(num_of_walks_vec));
% colors = cool(length(num_of_walks_vec))*0.8;
subplot(1,2,1) % #######################################################
hold all
lg  = {};
for index = 1:length(num_of_walks_vec)
    Q = Q_ensemble{index,1};
    plot(1:length(distances), Q, '-', 'LineWidth', 2, 'Color', colors(index,:))
    lg{end+1} = sprintf('# of realization: %d', num_of_walks_vec(index));
end

xlabel('\# of steps','FontSize',14,'Interpreter','latex')
ylabel('autocovariance [m]','FontSize',14,'Interpreter','latex')
title('Autocovariance as a Function of Steps - Full Ensemble')
legend(lg)
box on
grid on
grid minor

zoom = axes('position',[0.25 0.36 0.2 0.2]);
box on % put box around new pair of axes
hold all
for index = 1:length(num_of_walks_vec)
    Q = Q_ensemble{index,1};
    plot(1:length(distances), Q, '-', 'LineWidth', 2, 'Color', colors(index,:))
end
zoom.XLim = [27.8, 28.6];
zoom.YLim = [2.29, 2.36];
grid on
grid minor

subplot(1,2,2) % #######################################################
hold all
lg  = {};
for index = 1:length(num_of_walks_vec)
    Q = Q_ensemble{index,1};
    plot(1:length(distances), Q/Q(1), '-', 'LineWidth', 2, 'Color', colors(index,:))
    hold on
    lg{end+1} = sprintf('# of realization: %d', num_of_walks_vec(index));
end

xlabel('\# of steps','FontSize',14,'Interpreter','latex')
ylabel('autocovariance [m]','FontSize',14,'Interpreter','latex')
title('Autocovariance as a Function of Steps - Full Ensemble')
legend(lg)
box on
grid on
grid minor

zoom = axes('position',[0.69 0.36 0.2 0.2]);
box on % put box around new pair of axes
hold all
for index = 1:length(num_of_walks_vec)
    Q = Q_ensemble{index,1};
    plot(1:length(distances), Q/Q(1), '-', 'LineWidth', 2, 'Color', colors(index,:))
end
zoom.XLim = [27.8, 28.6];
zoom.YLim = [0.255, 0.258];
grid on
grid minor

% exportgraphics(fig5, 'images/graph5.png','Resolution',300);
% exportgraphics(fig5, 'images/graph5.png','Resolution',300); exportgraphics(fig4, 'images/graph4.png','Resolution',300); exportgraphics(fig3, 'images/graph3.png','Resolution',300); exportgraphics(fig1, 'images/graph1.png','Resolution',300); exportgraphics(fig2, 'images/graph2.png','Resolution',300);







% FUNCTIONS ############################################################
function [des_x, des_y] = one_step(src_x, src_y, step_size)
    theta = rand()*2*pi;
    des_x = src_x + step_size * cos(theta);
    des_y = src_y + step_size * sin(theta);    
end

function [end_x, end_y, distances] = one_walk(start_x, start_y, step_size, num_of_steps)
    x = start_x;
    y = start_y;
    distances = [];
    for i=1:num_of_steps
        [x,y] = one_step(x, y, step_size);
        distances(i) = sqrt(x^2+y^2);
    end
    end_x = x;
    end_y = y;
end

function [end_x_vec, end_y_vec, distances_vec] = one_run(start_x, start_y, step_size, num_of_steps, num_of_walks)
    fprintf('preforming a run ...\n');
    steps = [];
    distances_vec = {};
    for i=1:num_of_walks
        if ~mod(i, num_of_walks/10)
            fprintf('   complited: %2.0f%%\n', i/num_of_walks*100);
        end
        [x,y,distances] = one_walk(start_x, start_y, step_size, num_of_steps);
        steps(end+1,:) = [x,y];
        distances_vec{i, 1} = distances;
    end
    end_x_vec = steps(:,1);
    end_y_vec = steps(:,2);
end


function [count_of_band, count_of_outer_circ] = num_points_in_band_and_outer_circ(r, dr, x_vec, y_vec)
    if length(x_vec) ~= length(y_vec)
        fprintf('x and y are not the same length\n');
        return
    end

    count_of_inner_circ = 0;
    count_of_outer_circ = 0;
    for i = 1:length(x_vec)
        dist_squre = x_vec(i)^2+y_vec(i)^2;
        if dist_squre <= r^2
            count_of_inner_circ = count_of_inner_circ+1;
        end
        if dist_squre <= (r+dr)^2
            count_of_outer_circ = count_of_outer_circ+1;
        end
    end

    count_of_band = count_of_outer_circ - count_of_inner_circ;
end

function draw_circ(center_x, center_y, r, color, line_width)
    pos = [[center_x, center_y]-r, 2*r, 2*r];
    rectangle('Position',pos,'Curvature',[1 1], 'EdgeColor', color, 'LineWidth', line_width)
end

function pdf_vec = calc_pdf(rs, dr, x_vec, y_vec)
    if length(x_vec) ~= length(y_vec)
        fprintf('x and y are not the same length\n');
        return
    end
    fprintf('calc pdf, dr = %4f\n', dr);
    
    pdf_vec = [];
    for i=1:length(rs)
        r = rs(i);
        [count_of_band, count_of_outer_circ] = num_points_in_band_and_outer_circ(r, dr, x_vec, y_vec);
        cdf_inner = (count_of_outer_circ-count_of_band) / length(x_vec); % points inside circ
        cdf_outer = (count_of_outer_circ) / length(x_vec); % points outside circ
        pdf_vec(i) = (cdf_outer - cdf_inner) / dr;
    end
end
