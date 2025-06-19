clc; clear; close all;

%% NUMERICAL SOLUTION =====================================================

n_max = 1e4;
N = 1e3;
r_min = 1+1e-2;
r_max = 4;
num_of_r = 2e2;
rs = linspace(r_min, r_max, num_of_r);
a = 0;
b = 1;

res.r = [];
res.mu_numerical = [];
res.mu_anal = [];
for r_index = 1:length(rs)
    r = rs(r_index);
    presenteg = r_index / length(rs) * 100;
    fprintf('r: %0.4f | n_max: %d | N: %d | done: %d%%\n', r, n_max, N, floor(presenteg))
    % numerical solution
    sum = 0;
    for realization = 1:N
        % get final u for specific r and unsamble
        u = [];
        u(1) = rand() * (b - a) + a;
        for i = 1:N
            u(i+1) = r * u(i) * (1 - u(i));
        end
        sum = sum + u(end);
    end
    mu_numerical = sum / N;
    res.r(end+1) = r;
    res.mu_numerical(end+1) = mu_numerical;
    
    syms mu sig_squ
    eq1 = mu^2 + sig_squ == r^2 / (1 - r^2) * (-2 * (mu^3 + 3 * mu * sig_squ) + mu^4 + 6 * mu^2 * sig_squ + 3 * sig_squ);
    eq2 = mu == r / (r - 1) * (mu^2 + sig_squ);

    anal_solution = solve([eq1, eq2], [mu, sig_squ]);
    real_mu = [];
    for i = 1:length(anal_solution.mu)
        if imag(anal_solution.mu(i)) ~= 0
            continue
        end
        real_mu(end+1) = double(anal_solution.mu(i));
    end
    % double(anal_solution.mu)
    % real_mu
    % mu_numerical
    res.mu_anal(end+1) = max(real_mu);
end

fig1 = figure('Name','1', 'Position', [0, 250, 900, 600]);
hold all
size = 20;

plot(res.r, res.mu_numerical,'LineStyle','-','LineWidth',1.5)
plot(res.r, res.mu_anal,'LineStyle','--','LineWidth',1.5)

title('$\langle u\rangle$ as a Function of r', 'FontSize', size,'Interpreter','latex')
ylabel('$\langle u\rangle$', 'FontSize', size,'Interpreter','latex')
xlabel('$r$','FontSize', size,'Interpreter', 'latex')
legend({'numerical solution', 'analytical solution'}, 'Location', 'northwest', 'FontSize', size-4, 'Interpreter','latex')
grid on
grid minor
box on
% exportgraphics(fig1, 'images/Q2.4.png','Resolution',400);

