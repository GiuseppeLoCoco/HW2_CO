%% HW2 Continuous Optimization
clear; clc; close all;

%% Setup Output Folder for Figures

code_dir = fileparts(mfilename('fullpath'));
fig_dir = fullfile(code_dir, '..', 'results');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% Step 1: Generate the data
rng(0, 'twister'); % Fix random seed for global reproducibility
n = 200;
m = 6 * n;
sigma = 1e-3;
A = rand(m, n) / sqrt(n);
x_true = randn(n, 1);
y = abs(A * x_true) + sigma * randn(m, 1);
x0 = randn(n, 1);

%% Question 1 : Implement tCG
disp('Running tCG tests...');

%%%%% 1D Positive Definite (Interior Solution)
H1 = 4;
g1 = 2;
delta1 = 5;
H_fun1 = @(v) H1 * v;
v0_1 = 0;
max_it = 100;
v1 = tCG(H_fun1, -g1, delta1, max_it);
v1_star = -g1 / H1; 
assert(norm(v1 - v1_star) < 1e-5, 'Test 1 Failed: Incorrect interior solution for 1D.');
assert(norm(v1) < delta1, 'Test 1 Failed: Solution incorrectly hit the boundary.');
disp('Test 1 (1D Interior Solution) PASSED.');

%%%% 2D Positive Definite (Interior Solution)
H2 = [2 0; 0 4];
g2 = [1; 2];
delta2 = 5;
H_fun2 = @(v) H2 * v;
v0_2 = zeros(2,1);
v2 = tCG(H_fun2, -g2, delta2, max_it);
v2_star = -H2 \ g2; 
assert(norm(v2 - v2_star) < 1e-5, 'Test 2 Failed: Incorrect interior solution for 2D.');
assert(norm(v2) < delta2, 'Test 2 Failed: Solution exceeded the trust region radius.');
disp('Test 2 (2D Interior Solution) PASSED.');

%%%%% 2D Positive Definite (Boundary Hit / Out-of-bounds Solution)
H3 = [2 0; 0 2];
g3 = [-4; 0];
delta3 = 1;
H_fun3 = @(v) H3 * v;
v0_3 = zeros(2,1);
v3 = tCG(H_fun3, -g3, delta3, max_it);
v3_expected = [1; 0]; 
assert(abs(norm(v3) - delta3) < 1e-5, 'Test 3 Failed: Solution is not on the trust region boundary.');
assert(norm(v3 - v3_expected) < 1e-5, 'Test 3 Failed: Solution vector direction is incorrect.');
disp('Test 3 (2D Boundary Hit) PASSED.');

%%%%% 2D Initial Negative Curvature (p0' * H * p0 <= 0)
H4 = [-2 0; 0 2];
g4 = [1; 0];
delta4 = 1;
H_fun4 = @(v) H4 * v;
v0_4 = zeros(2,1);
v4 = tCG(H_fun4, -g4, delta4, max_it);
v4_expected = [-1; 0]; 
assert(abs(norm(v4) - delta4) < 1e-5, 'Test 4 Failed: Solution did not hit boundary on negative curvature.');
assert(norm(v4 - v4_expected) < 1e-5, 'Test 4 Failed: Incorrect vector direction on negative curvature.');
disp('Test 4 (Initial Negative Curvature) PASSED.');

%%%%% 3D Indefinite Matrix (Negative Curvature in Later Iterations)
H5 = [1 0 0; 0 -2 0; 0 0 3];
g5 = [1; 1; 0];
delta5 = 2;
H_fun5 = @(v) H5 * v;
v0_5 = zeros(3,1);
v5 = tCG(H_fun5, -g5, delta5, max_it);
assert(abs(norm(v5) - delta5) < 1e-5, 'Test 5 Failed: Solution did not terminate on boundary.');
disp('Test 5 (3D Indefinite Matrix) PASSED.');

%%%%% Zero Gradient Corner Case (g = 0)
H6 = [2 1; 1 3];
g6 = [0; 0];
delta6 = 1;
H_fun6 = @(v) H6 * v;
v0_6 = zeros(2,1);
v6 = tCG(H_fun6, -g6, delta6, max_it);
assert(norm(v6) < 1e-8, 'Test 6 Failed: Solution should remain zero vector when g = 0.');
disp('Test 6 (Zero Gradient Corner Case) PASSED.');

disp('===================================================');
disp('All tCG tests completed successfully!');
disp('===================================================');

%% Question 2: Implement the trust region method
delta_0   = 1.0;
rho_prime = 0.1;
delta_max = 5.0;
max_it    = 50000;

%%%% Test 1: 2D Rosenbrock Benchmark (Official Functions)
n1 = 2;
x0_1 = [-1.2; 1.0];
fun1   = @(x) rosenbrockf(x);
b_fun1 = @(x) -rosenbrockgrad(x); 
H_fun1 = @(x, u) rosenbrockhess(x, u);
x1_opt = TR_tCG(x0_1, delta_0, fun1, H_fun1, b_fun1, rho_prime, delta_max, max_it);
x1_star = ones(n1, 1);
assert(norm(x1_opt - x1_star) < 1e-2, 'Test 1 Failed: 2D Rosenbrock did not converge to [1; 1].');
disp('Test 1 (2D Rosenbrock) PASSED.');

%%%% Test 2: Higher-Dimensional Rosenbrock (n = 5)
n2 = 5;
x0_2 = [-1.2; 1.0; -1.2; 1.0; 1.0]; 
fun2   = @(x) rosenbrockf(x);
b_fun2 = @(x) -rosenbrockgrad(x);
H_fun2 = @(x, u) rosenbrockhess(x, u);
x2_opt = TR_tCG(x0_2, delta_0, fun2, H_fun2, b_fun2, rho_prime, delta_max, max_it);
x2_star = ones(n2, 1);
assert(norm(x2_opt - x2_star) < 1e-2, 'Test 2 Failed: 5D Rosenbrock did not converge to ones.');
disp('Test 2 (5D Rosenbrock) PASSED.');

%%%% Test 3: Corner Case - Starting directly at the exact minimizer (n = 3)
n3 = 3;
x0_3 = ones(n3, 1);
fun3   = @(x) rosenbrockf(x);
b_fun3 = @(x) -rosenbrockgrad(x);
H_fun3 = @(x, u) rosenbrockhess(x, u);
x3_opt = TR_tCG(x0_3, delta_0, fun3, H_fun3, b_fun3, rho_prime, delta_max, max_it);
assert(norm(x3_opt - x0_3) < 1e-6, 'Test 3 Failed: Algorithm moved away from the exact minimum.');
disp('Test 3 (Already at Minimum) PASSED.');

%%%% Test 4: Convex 2D Quadratic Function
fun4   = @(x) (x(1) - 3)^2 + 2*(x(2) + 1)^2;
b_fun4 = @(x) -[2*(x(1) - 3); 4*(x(2) + 1)];
H_fun4 = @(x, u) [2 0; 0 4] * u;
x0_4   = [0; 0];
x4_opt = TR_tCG(x0_4, delta_0, fun4, H_fun4, b_fun4, rho_prime, delta_max, max_it);
x4_star = [3; -1];
assert(norm(x4_opt - x4_star) < 1e-4, 'Test 4 Failed: Simple quadratic optimization failed.');
disp('Test 4 (Convex Quadratic) PASSED.');

disp('===================================================');
disp('All Question 2 (TR-tCG) tests completed successfully!');
disp('===================================================');

%% Question 3: Test Gradient and Hessian-Vector product accuracy
x_test = randn(n, 1);
v_test = randn(n, 1);
v_test = v_test / norm(v_test);
fx = phase_retrieval_f(x_test, A, y);
gradfx = phase_retrieval_grad(x_test, A, y);
hessfxv = phase_retrieval_hess(x_test, v_test, A, y);

t_vals = logspace(-8, -1, 50);
err_grad = zeros(size(t_vals));
err_hess = zeros(size(t_vals));

for k = 1:length(t_vals)
    t = t_vals(k);
    fx_tv = phase_retrieval_f(x_test + t * v_test, A, y);
    
    err_grad(k) = abs(fx_tv - fx - t * (gradfx' * v_test));
    err_hess(k) = abs(fx_tv - fx - t * (gradfx' * v_test) - 0.5 * t^2 * (v_test' * hessfxv));
end

idx = (t_vals >= 1e-4) & (t_vals <= 1e-1);
p_grad = polyfit(log10(t_vals(idx)), log10(err_grad(idx)), 1);
p_hess = polyfit(log10(t_vals(idx)), log10(err_hess(idx)), 1);

fprintf('Gradient Taylor Error Slope (Expected ~2.0): %.4f\n', p_grad(1));
fprintf('Hessian Taylor Error Slope  (Expected ~3.0): %.4f\n', p_hess(1));

idx_anchor = find(idx, 1, 'last'); 
ref_O2 = err_grad(idx_anchor) * (t_vals / t_vals(idx_anchor)).^2;
ref_O3 = err_hess(idx_anchor) * (t_vals / t_vals(idx_anchor)).^3;

f1 = figure(1);
set(f1, 'Name', 'Question 3: Taylor Derivatives Check', 'Position', [100, 100, 550, 400]);
movegui(f1, 'center');
loglog(t_vals, err_grad, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Gradient error E_1(t)');
hold on;
loglog(t_vals, err_hess, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Hessian error E_2(t)');
loglog(t_vals, ref_O2, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Reference O(t^2) slope 2');
loglog(t_vals, ref_O3, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Reference O(t^3) slope 3');
grid on;
xlabel('Step size t');
ylabel('Absolute Taylor Error');
title('Taylor Expansion Derivative Check');
legend('Location', 'northwest');
drawnow;
exportgraphics(f1, fullfile(fig_dir, 'taylor_check_plot.pdf'));
exportgraphics(f1, fullfile(fig_dir, 'taylor_check_plot.png'));

%% Question 4: Run TR-tCG on Real Phase Retrieval Data
disp('===================================================');
disp('Running TR-tCG on Real Phase Retrieval dataset...');
fun_pr   = @(x) phase_retrieval_f(x, A, y);
b_fun_pr = @(x) -phase_retrieval_grad(x, A, y);
H_fun_pr = @(x, u) phase_retrieval_hess(x, u, A, y);
delta_0_pr   = 1.0;
rho_prime_pr = 0.1;
delta_max_pr = 5.0;
max_it_pr    = 50000;

[x_opt, history] = TR_tCG(x0, delta_0_pr, fun_pr, H_fun_pr, b_fun_pr, ...
                          rho_prime_pr, delta_max_pr, max_it_pr, x_true);
disp('TR-tCG Phase Retrieval completed.');

f2 = figure(2);
set(f2, 'Name', 'Question 4: TR-tCG Standalone Performance', 'Position', [100, 100, 750, 350]);
movegui(f2, 'center');
subplot(1, 2, 1);
semilogy(0:length(history.grad_norm)-1, history.grad_norm, 'b-o', ...
         'LineWidth', 1.5, 'MarkerSize', 3);
grid on;
xlabel('Iteration k');
ylabel('||\nabla f(x_k)||_2');
title('Gradient Norm Convergence');

subplot(1, 2, 2);
semilogy(0:length(history.recon_err)-1, history.recon_err, 'r-s', ...
         'LineWidth', 1.5, 'MarkerSize', 3);
grid on;
xlabel('Iteration k');
ylabel('min(||x_k - x_{true}||, ||x_k + x_{true}||)');
title('Reconstruction Error');
drawnow;
exportgraphics(f2, fullfile(fig_dir, 'tr_tcg_performance.pdf'));
exportgraphics(f2, fullfile(fig_dir, 'tr_tcg_performance.png'));

%% Question 4: Run TR-tCG and GD Line-Search Comparison
disp('===================================================');
disp('Running TR-tCG and GD (Line-Search) on Phase Retrieval dataset...');
grad_pr  = @(x) phase_retrieval_grad(x, A, y);
tol_pr   = 1e-5;

t_start_tr = tic;
[x_opt_tr, history_tr] = TR_tCG(x0, delta_0_pr, fun_pr, H_fun_pr, b_fun_pr, ...
                                rho_prime_pr, delta_max_pr, max_it_pr, x_true);
time_tr = toc(t_start_tr);

alpha_bar_pr = 1.0;
t_start_gd = tic;
[x_opt_gd, history_gd] = GD_linesearch(x0, fun_pr, grad_pr, ...
                                       alpha_bar_pr, tol_pr, max_it_pr, x_true);
time_gd = toc(t_start_gd);

fprintf('\n--- Computational Summary ---\n');
fprintf('TR-tCG:  %d iterations, Time: %.4f s\n', length(history_tr.grad_norm)-1, time_tr);
fprintf('GD (LS): %d iterations, Time: %.4f s\n', length(history_gd.grad_norm)-1, time_gd);
fprintf('===================================================\n');

f3 = figure(3);
set(f3, 'Name', 'Question 4: TR-tCG vs GD Line-Search Performance', 'Position', [100, 100, 800, 550]);
movegui(f3, 'center');

% Top-Left: TR-tCG Gradient Norm
subplot(2, 2, 1);
semilogy(0:length(history_tr.grad_norm)-1, history_tr.grad_norm, 'b-o', ...
    'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
xlabel('Iteration k');
ylabel('||\nabla f(x_k)||_2');
title('TR-tCG: Gradient Norm Convergence');

% Top-Right: TR-tCG Reconstruction Error
subplot(2, 2, 2);
semilogy(0:length(history_tr.recon_err)-1, history_tr.recon_err, 'b-o', ...
    'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
xlabel('Iteration k');
ylabel('min(||x_k - x_{true}||, ||x_k + x_{true}||)');
title('TR-tCG: Reconstruction Error');

% Bottom-Left: GD (Line-Search) Gradient Norm
subplot(2, 2, 3);
semilogy(0:length(history_gd.grad_norm)-1, history_gd.grad_norm, 'r-s', ...
    'LineWidth', 1.2, 'MarkerSize', 3);
grid on;
xlabel('Iteration k');
ylabel('||\nabla f(x_k)||_2');
title('GD (Line-Search): Gradient Norm Convergence');

% Bottom-Right: GD (Line-Search) Reconstruction Error
subplot(2, 2, 4);
semilogy(0:length(history_gd.recon_err)-1, history_gd.recon_err, 'r-s', ...
    'LineWidth', 1.2, 'MarkerSize', 3);
grid on;
xlabel('Iteration k');
ylabel('min(||x_k - x_{true}||, ||x_k + x_{true}||)');
title('GD (Line-Search): Reconstruction Error');

drawnow;
exportgraphics(f3, fullfile(fig_dir, 'tr_vs_gd_comparison.pdf'));
exportgraphics(f3, fullfile(fig_dir, 'tr_vs_gd_comparison.png'));

%% Question 4: Test GD Near x_true Initialization
disp('===================================================');
disp('Running GD Line-Search Near-Optimal Initialization ');

% Perturb x_true to start inside its local basin of attraction
perturbation_scale = 0.1;
x0_near = x_true + perturbation_scale * randn(n, 1);

% Run GD Line-Search from x0_near
t_start_gd_near = tic;
[x_opt_gd_near, history_gd_near] = GD_linesearch(x0_near, fun_pr, grad_pr, ...
    alpha_bar_pr, tol_pr, max_it_pr, x_true);
time_gd_near = toc(t_start_gd_near);

fprintf('\n--- GD Near-Initialization Summary ---\n');
fprintf('GD LS (Near): %d iterations, Final Err: %.4e, Time: %.4f s\n', ...
    length(history_gd_near.grad_norm)-1, history_gd_near.recon_err(end), time_gd_near);
fprintf('===================================================\n');

% Plot Figure 4: GD Near-Initialization Convergence
f4 = figure(4);
set(f4, 'Name', 'Question 4: GD Near x_true Initialization', 'Position', [100, 100, 750, 350]);
movegui(f4, 'center');
subplot(1, 2, 1);
semilogy(0:length(history_gd_near.grad_norm)-1, history_gd_near.grad_norm, 'r-s', ...
    'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'GD (Line-Search, Near)');
grid on;
xlabel('Iteration k');
ylabel('||\nabla f(x_k)||_2');
title('GD Gradient Norm (Near x_{true})');
legend('Location', 'northeast');

subplot(1, 2, 2);
semilogy(0:length(history_gd_near.recon_err)-1, history_gd_near.recon_err, 'r-s', ...
    'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'GD (Line-Search, Near)');
grid on;
xlabel('Iteration k');
ylabel('min(||x_k - x_{true}||, ||x_k + x_{true}||)');
title('GD Reconstruction Error (Near x_{true})');
legend('Location', 'northeast');
drawnow;
exportgraphics(f4, fullfile(fig_dir, 'gd_near_init.pdf'));
exportgraphics(f4, fullfile(fig_dir, 'gd_near_init.png'));

