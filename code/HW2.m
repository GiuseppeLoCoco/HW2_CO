%% HW2 Continuous Optimization

%% Step 1: Generate the data

rng(0, 'twister'); % fix random seed for reproducibility
n = 200;
m = 6*n;
sigma = 1e-3;
A = rand(m, n) / sqrt(n);
x_true = randn(n, 1);
y = abs(A*x_true) + sigma*randn(m, 1);
x0 = randn(n, 1);

%% Question 1 : Implement tCG

H_mat = [4, 0; 0, 1];
b = [2; -1];
delta = 0.5;

H_fun = @(v) H_mat * v;