clear; clc; close all;

%% Parameters
L  = 1.0;
c  = 1.0;
nu = 1.0e-3;
f  = 1.0;
E  = 20;

%% ------------------------------------------------------------
% Part 1: trial-and-error with a few representative s values
%% ------------------------------------------------------------
s_try = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75];
err_try = zeros(size(s_try));

fprintf('Trial-and-error search for problem (c)\n');
fprintf('--------------------------------------\n');
fprintf('   s          max pointwise relative error\n');

for k = 1:length(s_try)
    s = s_try(k);
    [maxRelErr, ~, ~, ~, ~, ~, ~, ~] = solve_problem_c(E, s, L, c, nu, f);
    err_try(k) = maxRelErr;
    fprintf('%6.3f      %.12e\n', s, maxRelErr);
end

[bestErr_try, idx_try] = min(err_try);
bestS_try = s_try(idx_try);

fprintf('\nBest s among the trial values = %.3f\n', bestS_try);
fprintf('Corresponding max relative error = %.12e\n\n', bestErr_try);

%% ------------------------------------------------------------
% Part 2: refined search near the best trial value
%% ------------------------------------------------------------
% Since the best trial value is around 0.65, refine around it.
s_refined = 0.62:0.001:0.68;
err_refined = zeros(size(s_refined));

for k = 1:length(s_refined)
    s = s_refined(k);
    [maxRelErr, ~, ~, ~, ~, ~, ~, ~] = solve_problem_c(E, s, L, c, nu, f);
    err_refined(k) = maxRelErr;
end

[bestErr_refined, idx_refined] = min(err_refined);
bestS_refined = s_refined(idx_refined);

[maxRelErr, xmax, uh_at_max, uex_at_max, x_best, u_best, uex_best, h_best] = ...
    solve_problem_c(E, bestS_refined, L, c, nu, f);

fprintf('Refined search result\n');
fprintf('---------------------\n');
fprintf('Best s (refined) = %.6f\n', bestS_refined);
fprintf('Minimum max pointwise relative error = %.12e\n', bestErr_refined);
fprintf('Occurs at x = %.12f\n', xmax);
fprintf('u_h(x)     = %.12e\n', uh_at_max);
fprintf('u_exact(x) = %.12e\n\n', uex_at_max);

fprintf('Largest element size  h1 = %.12e\n', h_best(1));
fprintf('Smallest element size hE = %.12e\n', h_best(end));

%% ------------------------------------------------------------
% Plot 1: error vs s
%% ------------------------------------------------------------
figure;
plot(s_try, err_try, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
plot(s_refined, err_refined, 'b-', 'LineWidth', 1.5);
plot(bestS_refined, bestErr_refined, 'ks', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('s');
ylabel('Maximum pointwise relative error');
title('Problem (c): maximum pointwise relative error vs. grading factor s');
legend('Trial values', 'Refined search', 'Best s', 'Location', 'best');

%% ------------------------------------------------------------
% Plot 2: exact solution vs FEM solution at best s
%% ------------------------------------------------------------
figure;
plot(x_best, uex_best, 'k-', 'LineWidth', 2); hold on;
plot(x_best, u_best, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on;
xlabel('x');
ylabel('u(x)');
title(sprintf('Problem (c): best geometric mesh, E = %d, s = %.3f', E, bestS_refined));
legend('Exact solution', 'FEM solution', 'Location', 'best');

%% ------------------------------------------------------------
% Plot 3: pointwise relative error at best s
%% ------------------------------------------------------------
I = 2:length(x_best)-1;
rel_err_best = abs(u_best(I) - uex_best(I)) ./ abs(uex_best(I));

figure;
plot(x_best(I), rel_err_best, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 14);
grid on;
xlabel('x');
ylabel('Relative error');
title(sprintf('Problem (c): pointwise relative error at best s = %.3f', bestS_refined));

%% ------------------------------------------------------------
% Plot 4: mesh spacing for best s
%% ------------------------------------------------------------
figure;
stairs(1:E, h_best, 'LineWidth', 1.5);
grid on;
xlabel('Element index e');
ylabel('h_e');
title(sprintf('Problem (c): geometric element sizes for best s = %.3f', bestS_refined));

%% ============================================================
% Local function
%% ============================================================
function [maxRelErr, xmax, uh_at_max, uex_at_max, x, u, u_exact, h] = ...
         solve_problem_c(E, s, L, c, nu, f)

    % Construct geometric mesh
    if abs(s - 1.0) < 1e-14
        h = (L/E) * ones(E,1);
    else
        L1 = L * (1 - s) / (1 - s^E);
        h = zeros(E,1);
        for e = 1:E
            h(e) = L1 * s^(e-1);
        end
    end

    x = zeros(E+1,1);
    for e = 1:E
        x(e+1) = x(e) + h(e);
    end

    % Assemble global matrices
    N = E + 1;
    A = sparse(N, N);   % diffusion
    C = sparse(N, N);   % advection
    F = zeros(N, 1);    % load

    for e = 1:E
        he = h(e);
        idx = [e, e+1];

        Ae = (1/he) * [ 1, -1;
                       -1,  1];

        Ce = (c/2) * [-1,  1;
                      -1,  1];

        Fe = f * (he/2) * [1; 1];

        A(idx, idx) = A(idx, idx) + Ae;
        C(idx, idx) = C(idx, idx) + Ce;
        F(idx)      = F(idx)      + Fe;
    end

    % Homogeneous Dirichlet BC
    I = 2:N-1;
    K = nu * A(I, I) + C(I, I);
    rhs = F(I);

    u = zeros(N,1);
    u(I) = K \ rhs;

    % Exact solution (stable form)
    u_exact = (f/c) * ( x ...
              - (exp(-c*(L - x)/nu) - exp(-c*L/nu)) / (1 - exp(-c*L/nu)) );

    % Pointwise relative error at interior nodes
    rel_err = abs(u(I) - u_exact(I)) ./ abs(u_exact(I));
    [maxRelErr, kmax] = max(rel_err);

    xmax = x(I(kmax));
    uh_at_max = u(I(kmax));
    uex_at_max = u_exact(I(kmax));
end