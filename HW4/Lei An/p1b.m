clear; clc; close all;

%% Parameters
L  = 1.0;
c  = 1.0;
nu = 1.0e-3;
f  = 1.0;
E  = 20;          % number of elements
s  = 0.7;         % geometric scale factor

%% Construct nonuniform geometric mesh
% Element lengths: L_e = L1 * s^(e-1), e=1,...,E
% with sum_{e=1}^E L_e = L
L1 = L * (1 - s) / (1 - s^E);

h = zeros(E,1);
for e = 1:E
    h(e) = L1 * s^(e-1);
end

x = zeros(E+1,1);
for e = 1:E
    x(e+1) = x(e) + h(e);
end

%% Assemble global matrices
N = E + 1;                 % total number of nodes
A = sparse(N, N);          % diffusion stiffness matrix
C = sparse(N, N);          % advection matrix
F = zeros(N, 1);           % load vector

for e = 1:E
    he = h(e);
    idx = [e, e+1];        % global node indices of element e

    % Local matrices for linear FEM
    Ae = (1/he) * [ 1, -1;
                   -1,  1];

    Ce = (c/2) * [-1,  1;
                  -1,  1];

    Fe = f * (he/2) * [1; 1];

    % Assembly
    A(idx, idx) = A(idx, idx) + Ae;
    C(idx, idx) = C(idx, idx) + Ce;
    F(idx)      = F(idx)      + Fe;
end

%% Apply homogeneous Dirichlet BC: u(0)=u(1)=0
I = 2:N-1;                 % interior DOFs
K = nu * A(I, I) + C(I, I);
rhs = F(I);

%% Solve linear system
u = zeros(N,1);
u(I) = K \ rhs;

%% Exact solution
% Stable formula:
% u(x) = (f/c) * [ x - (exp(-c*(L-x)/nu) - exp(-c*L/nu)) / (1-exp(-c*L/nu)) ]
u_exact = (f/c) * ( x ...
          - (exp(-c*(L - x)/nu) - exp(-c*L/nu)) / (1 - exp(-c*L/nu)) );

%% Pointwise relative error at interior nodes
rel_err = abs(u(I) - u_exact(I)) ./ abs(u_exact(I));
[maxRelErr, kmax] = max(rel_err);
xmax = x(I(kmax));

%% Print results
fprintf('Problem (b): nonuniform geometric mesh, E = %d, s = %.2f\n', E, s);
fprintf('Largest element size  h1   = %.12e\n', h(1));
fprintf('Smallest element size hE   = %.12e\n', h(end));
fprintf('Maximum nodal pointwise relative error = %.12e\n', maxRelErr);
fprintf('Occurs at x = %.12f\n', xmax);
fprintf('u_h(x)     = %.12e\n', u(I(kmax)));
fprintf('u_exact(x) = %.12e\n', u_exact(I(kmax)));

%% Plot solution
figure;
plot(x, u_exact, 'k-', 'LineWidth', 2); hold on;
plot(x, u, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on;
xlabel('x');
ylabel('u(x)');
title('Problem (b): steady 1D advection-diffusion, geometric mesh, E = 20, s = 0.7');
legend('Exact solution', 'FEM solution', 'Location', 'best');

%% Plot pointwise relative error
figure;
plot(x(I), rel_err, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 14);
grid on;
xlabel('x');
ylabel('Relative error');
title('Pointwise relative error at interior nodes');

%% Plot mesh spacing distribution (optional but helpful)
figure;
stairs(1:E, h, 'LineWidth', 1.5);
grid on;
xlabel('Element index e');
ylabel('h_e');
title('Geometric element sizes for problem (b)');