clear; clc; close all;

%% Parameters
L  = 1.0;
c  = 1.0;
nu = 1.0e-3;
f  = 1.0;
E  = 100;              % number of elements

%% Uniform mesh
x = linspace(0, L, E+1).';   % node coordinates
h = diff(x);                 % element sizes
N = E + 1;                   % total number of nodes

%% Assemble global matrices
A = sparse(N, N);            % diffusion stiffness matrix
C = sparse(N, N);            % advection matrix
F = zeros(N, 1);             % load vector

for e = 1:E
    he = h(e);
    idx = [e, e+1];          % global node indices of element e

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
I = 2:N-1;                   % interior DOFs
K = nu * A(I, I) + C(I, I);
rhs = F(I);

%% Solve linear system
u = zeros(N,1);
u(I) = K \ rhs;

%% Exact solution (stable formula)
% u(x) = (f/c) * [ x - (exp(-c*(L-x)/nu) - exp(-c*L/nu)) / (1 - exp(-c*L/nu)) ]
u_exact = (f/c) * ( x ...
          - (exp(-c*(L - x)/nu) - exp(-c*L/nu)) / (1 - exp(-c*L/nu)) );

%% Pointwise relative error at interior nodes
rel_err = abs(u(I) - u_exact(I)) ./ abs(u_exact(I));
[maxRelErr, kmax] = max(rel_err);
xmax = x(I(kmax));

fprintf('Maximum nodal pointwise relative error = %.12e\n', maxRelErr);
fprintf('Occurs at x = %.6f\n', xmax);
fprintf('u_h(x)     = %.12e\n', u(I(kmax)));
fprintf('u_exact(x) = %.12e\n', u_exact(I(kmax)));

%% Plot solution
figure;
plot(x, u_exact, 'k-', 'LineWidth', 2); hold on;
plot(x, u, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 4);
grid on;
xlabel('x');
ylabel('u(x)');
title('Problem (a): steady 1D advection-diffusion, linear FEM, E = 100');
legend('Exact solution', 'FEM solution', 'Location', 'best');

%% Plot pointwise relative error at interior nodes
figure;
plot(x(I), rel_err, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 14);
grid on;
xlabel('x');
ylabel('Relative error');
title('Pointwise relative error at interior nodes');