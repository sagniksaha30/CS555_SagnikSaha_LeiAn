clear; clc; close all;

%% Parameters
L  = 1.0;
c  = 1.0;
nu = 1.0e-3;
f  = 0.0;            % here f = 0
E  = 100;            % number of uniform elements
N  = E + 1;          % number of nodes

T  = 1.5;            % final time
dt = 1.0e-3;         % time step
Nt = round(T/dt);

%% Uniform mesh
x = linspace(0, L, N).';
h = x(2) - x(1);

%% Assemble FEM matrices
M = sparse(N, N);    % mass matrix
A = sparse(N, N);    % diffusion stiffness matrix
C = sparse(N, N);    % advection matrix
F = zeros(N, 1);     % load vector (zero here)

for e = 1:E
    idx = [e, e+1];

    Me = (h/6) * [2, 1;
                  1, 2];

    Ae = (1/h) * [ 1, -1;
                  -1,  1];

    Ce = (c/2) * [-1,  1;
                  -1,  1];

    Fe = f * (h/2) * [1; 1];

    M(idx, idx) = M(idx, idx) + Me;
    A(idx, idx) = A(idx, idx) + Ae;
    C(idx, idx) = C(idx, idx) + Ce;
    F(idx)      = F(idx)      + Fe;
end

%% Boundary condition functions
gL = @(t) sin(pi * t);   % left Dirichlet boundary
% right boundary is Neumann: u_x(1,t)=0, so no gR needed

%% Initial condition
U0 = zeros(N,1);
U0(1) = gL(0.0);

%% ---- Step 1: IMEX Euler / EXT1 ----
t1 = dt;

K1   = (1/dt) * M + nu * A;
rhs1 = (1/dt) * (M * U0) - C * U0 + F;

[K1_bc, rhs1_bc] = apply_left_dirichlet_bc(K1, rhs1, t1, gL);
U1 = K1_bc \ rhs1_bc;

%% ---- Step 2: IMEX BDF2 / EXT2 ----
t2 = 2*dt;

K2   = (3/(2*dt)) * M + nu * A;
rhs2 = (2/dt) * (M * U1) ...
     - (1/(2*dt)) * (M * U0) ...
     - C * (2*U1 - U0) ...
     + F;

[K2_bc, rhs2_bc] = apply_left_dirichlet_bc(K2, rhs2, t2, gL);
U2 = K2_bc \ rhs2_bc;

%% ---- Main loop: IMEX BDF3 / EXT3 ----
Um2 = U0;   % U^{n-2}
Um1 = U1;   % U^{n-1}
Un  = U2;   % U^{n}

for n = 2:Nt-1
    tnp1 = (n+1) * dt;

    K = (11/(6*dt)) * M + nu * A;

    rhs = (3/dt)     * (M * Un) ...
        - (3/(2*dt)) * (M * Um1) ...
        + (1/(3*dt)) * (M * Um2) ...
        - C * (3*Un - 3*Um1 + Um2) ...
        + F;

    [K_bc, rhs_bc] = apply_left_dirichlet_bc(K, rhs, tnp1, gL);
    Up1 = K_bc \ rhs_bc;

    % Shift time levels
    Um2 = Um1;
    Um1 = Un;
    Un  = Up1;
end

U_final = Un;

%% Advection-dominated reference solution
% For pure advection: u(x,t) ≈ sin(pi(t-x))
u_adv_ref = sin(pi * (T - x));

%% Print diagnostics
fprintf('Problem 2(b): unsteady advection-diffusion with BDF3/EXT3\n');
fprintf('E = %d, dt = %.4e, T = %.2f\n', E, dt, T);
fprintf('Left boundary value at T: u(0,T) = %.6f\n', gL(T));
fprintf('Numerical min(U) = %.6f\n', min(U_final));
fprintf('Numerical max(U) = %.6f\n', max(U_final));
fprintf('Numerical value at x=1: u(1,T) = %.6f\n', U_final(end));

%% Plot 1: numerical solution and reference
figure;
plot(x, U_final, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
plot(x, u_adv_ref, 'k-', 'LineWidth', 1.8);
grid on;
xlabel('x');
ylabel('u(x,1.5)');
title('Problem 2(b): solution at t = 1.5, E = 100, right Neumann BC');
legend('FEM solution', 'Advection-dominated reference: sin(\pi(T-x))', ...
       'Location', 'best');

%% Plot 2: zoom near the right boundary
figure;
plot(x, U_final, 'ro-', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
plot(x, u_adv_ref, 'k-', 'LineWidth', 1.8);
grid on;
xlim([0.88, 1.0]);
xlabel('x');
ylabel('u(x,1.5)');
title('Problem 2(b): zoom near the right boundary');

%% Plot 3: FEM solution only
figure;
plot(x, U_final, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 12);
grid on;
xlabel('x');
ylabel('u(x,1.5)');
title('Problem 2(b): FEM solution only');

%% Optional: approximate derivative at x = 1
ux_right_approx = (U_final(end) - U_final(end-1)) / h;
fprintf('Approximate derivative at x=1: u_x(1,T) ≈ %.6e\n', ux_right_approx);

%% ------------------------------------------------------------
% Local function: impose only the LEFT Dirichlet BC strongly
%% ------------------------------------------------------------
function [Kbc, rhsbc] = apply_left_dirichlet_bc(K, rhs, t, gL)

    Kbc   = K;
    rhsbc = rhs;

    % left boundary node only
    Kbc(1, :) = 0;
    Kbc(1, 1) = 1;
    rhsbc(1)  = gL(t);
end