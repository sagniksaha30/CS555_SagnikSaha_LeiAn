function p1b_burgers_slope_bdf3ext3()
% p1b_burgers_slope_bdf3ext3.m
% HW2 P1(b): compute s(t)=max_x |u_x(x,t)| at each timestep for 1D viscous Burgers
% u_t + u u_x = nu u_xx, u(0,t)=u(1,t)=0, u(x,0)=sin(pi x), x in [0,1]
% Space: 3-point FD (uniform grid, N=200)
% Time: IMEX BDF3/EXT3 (bootstrapped with BDF1/EXT1, BDF2/EXT2)
%
% Outputs:
%  - saves figure: p1b_burgers_slope_N200.png
%  - prints s* = max_{t in [0,T]} s(t) and the argmax time t*

clear; close all; clc;

% ---------------- user parameters ----------------
N   = 200;          % (use N=200 to match 1a baseline)
nu  = 0.01/pi;         % viscosity (set to your intended value)
T   = 2.0;
CFL = 0.10;         % chosen because advection is explicit (EXT), so dt should satisfy CFL-type constraint
% -------------------------------------------------

dx = 1/N;

% Main run with dt from CFL
dt0 = CFL*dx;
[t, s, sstar, tstar] = run_burgers_slope(N, nu, T, dt0);

% Verification (simple): rerun with dt/2 and compare s*, t*
[t2, s2, sstar2, tstar2] = run_burgers_slope(N, nu, T, dt0/2);

fprintf('\n=== P1(b) results (convective form) ===\n');
fprintf('N=%d, nu=%.6g, T=%.2f\n', N, nu, T);
fprintf('dt  = %.6g  ->  s* = %.10g at t* = %.10g\n', t(2)-t(1),  sstar,  tstar);
fprintf('dt/2= %.6g  ->  s* = %.10g at t* = %.10g\n', t2(2)-t2(1), sstar2, tstar2);

rel_diff_s = abs(sstar2 - sstar) / max(abs(sstar2), 1e-14);
abs_diff_t = abs(tstar2 - tstar);
fprintf('Verification: rel diff in s* = %.3e, abs diff in t* = %.3e\n', rel_diff_s, abs_diff_t);

% Plot (use the finer dt/2 curve to be safer)
fname = sprintf('p1b_burgers_slope_N%d.png', N);
figure('Position',[100 100 800 450]); hold on; box on;
plot(t2, s2, 'LineWidth', 2);
plot(tstar2, sstar2, 'o', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('time t');
ylabel('s(t) = max_x |u_x(x,t)|');
title(sprintf('Burgers slope metric (IMEX BDF3/EXT3), N=%d, nu=%.4g', N, nu));
grid on;

text(tstar2, sstar2, sprintf('  s* = %.3g at t* = %.3g', sstar2, tstar2), ...
     'VerticalAlignment','bottom', 'HorizontalAlignment','left');

print(gcf, fname, '-dpng', '-r200');
fprintf('Saved figure: %s\n', fname);

end

% ============================================================
function [t, s, sstar, tstar] = run_burgers_slope(N, nu, T, dt_in)
% Runs IMEX BDFk/EXTk (k=1,2 bootstrap, then k=3) and returns s(t)=max|u_x|.

dx = 1/N;

% enforce integer steps so final time hits T exactly
nsteps = max(1, round(T/dt_in));
dt = T/nsteps;

% grid
x = linspace(0,1,N+1).';
n = N-1;                   % interior dofs

% operators on interior
e = ones(n,1);

% A approximates (-nu u_xx), SPD
A = (nu/dx^2) * spdiags([-e 2*e -e],[-1 0 1], n, n);

% D approximates u_x (2nd-order central); boundary values are 0 so omission is fine
D = (1/(2*dx)) * spdiags([-e 0*e e],[-1 0 1], n, n);

I = speye(n);

% initial condition on interior
u_full0 = sin(pi*x);
u_n  = u_full0(2:end-1);
z_n  = nonlinear_convective(u_n, D);

% history placeholders
u_nm1 = u_n; z_nm1 = z_n;
u_nm2 = u_n; z_nm2 = z_n;

% store s(t)
t = (0:nsteps).' * dt;
s = zeros(nsteps+1,1);
s(1) = max(abs(D*u_n));

for step = 1:nsteps
    if step == 1
        k = 1;
    elseif step == 2
        k = 2;
    else
        k = 3;
    end

    if k==1
        b0=1;    b1=-1;  b2=0;    b3=0;
        a1=1;    a2=0;   a3=0;
    elseif k==2
        b0=1.5;  b1=-2;  b2=0.5;  b3=0;
        a1=2;    a2=-1;  a3=0;
    else
        b0=11/6; b1=-3;  b2=1.5;  b3=-1/3;
        a1=3;    a2=-3;  a3=1;
    end

    H = b0*I + dt*A;

    rhs = -(b1*u_n + b2*u_nm1 + b3*u_nm2) ...
          + dt*(a1*z_n + a2*z_nm1 + a3*z_nm2);

    u_np1 = H \ rhs;

    % shift history
    u_nm2 = u_nm1; z_nm2 = z_nm1;
    u_nm1 = u_n;   z_nm1 = z_n;
    u_n   = u_np1; z_n   = nonlinear_convective(u_n, D);

    s(step+1) = max(abs(D*u_n));
end

% max and argmax
[sstar, idx] = max(s);
tstar = t(idx);

end

% ============================================================
function z = nonlinear_convective(u, D)
% convective form: u_t + u u_x = nu u_xx  ->  u_t + A u = -(u .* (D u))
z = -(u .* (D*u));
end