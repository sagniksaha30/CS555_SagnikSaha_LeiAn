% p1a_burgers_bdf3ext3.m
% HW2 P1(a): 1D viscous Burgers with Dirichlet BCs on [0,1]
% u_t + u u_x = nu u_xx,  u(0,t)=u(1,t)=0,  u(x,0)=sin(pi x)
% Space: 3-pt FD (uniform grid, N=200)
% Time: IMEX BDF3/EXT3 (bootstrapped with BDF1/EXT1, BDF2/EXT2)
% Plot u(x,t) for t=0:0.1:2 on a single figure, save png.

clear; close all; clc;

% ---------------- user parameters ----------------
N  = 200;          % as required in 1a
nu = 0.01;         % viscosity (set this to the value you intend to use)
T  = 2.0;
CFL = 0.1;         % not required in 1a, but a safe choice (also used later in 1c)
% -------------------------------------------------

dx = 1/N;
dt = CFL*dx;                 % dt = 0.1 * (1/N)
nsteps = round(T/dt);        % make T hit exactly
dt = T/nsteps;

% grid including boundaries
x = linspace(0,1,N+1).';
n = N-1;                     % interior unknown count
xi = x(2:end-1);             % interior grid

% --- build A approximating (-nu * u_xx) on interior (SPD) ---
e = ones(n,1);
A = (nu/dx^2) * spdiags([-e 2*e -e],[-1 0 1], n, n);

% --- build D approximating u_x on interior (2nd-order central) ---
D = (1/(2*dx)) * spdiags([-e 0*e e],[-1 0 1], n, n);
% (boundary values u0=uN=0 are implicitly used because we only store interior)

I = speye(n);

% initial condition
u0_full = sin(pi*x);
u_n = u0_full(2:end-1);      % interior vector at t^0
z_n = nonlinear_convective(u_n, D);   % z^n = -u .* (D*u)

% storage for plot times
t_plot = (0:0.1:2.0).';
kplot  = numel(t_plot);
Uplot  = zeros(N+1, kplot);
Uplot(:,1) = u0_full;

plot_idx = round(t_plot/dt) + 1;  % indices in time grid (t=0 is idx=1)

% history placeholders (needed for k=2,3)
u_nm1 = u_n;   z_nm1 = z_n;
u_nm2 = u_n;   z_nm2 = z_n;

% time stepping
t = 0.0;
iplot = 2;  % next slot in Uplot (since t=0 stored)

for step = 1:nsteps
    % choose order k for bootstrapping
    if step == 1
        k = 1;  % BDF1/EXT1
    elseif step == 2
        k = 2;  % BDF2/EXT2
    else
        k = 3;  % BDF3/EXT3
    end

    % coefficients (from assignment snippet)
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

    % advance time
    t = t + dt;

    % update history (shift)
    u_nm2 = u_nm1;   z_nm2 = z_nm1;
    u_nm1 = u_n;     z_nm1 = z_n;
    u_n   = u_np1;   z_n   = nonlinear_convective(u_n, D);

    % save snapshot if needed
    if iplot <= kplot && (step+1) == plot_idx(iplot)
        u_full = zeros(N+1,1);
        u_full(2:end-1) = u_n;
        Uplot(:,iplot) = u_full;
        iplot = iplot + 1;
    end
end

% ---- plotting ----
fname = sprintf('p1a_burgers_uxt_N%d.png', N);

figure('Position',[100 100 900 500]); hold on; box on;
for j = 1:kplot
    plot(x, Uplot(:,j), 'LineWidth', 1.0, ...
         'DisplayName', sprintf('t=%.1f', t_plot(j)));
end
xlabel('x'); ylabel('u(x,t)');
title(sprintf('Burgers (IMEX BDF3/EXT3), N=%d, dt=%.4g, \\nu=%.4g', N, dt, nu));
legend('Location','eastoutside');
grid on;

print(gcf, fname, '-dpng', '-r200');
fprintf('Saved figure: %s\n', fname);

% ---- nonlinear term: convective form ----
function z = nonlinear_convective(u, D)
    z = -(u .* (D*u));   % z(u) = -u * u_x (pointwise)
end