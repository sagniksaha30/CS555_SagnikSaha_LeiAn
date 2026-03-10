function p1d_conservation_form()
% Problem 1(d): Use conservative form for nonlinear term and report s* at finest resolution.
%
% PDE (HW2): u_t + u u_x = nu u_xx on [0,1], u(0,t)=u(1,t)=0, u(x,0)=sin(pi x)
% Time: IMEX BDF3/EXT3 (bootstrap BDF1/2)
% Space: 3-pt FD operators on either uniform or Chebyshev grid
%
% Compare nonlinear discretizations:
%   convective:   N(u) = -(u .* (D*u))
%   conservative: N(u) = -(D*(0.5*u.^2))
%
% Report s* = max_t s(t), s(t)=max_x |u_x|.

clear; close all; clc;

% --------- set these to match your 1(c) setup ----------
nu     = 1/(100*pi);   % recommended: same as paper; change if your HW specifies otherwise
Tfinal = 2.0;
CFL    = 0.1;
Nfinest = 2048;         % <-- set to the finest N you used in 1(c)
% ------------------------------------------------------

gridTypes = {'uniform','chebyshev'};

fprintf('=== Problem 1(d): conservative vs convective nonlinearity ===\n');
fprintf('nu=%.8g, T=%.2f, CFL=%.3g, N=%d\n\n', nu, Tfinal, CFL, Nfinest);

for g = 1:numel(gridTypes)
    gridType = gridTypes{g};

    % CFL rule requested: dx=1/N and |c|=1 for both uniform and Chebyshev
    dx_ref = 1.0/Nfinest;
    dt0    = CFL*dx_ref;
    nsteps = ceil(Tfinal/dt0);
    dt     = Tfinal/nsteps;

    [sstar_conv, tstar_conv] = burgers_sstar(Nfinest, nu, dt, nsteps, gridType, 'convective');
    [sstar_cons, tstar_cons] = burgers_sstar(Nfinest, nu, dt, nsteps, gridType, 'conservative');

    fprintf('--- grid: %s ---\n', gridType);
    fprintf('dt=%.3e, nsteps=%d\n', dt, nsteps);
    fprintf('convective   : s* = %.8g at t* = %.6g\n', sstar_conv, tstar_conv);
    fprintf('conservative : s* = %.8g at t* = %.6g\n\n', sstar_cons, tstar_cons);
end

end

% ============================================================
function [sstar, tstar] = burgers_sstar(N, nu, dt, nsteps, gridType, nlType)

    [xb, A, D] = build_ops(N, nu, gridType);
    n = N-1;
    I = speye(n);

    % IC: sin(pi x) on interior
    x_in = xb(2:end-1);
    u1 = sin(pi*x_in);  % u^n
    u2 = u1;            % u^{n-1} (dummy)
    u3 = u1;            % u^{n-2} (dummy)

    Nu1 = nonlinear(u1, D, nlType);
    Nu2 = Nu1;
    Nu3 = Nu1;

    tvals = (0:nsteps)*dt;
    svals = zeros(nsteps+1,1);
    svals(1) = max(abs(D*u1));

    for step = 1:nsteps
        k = min(step,3);
        [b0,b1,b2,b3,a1,a2,a3] = bdf_ext_coeffs(k);

        H = b0*I + dt*A;

        N_ext = a1*Nu1 + a2*Nu2 + a3*Nu3;
        rhs   = -(b1*u1 + b2*u2 + b3*u3) + dt*N_ext;

        unew = H \ rhs;

        % shift
        u3 = u2; u2 = u1; u1 = unew;
        Nu3 = Nu2; Nu2 = Nu1; Nu1 = nonlinear(u1, D, nlType);

        svals(step+1) = max(abs(D*u1));
    end

    [sstar, idx] = max(svals);
    tstar = tvals(idx);
end

% ============================================================
function Nu = nonlinear(u, D, nlType)
    switch lower(nlType)
        case 'convective'
            Nu = -(u .* (D*u));           % -(u u_x)
        case 'conservative'
            Nu = -(D*(0.5*(u.^2)));       % - ( (u^2)/2 )_x
        otherwise
            error('Unknown nlType="%s".', nlType);
    end
end

% ============================================================
function [b0,b1,b2,b3,a1,a2,a3] = bdf_ext_coeffs(k)
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
end

% ============================================================
function [xb, A, D] = build_ops(N, nu, gridType)
    j = (0:N)';
    if strcmpi(gridType,'uniform')
        xb = j/N;
    elseif strcmpi(gridType,'chebyshev')
        xb = (1 - cos(pi*j/N))/2;
    else
        error('Unknown gridType="%s".', gridType);
    end

    n = N-1;
    A = spalloc(n,n,3*n);
    D = spalloc(n,n,2*n);

    for i = 1:n
        jx = i+1; % interior index -> grid index
        dm = xb(jx)   - xb(jx-1);
        dp = xb(jx+1) - xb(jx);

        % A approximates (-nu u_xx), SPD tri-diagonal
        A(i,i) = 2*nu/(dm*dp);
        if i>1
            A(i,i-1) = -2*nu/(dm*(dm+dp));
            D(i,i-1) = -1/(dm+dp);        % centered: (w_{j+1}-w_{j-1})/(x_{j+1}-x_{j-1})
        end
        if i<n
            A(i,i+1) = -2*nu/(dp*(dm+dp));
            D(i,i+1) =  1/(dm+dp);
        end
    end
end