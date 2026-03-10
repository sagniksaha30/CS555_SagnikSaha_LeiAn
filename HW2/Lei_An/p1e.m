function p1e_four_cases_tables()
% Problem 1(e):
% Produce tables for N=2^k, k=5..11 containing:
%   N, nsteps, s*, rel.err., ratio
% where rel.err = |s - s*|/|s| and s is analytical max|u_x| (Table 3 in [1]).
% ratio := err_{N/2} / err_{N}, expected ~4 for 2nd-order spatial convergence.
%
% Four cases:
%  (i)  Uniform spacing,   convective form
%  (ii) Chebyshev spacing, convective form
%  (iii)Uniform spacing,   conservation form
%  (iv) Chebyshev spacing, conservation form
%
% Time stepping: IMEX BDF3/EXT3 with BDF1/2 bootstrap.
% CFL fixed to 0.1 with dx=1/N and |c|=1 per assignment.

clear; close all; clc;

% ---------------- user settings ----------------
nu     = 1/(100*pi);      % recommended (paper-like)
Tfinal = 2.0;
CFL    = 0.1;
klist  = 5:11;
Nlist  = 2.^klist;

% Analytical s = max|u_x| from Table 3 in [1] (Basdevant et al.)
s_exact = 152.00516;   % <-- change if your handout specifies a different value
% -----------------------------------------------

cases = { ...
  struct('grid','uniform',   'nl','convective',   'tag','uni_conv',  'caption','Uniform spacing, convective form'), ...
  struct('grid','chebyshev', 'nl','convective',   'tag','cheb_conv', 'caption','Chebyshev spacing, convective form'), ...
  struct('grid','uniform',   'nl','conservative', 'tag','uni_cons',  'caption','Uniform spacing, conservation form'), ...
  struct('grid','chebyshev', 'nl','conservative', 'tag','cheb_cons', 'caption','Chebyshev spacing, conservation form') ...
};

fprintf('=== Problem 1(e): N=2^k, k=5..11, CFL=0.1, nu=%.8g, T=%.2f ===\n', nu, Tfinal);
fprintf('Analytical s = %.8f\n\n', s_exact);

for c = 1:numel(cases)
    C = cases{c};
    R = run_case(Nlist, nu, Tfinal, CFL, C.grid, C.nl, s_exact);

    % print to screen
    fprintf('--- Case %d: %s ---\n', c, C.caption);
    fprintf('%8s %8s %14s %12s %10s\n', 'N','nsteps','s*','rel.err','ratio');
    for i=1:numel(R.N)
        if isnan(R.ratio(i))
            fprintf('%8d %8d %14.8f %12.3e %10s\n', R.N(i), R.nsteps(i), R.sstar(i), R.relerr(i), '---');
        else
            fprintf('%8d %8d %14.8f %12.3e %10.3f\n', R.N(i), R.nsteps(i), R.sstar(i), R.relerr(i), R.ratio(i));
        end
    end
    fprintf('\n');

    % write LaTeX table
    texname = sprintf('p1e_%s.tex', C.tag);
    write_latex_table(texname, R, C.caption);
end

fprintf('Wrote LaTeX tables:\n');
for c = 1:numel(cases)
    fprintf('  p1e_%s.tex\n', cases{c}.tag);
end

end

% ============================================================
function R = run_case(Nlist, nu, Tfinal, CFL, gridType, nlType, s_exact)

m = numel(Nlist);

R.N      = zeros(m,1);
R.nsteps = zeros(m,1);
R.dt     = zeros(m,1);
R.sstar  = zeros(m,1);
R.tstar  = zeros(m,1);
R.relerr = zeros(m,1);
R.ratio  = nan(m,1);

prev_err = NaN;

for k = 1:m
    N = Nlist(k);

    % CFL rule per assignment (same for uniform and Chebyshev)
    dx_ref = 1/N;
    dt0    = CFL*dx_ref;         % |c|=1
    nsteps = ceil(Tfinal/dt0);
    dt     = Tfinal/nsteps;

    [sstar, tstar] = burgers_sstar(N, nu, dt, nsteps, gridType, nlType);

    relerr = abs(sstar - s_exact)/abs(s_exact);

    R.N(k)      = N;
    R.nsteps(k) = nsteps;
    R.dt(k)     = dt;
    R.sstar(k)  = sstar;
    R.tstar(k)  = tstar;
    R.relerr(k) = relerr;

    if k > 1
        R.ratio(k) = prev_err / relerr;
    end
    prev_err = relerr;
end

end

% ============================================================
function [sstar, tstar] = burgers_sstar(N, nu, dt, nsteps, gridType, nlType)
    [xb, A, D] = build_ops(N, nu, gridType);
    n = N-1;
    I = speye(n);

    % IC: sin(pi x), interior only
    u1 = sin(pi*xb(2:end-1));
    u2 = u1; u3 = u1;

    Nu1 = nonlinear(u1, D, nlType);
    Nu2 = Nu1; Nu3 = Nu1;

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
            Nu = -(D*(0.5*(u.^2)));       % -((u^2)/2)_x
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
        jx = i+1;
        dm = xb(jx)   - xb(jx-1);
        dp = xb(jx+1) - xb(jx);

        % A ~ (-nu u_xx) (SPD)
        A(i,i) = 2*nu/(dm*dp);
        if i>1
            A(i,i-1) = -2*nu/(dm*(dm+dp));
            D(i,i-1) = -1/(dm+dp);
        end
        if i<n
            A(i,i+1) = -2*nu/(dp*(dm+dp));
            D(i,i+1) =  1/(dm+dp);
        end
    end
end

% ============================================================
function write_latex_table(fname, R, captionStr)
    fid = fopen(fname,'w');
    if fid<0, error('Could not open %s for writing.', fname); end

    fprintf(fid, '%% Auto-generated by p1e_four_cases_tables.m\n');
    fprintf(fid, '\\begin{table}[t]\\centering\n');
    fprintf(fid, '\\caption{%s. Table entries: $N$, $n_{\\rm steps}$, $s^*$, relative error, and ratio $\\mathrm{err}_{N/2}/\\mathrm{err}_N$.}\\label{tab:p1e_%s}\n', ...
        captionStr, fname_to_label(fname));
    fprintf(fid, '\\begin{tabular}{r r r r r}\\hline\n');
    fprintf(fid, '$N$ & $n_{\\rm steps}$ & $s^*$ & rel.err. & ratio \\\\ \\hline\n');

    for k=1:numel(R.N)
        if isnan(R.ratio(k))
            fprintf(fid, '%d & %d & %.5g & %.3e & --- \\\\\n', R.N(k), R.nsteps(k), R.sstar(k), R.relerr(k));
        else
            fprintf(fid, '%d & %d & %.5g & %.3e & %.3f \\\\\n', R.N(k), R.nsteps(k), R.sstar(k), R.relerr(k), R.ratio(k));
        end
    end

    fprintf(fid, '\\hline\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end

function lab = fname_to_label(fname)
    lab = regexprep(fname, '\.tex$', '');
    lab = regexprep(lab, '[^a-zA-Z0-9]+', '_');
end