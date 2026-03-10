function p1c_sstar_vsN()
% Problem 1c (CS555 HW2):
%  - Burgers: u_t + u u_x = nu u_xx on [0,1], u(0,t)=u(1,t)=0, u(x,0)=sin(pi x)
%  - Spatial grids: (i) uniform, (ii) Chebyshev x_j=(1-cos(pi j/N))/2
%  - Operators: 3-pt stencils for diffusion A (implicit) and first-derivative D (explicit)
%  - Time: IMEX BDF3/EXT3, bootstrapped by BDF1/EXT1 then BDF2/EXT2
%  - Metric: s(t) = max_x |u_x|,  s* = max_t s(t)
%  - Study: fix CFL = |c| dt / dx = 0.1 with dx=1/N and |c|=1  => dt ~ 0.1/N
%  - Compare: relative error in s* vs analytical value from Table 3 in [1].

  % ----- user settings -----
  nu     = 1/(100*pi);     % match Basdevant et al. paper setting (recommended for Table 3)
  Tfinal = 2.0;
  CFL    = 0.1;
  Nlist  = 2.^(5:10);      % 32..1024 (edit if you want)

  % Analytical s* from Table 3 in Basdevant et al. (1986): "Analytical" row
  s_star_exact = 152.00516;  % see Table 3 in [1]

  % Run both grids
  results_uniform = run_grid_family('uniform',   Nlist, nu, Tfinal, CFL, s_star_exact);
  results_cheb    = run_grid_family('chebyshev', Nlist, nu, Tfinal, CFL, s_star_exact);

  % Print to screen
  fprintf('\n=== Problem 1c: relative error in s* vs N (nu=%.6g, CFL=%.3g, T=%.2f) ===\n', nu, CFL, Tfinal);
  print_table(results_uniform, 'Uniform grid');
  print_table(results_cheb,    'Chebyshev grid');

  % Write LaTeX tables
  write_latex_table('p1c_table_uniform.tex', results_uniform, 'Uniform spacing');
  write_latex_table('p1c_table_cheb.tex',    results_cheb,    'Chebyshev spacing');

  fprintf('\nWrote LaTeX tables:\n  p1c_table_uniform.tex\n  p1c_table_cheb.tex\n');
end

% ============================================================
function results = run_grid_family(gridType, Nlist, nu, Tfinal, CFL, s_star_exact)
  m = numel(Nlist);
  results = struct();
  results.gridType = gridType;
  results.N        = zeros(m,1);
  results.dt       = zeros(m,1);
  results.nsteps   = zeros(m,1);
  results.sstar    = zeros(m,1);
  results.tstar    = zeros(m,1);
  results.relerr   = zeros(m,1);
  results.ratio    = nan(m,1);

  prev_err = NaN;

  for k = 1:m
    N = Nlist(k);

    % CFL rule requested in HW: use dx=1/N and |c|=1 for BOTH grids
    dx_ref = 1.0/N;
    dt0    = CFL*dx_ref;

    % Choose integer steps; keep dt <= dt0 (so CFL not exceeded)
    nsteps = ceil(Tfinal/dt0);
    dt     = Tfinal/nsteps;

    [sstar, tstar] = burgers_imex_sstar(N, nu, dt, nsteps, gridType);

    relerr = abs(sstar - s_star_exact)/abs(s_star_exact);

    results.N(k)      = N;
    results.dt(k)     = dt;
    results.nsteps(k) = nsteps;
    results.sstar(k)  = sstar;
    results.tstar(k)  = tstar;
    results.relerr(k) = relerr;

    if k > 1
      results.ratio(k) = prev_err / relerr;
    end
    prev_err = relerr;
  end
end

% ============================================================
function [sstar, tstar] = burgers_imex_sstar(N, nu, dt, nsteps, gridType)
  % Build grid and operators on interior unknowns
  [xb, A, D] = build_ops(N, nu, gridType); %#ok<ASGLU>
  n = N-1;                 % interior DOFs
  I = speye(n);

  % Initial condition u0 = sin(pi x), with homogeneous Dirichlet boundaries
  x_in = xb(2:end-1);
  u1   = sin(pi*x_in);     % u^n
  u2   = zeros(n,1);       % u^{n-1} (dummy init; will be overwritten after step 1)
  u3   = zeros(n,1);       % u^{n-2}

  % Nonlinear history N(u) = -u .* (D u) (convective form)
  Nu1 = nonlinear_term(u1, D);
  Nu2 = Nu1*0;
  Nu3 = Nu1*0;

  % Track s(t)
  tvals = (0:nsteps)*dt;
  svals = zeros(nsteps+1,1);
  ux    = D*u1;
  svals(1) = max(abs(ux));

  % Time marching
  for step = 1:nsteps
    k = min(step, 3);              % bootstrap: 1, then 2, then 3 forever

    [b0,b1,b2,b3,a1,a2,a3] = bdf_ext_coeffs(k);

    % Implicit matrix (diffusion): (b0 I + dt A) u^{n+1} = ...
    H = b0*I + dt*A;

    % EXTrapolated nonlinear term
    N_ext = a1*Nu1 + a2*Nu2 + a3*Nu3;

    rhs = -(b1*u1 + b2*u2 + b3*u3) + dt*N_ext;

    unew = H \ rhs;

    % update histories
    u3 = u2;  u2 = u1;  u1 = unew;
    Nu3 = Nu2; Nu2 = Nu1; Nu1 = nonlinear_term(u1, D);

    ux = D*u1;
    svals(step+1) = max(abs(ux));
  end

  % s* and t*
  [sstar, idx] = max(svals);
  tstar = tvals(idx);
end

% ============================================================
function Nu = nonlinear_term(u, D)
  ux = D*u;
  Nu = -(u .* ux);     % convective form: -u u_x
end

% ============================================================
function [b0,b1,b2,b3,a1,a2,a3] = bdf_ext_coeffs(k)
  % Matches the HW2 snippet coefficients
  if k==1
    b0=1;    b1=-1;  b2=0;    b3=0;
    a1=1;    a2=0;   a3=0;          % EXT1
  elseif k==2
    b0=1.5;  b1=-2;  b2=0.5;  b3=0;
    a1=2;    a2=-1;  a3=0;          % EXT2
  else
    b0=11/6; b1=-3;  b2=1.5;  b3=-1/3;
    a1=3;    a2=-3;  a3=1;          % EXT3
  end
end

% ============================================================
function [xb, A, D] = build_ops(N, nu, gridType)
  % Grid includes boundaries: xb(1)=0, xb(end)=1
  j = (0:N)';
  if strcmpi(gridType,'uniform')
    xb = j / N;
  elseif strcmpi(gridType,'chebyshev')
    xb = (1 - cos(pi*j/N))/2;
  else
    error('Unknown gridType="%s". Use "uniform" or "chebyshev".', gridType);
  end

  n = N-1; % interior
  A = spalloc(n,n,3*n);
  D = spalloc(n,n,2*n);

  % Build A approximating (-nu u_xx) with variable-spacing 3-pt stencil (SPD)
  % and D approximating u_x by the provided centered formula: (w_{j+1}-w_{j-1})/(x_{j+1}-x_{j-1})
  for i = 1:n
    jx = i+1; % maps interior index i -> grid index jx in xb (since xb(1) is boundary)

    dm = xb(jx)   - xb(jx-1);
    dp = xb(jx+1) - xb(jx);

    % Diffusion: -nu*u_xx  (SPD tri-diagonal)
    A(i,i) = 2*nu/(dm*dp);
    if i > 1
      A(i,i-1) = -2*nu/(dm*(dm+dp));
      D(i,i-1) = -1/(dm+dp);
    end
    if i < n
      A(i,i+1) = -2*nu/(dp*(dm+dp));
      D(i,i+1) =  1/(dm+dp);
    end
  end
end

% ============================================================
function print_table(R, titleStr)
  fprintf('\n--- %s (%s) ---\n', titleStr, R.gridType);
  fprintf('%8s %10s %10s %14s %12s %10s\n', 'N','dt','nsteps','s*','rel.err','ratio');
  for k=1:numel(R.N)
    if isnan(R.ratio(k))
      fprintf('%8d %10.3e %10d %14.8f %12.3e %10s\n', R.N(k), R.dt(k), R.nsteps(k), R.sstar(k), R.relerr(k), '---');
    else
      fprintf('%8d %10.3e %10d %14.8f %12.3e %10.3f\n', R.N(k), R.dt(k), R.nsteps(k), R.sstar(k), R.relerr(k), R.ratio(k));
    end
  end
end

% ============================================================
function write_latex_table(fname, R, captionStr)
  fid = fopen(fname,'w');
  if fid<0, error('Could not open %s for writing.', fname); end

  fprintf(fid, '%% Auto-generated by p1c_sstar_vsN.m\n');
  fprintf(fid, '\\begin{table}[t]\\centering\n');
  fprintf(fid, '\\caption{%s: relative error in $s^*$ vs $N$ (IMEX BDF3/EXT3, CFL=0.1).}\\label{tab:p1c_%s}\n', ...
          captionStr, lower(R.gridType));
  fprintf(fid, '\\begin{tabular}{r r r r r}\\hline\n');
  fprintf(fid, '$N$ & $\\Delta t$ & $n_{\\rm steps}$ & $s^*$ & rel.err. \\\\ \\hline\n');

  for k=1:numel(R.N)
    fprintf(fid, '%d & %.3e & %d & %.5g & %.3e \\\\\n', ...
            R.N(k), R.dt(k), R.nsteps(k), R.sstar(k), R.relerr(k));
  end

  fprintf(fid, '\\hline\\end{tabular}\n');
  fprintf(fid, '\\end{table}\n');
  fclose(fid);
end