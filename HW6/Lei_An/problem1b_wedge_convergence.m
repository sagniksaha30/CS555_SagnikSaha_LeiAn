clear; clc; close all;
hdr;

% Suppress the mesh figure that box_elem draws each time
set(0,'DefaultFigureVisible','off');

%% Fixed parameters from the assignment / starter code
rho = 1.225;         % kg/m^3
nu  = 1.5e-5;        % m^2/s
mu  = rho*nu;        % dynamic viscosity
U   = 20;            % m/s
L   = 0.0041;        % m
W   = 0.15625*L;     % m
T   = 0.09375*L;     % m
Ta  = 0.0;           % wedge test: no leading taper

h1  = 3.70e-7;       % m
h2  = 2.50e-7;       % m
Ey  = 40;            % required by the assignment

%% Analytical parameters
alpha  = h1/h2;
Lambda = 6*mu*U*L/h2^2;

%% Sweep in Ex
Ex_list = [10 20 40 80 160];
err_inf = zeros(size(Ex_list));

for k = 1:length(Ex_list)

    Ex = Ex_list(k);

    %% --- starter-code structure ---
    [x,y,t] = box_elem(Ex,Ey,L,W);
    E  = 2*Ex*Ey;

    nv = 3;
    xL = x(t'); 
    yL = y(t');
    xe = sum(xL,1)'/nv;
    ye = sum(yL,1)'/nv; %#ok<NASGU>

    [AL,BL,Q,t,areaL] = abqfem([x y],t); %#ok<ASGLU>
    Ab = Q'*AL*Q; %#ok<NASGU>
    Bb = Q'*BL*Q; %#ok<NASGU>
    nb = size(Q,2);

    he = profile_taper_flat(xe,L,T,Ta,h1,h2);
    h3 = he.*he.*he;
    nu_loc = (1./(12*mu))*h3;
    nu_loc = spdiags(nu_loc,0,E,E);

    Iv = speye(nv);
    nu_loc = kron(nu_loc,Iv);
    An = nu_loc*AL;

    order = 1;
    nvq   = 3*order;
    Nq    = 3;
    [z,w] = trigausspoints(Nq);
    rq = z(:,1); 
    sq = z(:,2);  
    Bh = diag(w);

    nq = length(rq); 
    Dr = zeros(nq,nvq); 
    Ds = Dr; 
    Jq = zeros(nq,nvq);

    for j = 1:nvq
        [Dr(:,j), Ds(:,j)] = basis_deriv_12(rq,sq,j,order);
        Jq(:,j)            = basis_tri_12(rq,sq,j,order);
    end

    Xr  = Dr*xL; 
    Yr  = Dr*yL; 
    Xs  = Ds*xL; 
    Ys  = Ds*yL;
    Jac = Xr.*Ys - Xs.*Yr;

    Rx  = Ys ./ Jac; 
    Ry  = -Xs ./ Jac; %#ok<NASGU>
    Sy  = Xr ./ Jac; %#ok<NASGU>
    Sx  = -Yr ./ Jac;
    Bq  = Bh*Jac;

    hL  = profile_taper_flat(xL,L,T,Ta,h1,h2);
    hq  = Jq*hL;
    Uh  = (U/2)*(Bq.*hq);

    Ie  = speye(E);
    DLr = kron(Ie,Dr');
    DLs = kron(Ie,Ds');
    ULr = spdiags(reshape(Rx.*Uh,nq*E,1),0,nq*E,nq*E);
    ULs = spdiags(reshape(Sx.*Uh,nq*E,1),0,nq*E,nq*E);
    JL  = kron(Ie,Jq);

    CL  = (DLr*ULr + DLs*ULs)*JL;
    Cb  = Q'*CL*Q;

    %% Dirichlet only at x=0 and x=L
    boundary_nodes = find(abs(x) < 1e-14 | abs(x-L) < 1e-14);

    R = restriction(nb,boundary_nodes);
    A = R*(Q'*An*Q)*R';
    rhs = R*(Cb*ones(nb,1));

    P  = A \ rhs;
    Pb = R'*P;

    %% Exact 1D wedge solution evaluated at node x-coordinates
    H      = alpha + (1-alpha)*(x/L);
    p_exact = alpha*Lambda/(1-alpha^2) .* (1./(H.^2) - 1/alpha^2) ...
            - Lambda/(1-alpha)       .* (1./H      - 1/alpha);

    %% Max nodal pointwise error
    err_inf(k) = max(abs(Pb - p_exact));

    fprintf('Ex = %4d, max nodal error = %.12e\n', Ex, err_inf(k));
end

%% Fit slope on log-log scale
pp    = polyfit(log(Ex_list), log(err_inf), 1);
slope = pp(1);

fprintf('\nEstimated slope on log-log plot = %.6f\n', slope);
fprintf('Estimated convergence order      = %.6f\n', -slope);

%% Plot convergence
set(0,'DefaultFigureVisible','on');

figure;
loglog(Ex_list, err_inf, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
loglog(Ex_list, exp(pp(2))*Ex_list.^pp(1), '--', 'LineWidth', 1.5);
grid on;
xlabel('E_x');
ylabel('max pointwise error in p(x)');
title('Problem 1(b): 1D wedge-bearing verification');
legend('FEM error', sprintf('fit slope = %.3f', slope), 'Location', 'southwest');