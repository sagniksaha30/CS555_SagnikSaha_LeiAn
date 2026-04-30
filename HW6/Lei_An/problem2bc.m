clear; clc; close all;
hdr;

%% Common parameters
Ex = 120;
Ey = 40;

rho = 1.225;         % kg/m^3 density of air
nu_air = 1.5e-5;     % m^2/s kinematic viscosity
mu = rho*nu_air;     % dynamic viscosity
U = 20;              % m/s speed of plate
L = .0041;           % m slider length
W = .15625*L;        % m slider width
T = .09375*L;        % m slider taper length
Ta = pi/180;         % rad taper angle (1 deg)

h1 = 3.70e-7;        % m leading edge gap (without taper)
h2 = 2.50e-7;        % m trailing edge gap
patm = 101325;       % N/m^2 atmospheric pressure

gamma = (h1-h2)/L;   % small-angle pitch

%% Case 1a: incompressible taper-flat
case1a = solve_slider_case(Ex,Ey,L,W,T,Ta,h1,h2,U,mu,patm,false,false);

%% Case 1b: incompressible wedge (Ta=0, Neumann top/bottom)
case1b = solve_slider_case(Ex,Ey,L,W,T,0.0,h1,h2,U,mu,patm,false,true);

%% Case 2a: compressible taper-flat
case2a = solve_slider_case(Ex,Ey,L,W,T,Ta,h1,h2,U,mu,patm,true,false);

%% Print results
fprintf('\nComparison of the three cases\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%-28s %-16s %-16s %-16s\n','Case','F (N)','xp (m)','xp/L');
fprintf('---------------------------------------------------------------\n');
fprintf('%-28s %-16.12e %-16.12e %-16.12f\n','1a: incompressible taper-flat',case1a.F,case1a.xp,case1a.xp/L);
fprintf('%-28s %-16.12e %-16.12e %-16.12f\n','1b: incompressible wedge',case1b.F,case1b.xp,case1b.xp/L);
fprintf('%-28s %-16.12e %-16.12e %-16.12f\n','2a: compressible taper-flat',case2a.F,case2a.xp,case2a.xp/L);

fprintf('\nPeak relative pressures\n');
fprintf('1a: pmax = %.12e Pa\n', max(case1a.Pb));
fprintf('1b: pmax = %.12e Pa\n', max(case1b.Pb));
fprintf('2a: pmax = %.12e Pa\n', max(case2a.Pb));

%% Side-by-side pressure plots
figure;
subplot(1,3,1);
scale1 = W/max(abs(case1a.Pb));
trimesh(case1a.t,case1a.x,case1a.y,scale1*case1a.Pb);
axis([0 L -W/1.9 W/1.9 0 W]);
axis equal;
view(35,25);
title('1a: incompressible taper-flat',fs,13);
xlabel('x',fs,16); ylabel('y',fs,16); zlabel('scaled p',fs,16);

subplot(1,3,2);
scale2 = W/max(abs(case1b.Pb));
trimesh(case1b.t,case1b.x,case1b.y,scale2*case1b.Pb);
axis([0 L -W/1.9 W/1.9 0 W]);
axis equal;
view(35,25);
title('1b: incompressible wedge',fs,13);
xlabel('x',fs,16); ylabel('y',fs,16); zlabel('scaled p',fs,16);

subplot(1,3,3);
scale3 = W/max(abs(case2a.Pb));
trimesh(case2a.t,case2a.x,case2a.y,scale3*case2a.Pb);
axis([0 L -W/1.9 W/1.9 0 W]);
axis equal;
view(35,25);
title('2a: compressible taper-flat',fs,13);
xlabel('x',fs,16); ylabel('y',fs,16); zlabel('scaled p',fs,16);

sgtitle('Problem 2(b): pressure distributions for cases 1a, 1b, and 2a',fs,16);

%% Optional table display in MATLAB
Tcomp = table( ...
    [case1a.F; case1b.F; case2a.F], ...
    [case1a.xp; case1b.xp; case2a.xp], ...
    [case1a.xp/L; case1b.xp/L; case2a.xp/L], ...
    [max(case1a.Pb); max(case1b.Pb); max(case2a.Pb)], ...
    'RowNames', {'1a_incompressible_taper_flat','1b_incompressible_wedge','2a_compressible_taper_flat'}, ...
    'VariableNames', {'F_N','xp_m','xp_over_L','pmax_Pa'} );
disp(Tcomp);

%% ------------------------------------------------------------
function out = solve_slider_case(Ex,Ey,L,W,T,Ta,h1,h2,U,mu,patm,isCompressible,isWedge)

[x,y,t] = box_elem(Ex,Ey,L,W);
E  = 2*Ex*Ey;
nv = 3;

xL = x(t');
yL = y(t');
xe = sum(xL,1)'/nv;

[AL,BL,Q,t,areaL] = abqfem([x y],t); %#ok<ASGLU>
Bb = Q'*BL*Q;
nb = size(Q,2);

%% Build Cb
order = 1;
nvq   = 3*order;
Nq    = 3;
[z,w] = trigausspoints(Nq);
rq = z(:,1);
sq = z(:,2);
Bh = diag(w);

nq = length(rq);
Dr = zeros(nq,nvq);
Ds = zeros(nq,nvq);
Jq = zeros(nq,nvq);
for k = 1:nvq
    [Dr(:,k), Ds(:,k)] = basis_deriv_12(rq,sq,k,order);
    Jq(:,k)            = basis_tri_12(rq,sq,k,order);
end

Xr  = Dr*xL;
Yr  = Dr*yL;
Xs  = Ds*xL;
Ys  = Ds*yL;
Jac = Xr.*Ys - Xs.*Yr;
if min(Jac(:)) <= 0
    error('Vanishing or negative Jacobian detected.');
end

Rx =  Ys ./ Jac;
Sx = -Yr ./ Jac;
Bq =  Bh * Jac;

hL = profile_taper_flat(xL,L,T,Ta,h1,h2);
hq = Jq*hL;
Uh = (U/2) * (Bq .* hq);

Ie  = speye(E);
DLr = kron(Ie,Dr');
DLs = kron(Ie,Ds');
ULr = spdiags(reshape(Rx.*Uh,nq*E,1),0,nq*E,nq*E);
ULs = spdiags(reshape(Sx.*Uh,nq*E,1),0,nq*E,nq*E);
JL  = kron(Ie,Jq);
CL  = (DLr*ULr + DLs*ULs) * JL;
Cb  = Q' * CL * Q;

%% Restriction matrix
if isWedge
    boundary_nodes = find(abs(x) < 1e-14 | abs(x-L) < 1e-14);
else
    boundary_nodes = boundedges([x,y],t);
end
R = restriction(nb,boundary_nodes);

%% Solve
if ~isCompressible
    he = profile_taper_flat(xe,L,T,Ta,h1,h2);
    h3 = he.^3;
    nu_elem = (1./(12*mu))*h3;
    NuElem = spdiags(nu_elem,0,E,E);
    NuLoc  = kron(NuElem, speye(3));
    An     = NuLoc*AL;

    A   = R*(Q'*An*Q)*R';
    rhs = R*(Cb*ones(nb,1));
    P   = A \ rhs;
    Pb  = R'*P;
else
    Pb_old = zeros(nb,1);
    tol = 1e-10;
    maxit = 30;

    for it = 1:maxit
        pa_nodes = Pb_old + patm*ones(nb,1);
        pa_elem  = sum(pa_nodes(t'),1)'/3;

        he = profile_taper_flat(xe,L,T,Ta,h1,h2);
        h3 = he.^3;
        nu_elem = (pa_elem .* h3)/(12*mu);
        NuElem = spdiags(nu_elem,0,E,E);
        NuLoc  = kron(NuElem, speye(3));
        Abar   = Q'*(NuLoc*AL)*Q;

        A   = R*(Abar - Cb)*R';
        rhs = R*(Cb*(patm*ones(nb,1)));

        P  = A \ rhs;
        Pb = R'*P;

        relchg = norm(Pb - Pb_old)/max(norm(Pb),1);
        if relchg < tol
            break;
        end
        Pb_old = Pb;
    end
end

%% Load and center-of-pressure
Bp = Bb*Pb;
F  = sum(Bp);
xp = (x' * Bp)/F;

out.x  = x;
out.y  = y;
out.t  = t;
out.Pb = Pb;
out.F  = F;
out.xp = xp;
end
