% problem2a.m
% Compressible Reynolds equation for the taper-flat slider (HW4, Part 2a)
%
% Place this file in the same folder as the starter-code helper functions:
%   hdr.m, box_elem.m, abqfem.m, profile_taper_flat.m, trigausspoints.m,
%   basis_deriv_12.m, basis_tri_12.m, boundedges.m, restriction.m
%
% This script solves the compressible Reynolds equation
%   -div( pa*h^3/(12*mu) grad p ) + (1/2) div(u h p) = -(patm/2) div(u h)
% with
%   pa = p + patm,   p = 0 on boundary,
% by fixed-point iteration.
%
% It reports:
%   h2, gamma, load F, center-of-pressure xp, and xp/L,
% and makes a mesh plot of the pressure distribution.

clear; clc; close all;
hdr;

%% -------------------- Parameters --------------------
Ex = 120;
Ey = 40;

rho = 1.225;          % kg/m^3      density of air
nu_air = 1.5e-5;      % m^2/s       kinematic viscosity
mu = rho*nu_air;      %             dynamic viscosity
U = 20;               % m/s         disk speed
L = 0.0041;           % m           slider length
W = 0.15625*L;        % m           slider width
T = 0.09375*L;        % m           taper length
Ta = pi/180;          % rad         taper angle (1 degree)

h1 = 3.70e-7;         % m           leading-edge height (without taper)
h2 = 2.50e-7;         % m           trailing-edge height
patm = 101325;        % N/m^2       atmospheric pressure

gamma = (h1-h2)/L;    % rad         small-angle pitch approximation

%% -------------------- Mesh and FEM setup --------------------
[x,y,t] = box_elem(Ex,Ey,L,W);
E  = 2*Ex*Ey;         % number of triangles
nv = 3;               % linear triangles

% element-local coordinate arrays (3 x E)
xL = x(t');
yL = y(t');

% element centroids for coefficient evaluation
xe = sum(xL,1)'/nv;
ye = sum(yL,1)'/nv; %#ok<NASGU>

[AL,BL,Q,t,areaL] = abqfem([x y],t); %#ok<ASGLU>
Bb = Q'*BL*Q;                % global consistent mass matrix
nb = size(Q,2);              % total number of global nodes

%% -------------------- Build advection matrix Cb --------------------
order = 1;                   % only linear elements supported here
nvq   = 3*order;
Nq    = 3;                   % triangle quadrature rule
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

CL = (DLr*ULr + DLs*ULs) * JL;
Cb = Q' * CL * Q;

%% -------------------- Dirichlet boundary restriction --------------------
% For Part 2(a), p = 0 on the full boundary.
boundary_nodes = boundedges([x,y],t);
R = restriction(nb,boundary_nodes);

%% -------------------- Compressible fixed-point iteration --------------------
% Start from p = 0, i.e. pa = patm.
Pb_old = zeros(nb,1);        % full-node pressure vector (relative pressure)
tol    = 1e-10;
maxit  = 30;
relchg_hist = zeros(maxit,1);

for it = 1:maxit

    % Absolute pressure at nodes
    pa_nodes = Pb_old + patm*ones(nb,1);

    % Elementwise absolute pressure: use nodal average on each triangle
    pa_elem = sum(pa_nodes(t'),1)'/3;

    % Bearing height coefficient on elements
    he = profile_taper_flat(xe,L,T,Ta,h1,h2);
    h3 = he.^3;

    % nu_coeff = pa * h^3 / (12*mu) on each element
    nu_coeff_elem = (pa_elem .* h3) / (12*mu);
    NuElem = spdiags(nu_coeff_elem,0,E,E);
    NuLoc  = kron(NuElem, speye(3));

    % Diffusion-like operator Abar(pa)
    Abar = Q' * (NuLoc * AL) * Q;

    % Solve: R(Abar - Cbar)R^T p = R Cbar p_atm
    A = R * (Abar - Cb) * R';
    rhs = R * (Cb * (patm*ones(nb,1)));

    P  = A \ rhs;            % reduced unknowns
    Pb = R' * P;             % full-node relative pressure

    relchg = norm(Pb - Pb_old) / max(norm(Pb),1);
    relchg_hist(it) = relchg;
    fprintf('it = %2d, relchg = %.3e, pmax = %.6e Pa\n', it, relchg, max(Pb));

    if relchg < tol
        relchg_hist = relchg_hist(1:it);
        break;
    end

    Pb_old = Pb;

    if it == maxit
        warning('Fixed-point iteration reached maxit before hitting tol.');
        relchg_hist = relchg_hist(1:it);
    end
end

%% -------------------- Load and center-of-pressure --------------------
% In FEM form:
%   F  = integral_Omega p dOmega  = 1^T B p
%   xp = (integral_Omega x p dOmega) / F = x^T B p / F
Bp = Bb * Pb;
F  = sum(Bp);
xp = (x' * Bp) / F;

%% -------------------- Report results --------------------
fprintf('\nCompressible taper-flat bearing results (Problem 2a)\n');
fprintf('---------------------------------------------------\n');
fprintf('h2       = %.12e m\n', h2);
fprintf('gamma    = %.12e rad\n', gamma);
fprintf('Load F   = %.12e N\n', F);
fprintf('xp       = %.12e m\n', xp);
fprintf('xp/L     = %.12f\n', xp/L);

%% -------------------- Plots --------------------
Pmax = max(abs(Pb));
if Pmax == 0
    scale = 1;
else
    scale = W / Pmax;        % make the pressure mesh roughly as tall as bearing width
end

figure;
trimesh(t,x,y,scale*Pb);
axis([0 L -W/1.9 W/1.9 0 W]);
axis equal;
xlabel('x',fs,20);
ylabel('y',fs,20);
zlabel('scaled pressure',fs,20);
title('Problem 2(a): Compressible taper-flat bearing',fs,15);

figure;
semilogy(1:length(relchg_hist), relchg_hist, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('iteration');
ylabel('relative change');
title('Fixed-point convergence history',fs,15);
