hdr;

Ex = 120;
Ey = 40;

rho = 1.225;         % kg/m^3      density of air
nu  = 1.5e-5;        % m^2/s       kinematic viscosity
mu  = rho*nu;        %             dynamic viscosity
U   = 20;            % m/s         speed of plate
L   = .0041;         % m           slider length
W   = .15625*L;      % m           slider width (not used in 1D)
T   = .09375*L;      % m           slider taper length
Ta  = 0;        % m           Taper angle (1-deg)
h1  = 3.70e-7;        % m           1/2-micron leading  edge gap
h2  = 2.50e-7;       % m           1/4-micron trailing edge gap
gamma = (h1-h2)/L;   %             pitch angle of slider << 1
patm= 101325;        %             Pressure in N/m^2
n_per_gram= 9.81e-3; % N/g         Weight of a 1 gram mass

[x,y,t] = box_elem(Ex,Ey,L,W);
E  = 2*Ex*Ey;

nv = 3;            %% Linear triangles for elements

xL=x(t'); yL=y(t');
xe=sum(xL,1)'/nv;   %% Use Centroids for Laplacian quadrature points
ye=sum(yL,1)'/nv;

[AL,BL,Q,t,areaL]=abqfem([x y],t);
Ab=Q'*AL*Q;
Bb=Q'*BL*Q;
nb=size(Q,2);

he = profile_taper_flat(xe,L,T,Ta,h1,h2); %% Evaluate on quad points
h3 = he.*he.*he;
nu = (1./(12*mu))*h3;  
nu = spdiags(nu,0,E,E);

Iv = speye(nv);
nu = kron(nu,Iv);
An = nu*AL;        %% AL, each element scaled by local nu


order=1;  % order=2 is not yet working
nv=3*order;

Nq=3; % 8 is Max allowable for trigausspoints
[z,w]=trigausspoints(Nq); rq=z(:,1);sq=z(:,2);  Bh=diag(w);

nq=length(rq); Dr=zeros(nq,nv); Ds=Dr; Jq=zeros(nq,nv);
for k=1:nv; 
 [Dr(:,k) Ds(:,k)]=basis_deriv_12(rq,sq,k,order);  % Differentiate to fine mesh
  Jq(:,k)         =basis_tri_12  (rq,sq,k,order);  % Interpolate to fine mesh
end;

Xr  = Dr*xL; Yr = Dr*yL; Xs  = Ds*xL; Ys = Ds*yL;
Jac = Xr.*Ys - Xs.*Yr; Jmin =min(min(Jac)); if Jmin <= 0; error('vanishing Jacobian'); end;
Rx  = Ys ./ Jac; Ry = -Xs ./ Jac; Sy  = Xr ./ Jac; Sx = -Yr ./ Jac;
Bq  = Bh*Jac;

hL  = profile_taper_flat(xL,L,T,Ta,h1,h2); %% Evaluate on quad points - this is the h(x,y) (gap b/w slider and disk)
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

boundary_nodes = boundedges([x,y],t);
%R = restriction(nb,boundary_nodes); %Fully Dirichlet.
R = restriction(nb, find(x==0 | x==L)); %von Neumann top and bottom.

A = R*(Q'*An*Q)*R';
B = R*Bb*R';
C = R*Cb*R';

rhs = R*(Cb*ones(nb,1));
P   = A \ rhs;
Pb  = R'*P;
Pmax=max(max(abs(Pb)));

%Compute the force.
F = ones(nb,1)' * Bb * Pb;
%Compute center-of-pressure.
xp = x(:);
Mx = xp' * Bb * Pb; 
Mx = Mx/F;

disp(['F = ', num2str(F)])
disp(['xp/L = ', num2str(Mx/L)])

%1a verification.
alpha = h1/h2;
lambda = 6*mu*U*L/(h2*h2);


%Analytic force (1D)
Ftild = ((lambda * L * W)/(1-alpha)^2)*(log(alpha) + 2*(1-alpha)/(1+alpha));

%disp(F)
%disp(Ftild)
%disp(Mx)

%Analytic pressure (1D) - function of x.
function p = ptild(xval,alpha,lambda,L)
H = alpha + (1-alpha)*xval/L;
p = ((alpha * lambda)/(1-alpha^2))*(1/H^2 - 1/alpha^2) - (lambda/(1-alpha))*(1/H - 1/alpha);
end

%Compute analytic pressure over the range of xvalues.
p_an = zeros(length(x),1);
for i = 1:length(x)
    p_an(i) = ptild(x(i), alpha, lambda, L);
end
mpe = max(abs(p_an - Pb));

%%%%Plotting stuff%%%%%%%
%axis([0 L -W/1.9 W/1.9 0 W]);
%axis equal;
%xlabel('x',fs,20);
%ylabel('y',fs,20);
%zlabel('W*p/p_{\max}',fs,20);
%title('Incompressible Case',fs,15);
%pause()

%trimesh(t,x,y,W*Pb/Pmax);
%xlabel('x',fs,20);
%ylabel('y',fs,20);
%zlabel('Scaled p',fs,20);
%pause()