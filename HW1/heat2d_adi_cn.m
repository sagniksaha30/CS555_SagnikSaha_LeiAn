%
% HEAT EQUATION, u_t = nu \nabla^2 u, USING m-POINT FINITE DIFFERENCE AND ADI
%

hdr

%N=150;
m=N-1;

Lx=1; Ly=1;
T = 1.2; nstep = ceil(T/dt); dt=T/nstep;
nu = .05;

dx=Lx/N; xb = dx*[0:N]'; x=xb(2:end-1);
e = ones(m,1); A = spdiags([-e 2*e -e],-1:1, m,m);
Ax = nu*A/(dx*dx);
Ix = speye(m);


dy=Ly/N; yb = dy*[0:N]'; y=yb(2:end-1);
e = ones(m,1); A = spdiags([-e 2*e -e],-1:1, m,m);
Ay = nu*A/(dy*dy);
Iy = speye(m);

[X,Y]=ndgrid(x,y);
[Xb,Yb]=ndgrid(xb,yb);

%
%  Set up standard CN operators
%

%Don't need the full matrices anymore. Comment out, keep for reference.
%HL = kron(Iy,Ix) + (dt/2)*(kron(Iy,Ax)+kron(Ay,Ix)); 
%HR = kron(Iy,Ix) - (dt/2)*(kron(Iy,Ax)+kron(Ay,Ix));

%
%  Set up Exact Solution + RHS
%

kx = 1; ky = 3;

U0 = sin(kx*pi*X/Lx).*sin(ky*pi*Y/Ly);
U  = U0;

lamxy = -nu*( (kx*pi/Lx)^2 + (ky*pi/Ly)^2 );
thx   = pi*kx*dx/Lx;
thy   = pi*ky*dy/Ly;
lamx  = 2*(1-cos(thx))/(dx*dx);
lamy  = 2*(1-cos(thy))/(dy*dy);
lamxy = -nu*(lamx+lamy);              % Judge accuracy by discrete lambdas

%u = reshape(U,m*m,1); %Don't need to reshape into a vector here. Can work directly with U. Keeping this line for reference.

%x matrices. Note that the A's already have \nu multiplied to them so I don't need to include \nu here.
Tlx=Ix + dt/2 * Ax;
Tly= Iy + dt/2 * Ay;
Trx=Ix - dt/2 * Ax;
Try=Iy - dt/2 * Ay;

%Solving.
for istep=1:nstep; time=istep*dt;
  rhs = Trx*U*Try';
  temp = Tly\(rhs');
  U = Tlx\temp';
end;

%U = reshape(u,m,m); %Again, already a matrix, no need to reshape.
Uex = exp(lamxy*time)*U0;
Err = Uex-U;
eo  = e2;
e2  = norm(Err,"fro") / norm(Uex,"fro");
ratio = eo/e2;

format shorte;
disp([N dt nstep time e2 ratio])
%disp([e2 ratio])
% semilogy(t2k,e2k,'r-',lw,2,t2k,e2k,'k.',ms,11); drawnow; hold on;
