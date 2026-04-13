%Define the problem parameters.
L = 1;
c = 1;
nu = 0.001;

V = 1; %No. of vertices = 1 due to linear interpolation.
E = 100; %No. of elements. 
N = E; %No. of global nodes. Equal to E due to indexing offset (0,1,...,N as opposed to 1,2,...,E).


%%%%%%%Define the Q matrix%%%%%%%%%%.
%We have to do some index gymnastics since v is 0 indexed, but everything else (including Octave itself) is 1 indexed.
V_plus_1 = V+1;
% We need to define the coordinates for every '1' in the Q matrix
% There are (V+1)*E total ones (2 per element)
rows = zeros((V_plus_1)*E, 1);
cols = zeros((V_plus_1)*E, 1);
vals = ones((V_plus_1)*E, 1);

counter = 1;
for e = 1:E
    for v = 0:1
        % Map to Octave 1-based row.
        l = v+2*e-1; % = v+1 + 2(e-1). Shift e to be 0 based, do the calculation, and revert to 1 based.

        % Map to Octave 1-based column
        i_octave = e + v; 
        
        rows(counter) = l;
        cols(counter) = i_octave;
        counter = counter + 1;
    end
end

Q = sparse(rows, cols, vals, 2*E, E+1); % Size: (Rows in local vector) x (Rows in global vector)

%%%%%%%Define the R matrix%%%%%%%%%%.
R = [zeros(N, 1), speye(N)]; %This is for Dirichlet on left end, Neumann on right end.

%%%%%%%Define the x values (global nodes) and the Le's. %%%%%%
%Uniform spacing.
xvals = linspace(0, 1, E+1);
Le = xvals(2)-xvals(1);

%%%%%%%Define A^e, B^e, C^e%%%%%%
Ae = zeros(2,2);
Be = zeros(2,2);
Ce = zeros(2,2);

Ae(1,1)=1/Le;
Ae(1,2)=-1/Le;
Ae(2,1)=-1/Le;
Ae(2,2)=1/Le;

Be(1,1)=Le/3;
Be(1,2)=Le/6;
Be(2,1)=Le/6;
Be(2,2)=Le/3;

Ce(1,1)=-1/2;
Ce(1,2)=1/2;
Ce(2,1)=-1/2;
Ce(2,2)=1/2;

%%%%%%%%%Define A_L, B_L, C_L%%%%%%%%%%
AL = kron(eye(E), Ae);
BL = kron(eye(E), Be);
CL = kron(eye(E), Ce);

%%%%%%%%%%Define Abar, Bbar, Cbar.%%%%%%
Abar = Q' * AL * Q;
Bbar = Q' * BL * Q;
Cbar = Q' * CL * Q;

%%%%%%%%%% Construct A, B, C %%%%%%
A = R * Abar * R';
B = R * Bbar * R';
C = R * Cbar * R';


%%%%%%%Implement BDF3/EXT3 from here%%%%%%%%
T = 1.5;
dt = 0.001;
nstep = T/dt;

%Make lists of all the BDFk coefficients. Store as column vectors.
%Careful: bdf(k) corresponds to \beta_{k-1}.
bdf1  = (([ 1.  -1.  0.  0. 0.])/1.)';
bdf2  = (([ 3.  -4.  1.  0. 0.])/2.)';
bdf3  = (([11. -18.  9. -2. 0.])/6.)';

%Make lists of all the EXTk coefficients. Store as column vectors.
ext1  = (([ 1.  0.  0.  0. 0.]))';
ext2  = (([ 2.  -1.  0.  0. 0.]))';
ext3  = (([3. -3. 1. 0. 0.]))';

%Vectors that will be used to store the correct BDFk/EXTk coefficients according 
%to the situation.
bdf = zeros(size(bdf1));
ext = zeros(size(ext1));

x  = xvals(2:end)';
u0 = zeros(size(x)); %Since IC = 0.
u=u0;

%Temp arrays to store previous results.
ulm1=zeros(size(u));
ulm2=zeros(size(u));
ulm3=zeros(size(u));

%For the first time step, u0 is the u_{l-1}.
ulm1 = u0;

%Time evolution

for k=1:nstep;  time=k*dt;
ub = zeros(N+1,1);
ub(1,1) = sin(pi * time);

ubd = zeros(N+1,1);
ubd(1,1) = pi*cos(pi * time);

  if(k==1) %BDF1/EXT1
    bdf = bdf1;
    ext = ext1;

    %Left side operator.
    H = B*bdf(1) + nu * A * dt;
    %Right side
    rhsbeta = -B*(bdf(2)*ulm1); %BDF1
    rhsalpha = -dt*C* (ext(1)*(ulm1)); %EXT1
    rhsleftbound = -R * dt* ((nu * Abar + Cbar)*ub + Bbar *ubd); %From Dirichlet BC at x = 0.
    
    %Combine all terms to get RHS.
    rhs = rhsbeta + rhsalpha + rhsleftbound;

    %Solve to get current 'u'.
    u   = H\(rhs);

    %Update past results.
    ulm2=ulm1; 
    ulm1=u;
  
  elseif(k==2) %BDF2/EXT2
    bdf = bdf2;
    ext = ext2;

    %Left side operator.
    H = B*bdf(1) + nu * A * dt;
    %Right side
    rhsbeta = -B*(bdf(2)*ulm1 + bdf(3)*ulm2); %BDF2
    rhsalpha = -dt*C* (ext(1)*(ulm1) + ext(2)*(ulm2)); %EXT2
    rhsleftbound = -R * dt* ((nu * Abar + Cbar)*ub + Bbar *ubd); %From Dirichlet BC at x = 0.
    
    %Combine all terms to get RHS.
    rhs = rhsbeta + rhsalpha + rhsleftbound;

    %Solve to get current 'u'.
    u   = H\(rhs);

    %Update past results.
    ulm3 = ulm2;
    ulm2 = ulm1; 
    ulm1=u;
  
  else %BDF3/EXT3
    bdf = bdf3;
    ext = ext3;
    %Left side operator.
    H = B*bdf(1) + nu * A * dt;
    %Right side
    rhsbeta = -B*(bdf(2)*ulm1 + bdf(3)*ulm2 + bdf(4)*ulm3); %BDF3
    rhsalpha = -dt*C* (ext(1)*(ulm1) + ext(2)*(ulm2) + ext(3)*(ulm3)); %EXT3
    rhsleftbound = -R * dt* ((nu * Abar + Cbar)*ub + Bbar *ubd); %From Dirichlet BC at x = 0.
    
    %Combine all terms to get RHS.
    rhs = rhsbeta + rhsalpha + rhsleftbound;

    %Solve to get current 'u'.
    u   = H\(rhs);

    %Update past results.
    ulm3 = ulm2;
    ulm2 = ulm1; 
    ulm1=u;
    
    end
  end

%%%%%%%%%Plotting the numerical solution%%%%%%%%%%%
u_full = zeros(E+1,1);       % full solution including boundaries
u_full(1) = sin(pi*time);     % left Dirichlet
u_full(2:end) = u;          % interior + right Neumann    
plot(xvals, u_full, '-o','markersize',2)
%xlabel('x'); ylabel('u(x)');
%title('Steady-state FEM solution');
%grid on
pause

%%%%%%%Pointwise error%%%%%%%%
%Define the exact solution
%function y=exact_steadysol(x,nu)
 %   y=x - (exp(-1/nu) - exp((x-1)/nu))/(exp(-1/nu) - 1);
%end

%Evaluate the exact solution over the range of xvals. Note that we exclude the boundary terms x = 0 and 1.
%xvals = xvals'; %Convert xvals to column vector since we originally used linspace.
%exactsol = exact_steadysol(xvals(2:end-1), nu);

%Compute the pointwise error (again, ignoring boundary terms).
%err_pntwise = abs((exactsol - res)./exactsol); 

%Display the maximum pointwise error.
%disp(max(err_pntwise))

%Plot the pointwise error as a function of 'x'.
%plot(xvals(2:end-1),err_pntwise,'-o','markersize',2)
%pause

