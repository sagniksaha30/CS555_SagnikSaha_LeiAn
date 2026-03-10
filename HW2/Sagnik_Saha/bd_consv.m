%
%L = 1;
n  = N-1; % Number of points. Homogeneous Dirichlet boundary conditions.
nu = 0.003183;
%nu = 0.01;
ad_var; %Defining the diffusion and first derivative operator.

%dt = 0.0005; %Keeping dt constant as final T changed (in the driver script).
nstep = T/dt;
%dt = T/nstep;
%iostep = floor(nstep/10); For plotting purposes

%Make lists of all the BDFk coefficients. Store as column vectors.
%Careful: bdf(k) corresponds to \beta_{k-1}.
bdf1  = (([ 1.  -1.  0.  0. 0.])/1.)';
bdf2  = (([ 3.  -4.  1.  0. 0.])/2.)';
bdf3  = (([11. -18.  9. -2. 0.])/6.)';

%Make lists of all the EXTk coefficients. Store as column vectors.
ext1  = (([ 1.  0.  0.  0. 0.]))';
ext2  = (([ 2.  -1.  0.  0. 0.]))';
ext3  = (([3. -3. 1. 0. 0.]))';

%Create the nxn identity matrix.
I  = speye(n);

%Vectors that will be used to store the correct BDFk/EXTk coefficients according 
%to the situation.
bdf = zeros(size(bdf1));
ext = zeros(size(ext1));

%Specify the initial conditon.
sx=sin(pi*x);
u0 = sx;
u  = u0;

%Temp arrays to store previous results.
ulm1=zeros(size(u));
ulm2=zeros(size(u));
ulm3=zeros(size(u));

%For the first time step, u0 is the u_{l-1}.
ulm1 = u0;

%Time evolution
for k=1:nstep;  time=k*dt;
  if(k==1) %BDF1/EXT1
    bdf = bdf1;
    ext = ext1;

    %Left side operator.
    H = I*bdf(1) + A*dt;

    %Right side
    rhsbeta = -(bdf(2)*ulm1);
    rhsalpha = ext(1)*(0.5* D*(ulm1.^2));
    rhs = rhsbeta - dt*rhsalpha;

    %Solve to get current 'u'.
    u   = H\(rhs);

    %Update past results.
    ulm2=ulm1; 
    ulm1=u;
  
  elseif(k==2) %BDF2/EXT2
    bdf = bdf2;
    ext = ext2;

    %Left side operator.
    H = I*bdf(1) + A*dt;

    %Right side
    rhsbeta = -(bdf(2)*ulm1 + bdf(3)*ulm2);
    rhsalpha = ext(1)*(0.5* D*(ulm1.^2)) + ext(2)*(0.5* D*(ulm2.^2));
    rhs = rhsbeta - dt*rhsalpha;

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
    H = I*bdf(1) + A*dt;

    %Right side
    rhsbeta = -(bdf(2)*ulm1 + bdf(3)*ulm2 + bdf(4)*ulm3);
    rhsalpha = ext(1)*(0.5* D*(ulm1.^2)) + ext(2)*(0.5* D*(ulm2.^2)) + ext(3)*(0.5* D*(ulm3.^2));
    rhs = rhsbeta - dt*rhsalpha;

    %Solve to get current 'u'.
    u   = H\(rhs);

    %Update past results.
    ulm3 = ulm2;
    ulm2 = ulm1; 
    ulm1=u;
    
    end
  end

% PLOT AXIS and Initial Condition:
uinitial=[0;u0;0];
ufinal=[0;u;0];

%hold off;

%plot(xb,uinitial,'k-'); 
%hold on;
%plot(xb,ufinal,color); hold on;
%axis([0 1 -.5 1.3]); axis square;
%legend_entries{k} = ['T = ', num2str(T)];

%xlabel('-- x --','FontSize',14); ylabel('u(x,t)','FontSize',14);
%pause;
%for k=1:nstep;  time=k*dt; 

 % if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',2); 
  % title(['N = ' int2str(N) ',  Step = ' int2str([k])],'fontsize',15); pause(.1); end;

  %u   = H\u;              name='EB';  % Euler Backward

  %uex = exp(time*lam_ex)*sx;
  %err = uex-u;
  %emx = max(abs(err))/max(abs(uex));
  %el2 = sqrt( (err'*err) / (uex'*uex) );
  %elast = enew;

%end;

%enew  = el2; ratio = elast/enew;

%format longe
%disp([name '    ' num2str([N nstep el2 ratio ],' %6d  %6d  %e  %e')])

