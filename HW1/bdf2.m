
%
% HEAT EQUATION, u_t = nu u_xx, USING n-POINT FINITE DIFFERENCE 
%

n  = N-1;          % Number of (spatial) points
nu = .20;
Lx = 1;

T  = 1.0;          % Final time
dt = T/nstep;
iostep = floor(nstep/10); %floor rounds to lower integer.

L  = Lx; a_var;  %% Set A based on variable (here, Chebyshev) spacing

I  = speye(n);
G  = 3*I + 2*dt*nu*A;   %BDF2 
H  = I + dt*nu*A;   %EB (for bootstrap)   

t1 = pi/N;
sx = sin(pi*x);   lam_ex = -nu*pi*pi; %Eigenvalue of the derivative matrix.

u0 = (pi/4)+0*x;  % Constant IC
%u0 = sx;          % Smooth IC
u  = u0;

%Temp arrays to store previous results.
ukm1=zeros(size(u));
ukm2=zeros(size(u));

% PLOT AXIS and Initial Condition:
ub=[0;u0;0]; hold off; 
plot(xb, ub, 'k-');     % initial condition
hold on
plot([0 1], [0 0], 'k-');  % horizontal line
axis([0 1 -.5 1.3]); axis square;
xlabel('-- x --','FontSize',14); ylabel('u(x,t)','FontSize',14);

for k=1:nstep;  time=k*dt; 
    if (k==1) %Bootstrap with EB.
        u   = H\(u); %Obtain u_1 = H^{-1}u_0.
        ukm2=u0;%u0 is u_{k-2} for k=2.
        ukm1=u;%u0 is u_{k-1} for k=2.
    else
    u   = G\(4*I*ukm1-I*ukm2);     % BDF2           

    %Update the previous result arrays for the next iteration.
    ukm2=ukm1; 
    ukm1=u;
    end

  if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',1.5); #Choosing which times to plot.
   title(['N = ' int2str(N) ',  n_{steps} = ' int2str([nstep])],'fontsize',15); pause(.1); end;
 
  uex = exp(time*lam_ex)*sx;
  err = uex-u;
  emx = max(abs(err))/max(abs(uex));
  el2 = sqrt( (err'*err) / (uex'*uex) );
  elast = enew;

end;

enew  = el2; ratio = elast/enew;

format longe
disp([name '    ' num2str([N nstep el2 ratio ],' %6d  %6d  %e  %e')])

