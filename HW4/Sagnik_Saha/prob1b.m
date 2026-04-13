%Define the problem parameters.
L = 1;
c = 1;
nu = 0.001;
s=0.7;

V = 1; %No. of vertices = 1 due to linear interpolation.
E = 20; %No. of elements. 
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
Inm1 = speye(N-1);
R = [zeros(size(Inm1, 1), 1), Inm1, zeros(size(Inm1, 1), 1)]; %This is for Dirichlet on both sides.


%%%%Define Ce as before, since it is independent of Le.
Ce = zeros(2,2);
Ce(1,1)=-1/2;
Ce(1,2)=1/2;
Ce(2,1)=-1/2;
Ce(2,2)=1/2;

%%%%%%%%%Define A_L, f_L, which depend on Le. Must be set manually. C_L set as before.%%%%%%%%%%
AL = zeros(2*E,2*E);
fL=zeros(2*E,1);

L1 = (1-s)/(1-s^E);

for e = 1:E
    Le = L1 *s^(e-1);
    ind = 2*(e-1) + (1:2);

    Ae = [ 1  -1
      -1   1 ];
    AL(ind, ind) = Ae / Le;

    fe = [Le/2; Le/2];
    fL(ind) = fL(ind) + fe; 
end

%CL values are independent of Le so can be set as before.
CL = kron(eye(E), Ce);

%%%%%%%%%%Define Abar, Bbar, fbar%%%%%%
Abar = Q' * AL * Q;
Cbar = Q' * CL * Q;

fbar = Q'*fL;

%%%%%%%%%% Construct A, B, f %%%%%%
A = R * Abar * R';
C = R * Cbar * R';
f = R*fbar;

%%%%Construct linear system to be solved. Lu = f. %%%%%%%%
L = C + nu * A;
res = L\f;

%%%%%%%%Define the xvals for analysis%%%%%%%%%
xvals = zeros(E+1,1);
xvals(1) = 0;
for e = 1:E
    Le = L1 *s^(e-1);
    xvals(e+1) = xvals(e) + Le;
end

%%%%%%%%%Plotting%%%%%%%%%%%
%u_full = zeros(E+1,1);  % initialize full solution
%u_full(2:end-1) = res;  
%plot(xvals, u_full, '-o','markersize',3)
%xlabel('x'); ylabel('u(x)');
%title('Steady-state FEM solution');
%grid on
%pause

%%%%%%%Pointwise error%%%%%%%%
function y=exact_steadysol(x,nu)
    y=x - (exp(-1/nu) - exp((x-1)/nu))/(exp(-1/nu) - 1);
end

%Evaluate the exact solution over the range of xvals. Note that we exclude the boundary terms x = 0 and 1.
exactsol = exact_steadysol(xvals(2:end-1), nu);

%Compute the pointwise error (again, ignoring boundary terms).
err_pntwise = abs((exactsol - res)./exactsol); 

%Display the maximum pointwise error.
disp(max(err_pntwise))

%Plot the pointwise error as a function of 'x'.
%plot(xvals(2:end-1),err_pntwise,'-o','markersize',2)
%pause

