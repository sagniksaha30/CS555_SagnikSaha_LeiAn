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
Inm1 = speye(N-1);
R = [zeros(size(Inm1, 1), 1), Inm1, zeros(size(Inm1, 1), 1)]; %This is for Dirichlet on both sides.

%%%%%%%Define the x values (global nodes) and the Le's. %%%%%%
%Uniform spacing.
xvals = linspace(0, 1, E+1);
Le = xvals(2)-xvals(1);
%Levals = (xvals(2)-xvals(1)) * ones(1, length(xvals));

%Non-uniform spacing.

%%%%%%%Define A^e, B^e, f^e (C^e not needed for this part) %%%%%%
Ae = zeros(2,2);
Ce = zeros(2,2);

fe = (Le/2)*ones(2, 1);

Ae(1,1)=1/Le;
Ae(1,2)=-1/Le;
Ae(2,1)=-1/Le;
Ae(2,2)=1/Le;

%Be(1,1)=Le/3;
%Be(1,2)=-Le/6;
%Be(2,1)=-Le/6;
%Be(2,2)=Le/3;

Ce(1,1)=-1/2;
Ce(1,2)=1/2;
Ce(2,1)=-1/2;
Ce(2,2)=1/2;

%%%%%%%%%Define A_L, B_L, f_L%%%%%%%%%%
AL = kron(eye(E), Ae);
CL = kron(eye(E), Ce);

fL = repmat(fe, E, 1);

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

%%%%%%%%%Plotting the numerical solution%%%%%%%%%%%
%u_full = zeros(E+1,1);  % initialize full solution including boundary terms
%u_full(2:end-1) = res;  %
%plot(xvals, u_full, '-o','markersize',3)
%xlabel('x'); ylabel('u(x)');
%grid on
%pause

%%%%%%%Pointwise error%%%%%%%%
%Define the exact solution
function y=exact_steadysol(x,nu)
    y=x - (exp(-1/nu) - exp((x-1)/nu))/(exp(-1/nu) - 1);
end

%Evaluate the exact solution over the range of xvals. Note that we exclude the boundary terms x = 0 and 1.
xvals = xvals'; %Convert xvals to column vector since we originally used linspace.
exactsol = exact_steadysol(xvals(2:end-1), nu);

%Compute the pointwise error (again, ignoring boundary terms).
err_pntwise = abs((exactsol - res)./exactsol); 

%Display the maximum pointwise error.
disp(max(err_pntwise))

%Plot the pointwise error as a function of 'x'.
plot(xvals(2:end-1),err_pntwise,'-o','markersize',2)
pause

