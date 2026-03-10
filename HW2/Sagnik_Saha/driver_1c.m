%Make sure to specify Chebyshev/uniform in ad_var.m
Nvals = [32,64,128,256]

fprintf('%5s %10s %10s %10s %12s\n','N','dt','nsteps','s_*','rel err')
for loop=1:length(Nvals);
N = Nvals(loop);

Ts = 0.1:0.01:2;      % T values
L = 1;
dx = L/N;
c = 1;
cfl=0.1;
dt=0.1 * dx/c;

%s(t) array
st = zeros(length(Ts),1);
for j = 1:length(Ts)
    T = Ts(j);
    run('bd.m')
    %Take derivative
    w = D*u;
    st(j) = max(abs(w(:)));
end

[stmax, ind] = max(st);
stanalytic = 152.00516;
relerr = abs(stmax - stanalytic)/abs(stanalytic);

fprintf('%5d %10.4e %10d %10.4f %12.4e\n', N, dt, nstep, stmax, relerr);
end

