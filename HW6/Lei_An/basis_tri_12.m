function phi = basis_tri_12(r,s,k,order)

% basis_tri_12: triangular basis functions of order 1 or 2
%
% order = 1:
%   k = 1,2,3
%
% order = 2:
%   k = 1,2,3,4,5,6

if order == 1

    if k == 1
        phi = (1 - r - s);
    elseif k == 2
        phi = r;
    elseif k == 3
        phi = s;
    else
        error('Invalid local basis index k for order=1.');
    end

elseif order == 2

    if k == 1
        phi = (1 - r - s).*(1 - 2*r - 2*s);
    elseif k == 2
        phi = r.*(2*r - 1);
    elseif k == 3
        phi = s.*(2*s - 1);
    elseif k == 4
        phi = 4.*r.*(1 - r - s);
    elseif k == 5
        phi = 4.*r.*s;
    elseif k == 6
        phi = 4.*s.*(1 - r - s);
    else
        error('Invalid local basis index k for order=2.');
    end

else
    error('Stopping. order > 2 not supported in basis_tri_12.m')
end