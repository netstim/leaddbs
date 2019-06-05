function close = ea_allclose(a, b, rtol, atol)
% Determine if two arrays are element-wise equal within a tolerance.

if nargin < 3
    rtol = 1e-05;
end
if nargin < 4
    atol = 1e-08;
end

close = all( abs(a(:)-b(:)) <= atol+rtol*abs(b(:)) );
