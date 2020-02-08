function y = ea_sgolayfilt(x, order, framelen, weights, dim)
% varargin: x, order, framelen, weights, dim

has_sgolayfilt = ~isempty(which('sgolayfilt'));
if has_sgolayfilt
    if nargin < 4
        y = sgolayfilt(x, order, framelen);
    elseif nargin < 5
        y = sgolayfilt(x, order, framelen, weights);
    else
        y = sgolayfilt(x, order, framelen, dim);
    end
else
    % Reimplementation from scipy version
    % https://github.com/scipy/scipy/blob/v0.17.1/scipy/signal/_savitzky_golay.py#L228-L349
    if nargin > 3
        warning('`weights` ignored.');
    end
    if nargin < 5 || isempty(dim)
        dim = ndims(x);
    end
    deriv = 0;
    delta = 1.0;
    
    %coeffs = savgol_coeffs(framelen, order, deriv, delta);
    if order >= framelen
        error('order must be less than framelen');
    end
    
    if mod(framelen, 2) == 0
        error('framelen must be odd.');
    end
    
    halflen = (framelen - 1) / 2;
    design_cols = (-halflen:framelen - halflen - 1);
    design_cols = flip(design_cols)';
    design_rows = (0:order);
    design_mat = design_cols .^ design_rows;
    design_mat = design_mat';
    
    % z determines which order derivative is returned.
    z = zeros(order + 1, 1);
    % The coefficient assigned to z[deriv] scales the result to take into
    % account the order of the derivative and the sample spacing.
    z(deriv + 1) = factorial(deriv) / (delta .^ deriv);
    
    % Find the least-squares solution of design_mat*coeffs = z
    % coeffs = design_mat' \ z;
    coeffs = lsqminnorm(design_mat, z);
            
    %savgol_filter(x, window_length, polyorder, deriv=0, delta=1.0,
    %              axis=-1, mode='interp', cval=0.0)
    % y = convolve1d(x, coeffs, axis=axis, mode="constant")
    y = zeros(size(x));
    for col_ix = 1:size(x, 2)
        y(:, col_ix) = conv(x(:, col_ix), coeffs, 'same');
        
        poly_coeffs = polyfit((1:framelen), x(1:framelen, col_ix)', order);
        y(1:halflen, col_ix) = polyval(poly_coeffs, (1:halflen)) / (delta .^ deriv);
        
        poly_coeffs = polyfit((1:framelen), x(end-framelen+1:end, col_ix)', order);
        interp_x = (framelen - halflen:framelen - 1);
        y(end-halflen+1:end, col_ix) = polyval(poly_coeffs, interp_x) / (delta .^ deriv);
    end
end
