function lead_sm = LeG_alignElecs(lead,d)
%"lead" is xyz coordinates with deepest contact first.
%"d" is center-to-center spacing in mm.
%"d" may be one scalar or one value per contact interval.

nContacts = size(lead,1);

if nContacts < 2
    lead_sm = lead;
    return;
end

d = double(d(:));

% A scalar spacing applies to every interval.
if isscalar(d)
    d = repmat(d, nContacts-1, 1);
elseif numel(d) ~= nContacts-1
    error(['Spacing must contain either one value or exactly %d values ' ...
           'for an electrode with %d contacts.'], ...
           nContacts-1, nContacts);
end

% Validate spacing
if any(~isfinite(d)) || any(d <= 0)
    error('Contact spacing must be positive and finite.');
end

% Cumulative distances from the deepest contact
s = [0; cumsum(d)];

% Center original coordinates and cumulative distances
center = mean(lead,1);
lead_centered = lead - center;
s_centered = s - mean(s);

% Find optimal straight-line direction using least squares
direction = s_centered' * lead_centered;
directionLength = norm(direction);

if directionLength < 1e-10
    error('Cannot determine electrode trajectory.');
end

direction = direction ./ directionLength;

% Reposition contacts along the fitted trajectory
% with exact spacing and minimum total squared displacement
lead_sm = center + s_centered * direction;

end