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

lead_sm = zeros(size(lead));
lead_sm(1,:) = lead(1,:);

for j = 1:nContacts-1

    a = lead_sm(j,:);

    if j < nContacts-1
        b = (lead(j+1,:) + lead(j+2,:)) ./ 2;
    else
        b = lead(j+1,:);
    end

    direction = b-a;
    directionLength = norm(direction);

    if directionLength == 0
        error('Cannot align contacts because two trajectory points overlap.');
    end

    v = direction ./ directionLength;

    lead_sm(j+1,:) = a + d(j).*v;
end
end