function ind = ea_sub2ind(sz, sub)
% Wrapper for sub2ind to support sub input as an array
if numel(sz) ~= size(sub, 2)
    error('Size of the matrix doesn''t match dimension of the subscripts!');
end

if any(sub(:)<=0)
    error('Subscripts must be positive values!');
end

sub = mat2cell(sub, size(sub,1), ones(1, size(sub,2)));
ind = sub2ind(sz, sub{:});
