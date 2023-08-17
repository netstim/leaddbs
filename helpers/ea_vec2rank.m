function rank = ea_vec2rank(vector, method)
% Calculate rank of a vector
% Ties are adjusted accordingly based on the method (average or min)

if ~exist('method', 'var')
    method = 'average';
end

switch method
    case 'average' % Most commonly used ranking
        rank = tiedrank(vector);
    case 'min' % "competition" ranking
        sorted = sort(vector);
        [~, rank] = ismember(vector, sorted);
    otherwise
        error('Method is not supported!');
end
