function bbox = ea_get_bbox(input, mode)
% Get [largest] bounding box of input image[s]

arguments
    input   {mustBeText}
    mode    {mustBeTextScalar, mustBeMember(mode, {'all', 'max'})} = 'max' % Get the largest bbox in case of multiple inputs
end

% Make sure input is cell
if ischar(input)
    input = {input};
end

% Make sure input is column vector
if isrow(input)
    input = input';
end

% Get full path of the input
input = GetFullPath(input);

% Unzip input incase necessray
gzInputs = input(endsWith(input, '.gz'));
if ~isempty(gzInputs)
    gunzip(gzInputs);
    input = erase(input, ".gz" + textBoundary('end'));
end

% Get bounding box
bbox = cellfun(@(x) spm_get_bbox(x), input, 'Uni', 0);

if isscalar(input)
    bbox = cell2mat(bbox);
elseif strcmp(mode, 'max')
    bbox = cell2mat(bbox);
    bbox = [min(bbox(1:2:end, :)); max(bbox(2:2:end, :))];
end

% Delete unzipped input
if ~isempty(gzInputs)
    ea_delete(erase(gzInputs, ".gz" + textBoundary('end')));
end

