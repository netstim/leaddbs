function bbox = ea_get_bbox(input, opts)
% Get the bounding box of the input image
% In case of multiple input images, will return the largest bbox covering
% all the images.

arguments
    input       {mustBeText}
    opts.tight  {mustBeNumericOrLogical} = false % Thresholding the image to get a minimum (tight) bbox
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
if opts.tight
    bbox = cellfun(@(x) spm_get_bbox(x, 'nz'), input, 'Uni', 0);
else
    bbox = cellfun(@(x) spm_get_bbox(x), input, 'Uni', 0);
end

bbox = cell2mat(bbox);

if length(input) > 1
    bbox = [min(bbox(1:2:end, :)); max(bbox(2:2:end, :))];
end

% Delete unzipped input
if ~isempty(gzInputs)
    ea_delete(erase(gzInputs, ".gz" + textBoundary('end')));
end
