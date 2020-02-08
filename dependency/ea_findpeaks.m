function [pks, locs, varargout] = ea_findpeaks(data, varargin)

has_findpeaks = ~isempty(which('has_findpeaks'));
if has_findpeaks
    [pks, locs, w, p] = findpeaks(data, varargin{:});
    varargout = {w, p};
else
    slope_sign = sign(diff(data));
    locs = 1 + find(slope_sign(1:end-1) > 0 & slope_sign(2:end) < 0);
    if ~isempty(varargin) && strcmpi(varargin{1}, 'NPeaks')
        locs = locs(1:varargin{2});
    end
    pks = data(locs);
end