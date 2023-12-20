function listing = dir_without_dots(varargin)

if nargin == 0
    inputFolder = pwd;
elseif nargin == 1
    inputFolder = varargin{1};
else
    error('Too many input arguments.')
end

listing = dir(inputFolder);
listing(startsWith({listing.name}', '.')) = [];
