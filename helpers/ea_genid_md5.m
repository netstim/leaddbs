function ID = ea_genid_md5(ID, width)
% Generate patient ID based on MD5 hashtag.
arguments
    ID {mustBeText} % Original ID
    width {mustBeNumeric, mustBeLessThanOrEqual(width, 32)} = 4  % width of ID
end

if ~iscell(ID)
    ID = genid(ID, width);
else
    ID = cellfun(@(x) genid(x, width), ID, 'Uni', 0);
end


function ID = genid(ID, width)
Hash = GetMD5(ID, '8bit', 'HEX');
ID = Hash(1:width);
