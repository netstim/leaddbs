function dicm_val_rep(dcmFile, tag, val, newFile)
% Replace tag value with val in dicmFile, and save as newFile.
%
%  dicm_val_rep('myFile.dcm', 'PatientName', 'sub-002', 'newFile.dcm');
%  dicm_val_rep('myFile.dcm', '0177,0011', 'sub-002'); % use hex for tag
%
% The first three input argument are manditory. The forth one, if left out or
% empty, will be dcmFile appended with '_new'.
% 
% See also anonymize_dicm

% 230226 Wrote it (xiangrui.li@gmail.com)

if nargin<4 || isempty(newFile)
    [f, nam, ext] = fileparts(dcmFile);
    newFile = fullfile(f, [nam '_new' ext]);
end

a = regexp(tag, '[0-9a-fA-F]{4}', 'match');
if numel(a)==2, tag = hex2dec(a)'; end % private tag like (0011,001a)

[s, err] = dicm_hdr(dcmFile, 'Manufacturer');
if isempty(s), error(err); end
fid = fopen(s.Filename, 'r', 'l');
b8 = fread(fid, inf, 'uint8=>uint8')';
fclose(fid);
try n = s.PixelData.Start; catch, n = s.FileSize; end

try
    tsUID = s.TransferSyntaxUID;
    be = strcmp(tsUID, '1.2.840.10008.1.2.2');
    expl = ~strcmp(tsUID, '1.2.840.10008.1.2');
    if ~expl && isnumeric(tag)
        error('Cannot determine value type of private tag for implicit VR');
    end
catch
    be = false;
    expl = true;
end

if isnumeric(tag)
    dict = dicm_dict(s.Manufacturer);
    ind = tag(1)*65536+tag(2) == dict.tag;
    if any(ind), vr = dict.vr{ind}; else, vr = ''; end
    tag2 = uint16(tag);
    tag8 = typecast(tag2, 'uint8');
    if be && tag2(1)>2, tag8 = tag8([2 1 4 3]); end
else
    dict = dicm_dict(s.Manufacturer, tag);
    if isempty(dict.group), error('In valid tag: "%s"', tag); end
    for j = numel(dict.group):-1:1
        tag2 = [dict.group(j) dict.element(j)];
        tag8 = typecast(tag2, 'uint8');
        if be && tag2(1)>2, tag8 = tag8([2 1 4 3]); end
        i = strfind(b8(1:n), tag8); i = i(mod(i,2)==1);
        if ~isempty(i), vr = dict.vr{j}; break; end
    end
end

hasVR = expl || tag2(1)==2;
if be, ed = 'b'; else, ed = 'l'; end
i = strfind(b8(1:n), tag8); i = i(mod(i,2)==1);
if numel(i)>1 && ~isempty(vr)
    i = strfind(b8(1:n), [tag8 uint8(vr)]);
    i = i(mod(i,2)==1);
end
if n==s.FileSize || numel(i)>3, i = i(1); % take 1st for multiframe
else, i = i(end);  % last is effective
end
if hasVR, vr = char(b8(i+4+(0:1))); end % non-expl already error out
i = i + 4; % 4-byte tag

nLen = 4; % # of bytes for value length
len16 = 'AE AS AT CS DA DS DT FD FL IS LO LT PN SH SL SS ST TM UI UL US';
if ~hasVR % implicit, length irrevalent to vr
    fmtLen = 'uint32';
elseif ~isempty(strfind(len16, vr)) %#ok length in uint16
    i = i + 2; % VR
    fmtLen = 'uint16';
    nLen = 2;
else % skip 2 bytes then length in uint32
    i = i + 4; % VR & extra 2 bytes
    fmtLen = 'uint32';
end
nB = typecast(b8(i+(0:nLen-1)), fmtLen);
if be && tag2(1)~=2, nB = swapbytes(nB); end

if isstruct(val) && strcmp(vr, 'PN'), val = strjoin(struct2cell(val));
elseif iscellstr(val), val = sprintf('%s\n', val{:}); %#ok
end

switch vr
    case {'AE' 'AS' 'CS' 'DA' 'DT' 'LO' 'LT' 'PN' 'SH' 'ST' 'TM' 'UI' 'UT'}
               bpv = 1; fmt = 'char*1'; val = strtrim(val);
    case 'IS', bpv = 1; fmt = 'char*1'; val = sprintf('%.0f\\', val);  val(end) = '';
    case 'DS', bpv = 1; fmt = 'char*1'; val = sprintf('%.16g\\', val); val(end) = '';
    case {'US' 'AT' 'OW'},  bpv = 2; fmt = 'uint16';
    case 'SS',              bpv = 2; fmt = 'int16';
    case 'UL',              bpv = 4; fmt = 'uint32';
    case 'SL',              bpv = 4; fmt = 'int32';
    case {'FL' 'OF'},       bpv = 4; fmt = 'single';
    case {'FD' 'OD'},       bpv = 8; fmt = 'double';
    otherwise,              bpv = 1; fmt = 'uint8'; % 'OB' or unknown
end
if bpv == 1 % uint8 or char
    if mod(numel(val),2)==1, val = [val(:); 0]; end
elseif be && tag2(1)==2
    val = swapbytes(cast(val, fmt));
end
len = numel(val) * bpv;
if be && tag2(1)==2, len = swapbytes(cast(len, fmtLen)); end

fid = fopen(newFile, 'w', ed);
fwrite(fid, b8(1:i-1), 'uint8'); % till len
fwrite(fid, len, fmtLen); % len: endian related
fwrite(fid, val, fmt);
fwrite(fid, b8(i+nLen+double(nB):end), 'uint8'); % skip len & val
fclose(fid);
