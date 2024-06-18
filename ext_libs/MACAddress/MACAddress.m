function [mac, st] = MACAddress(allMac)
% [mac, st] = MACAddress()
% 
% The default is to return one MAC address, likely for ethernet adaptor. If the
% optional input is provided and true, all MAC address are returned in cellstr.
% No internet connection is required for this to work.
% 
% The optional 2nd output, if requested, is a struct with following fields:
%  st.FriendlyName (meaningful for Windows only)
%  st.Description  (OS dependent description)
%  st.MAC_address  (the same order as the 1st output)
%  st.IPv4_address (empty if not available)
%  st.IPv6_address (empty if not available)
% 
% Examples:
%  mac = MACAddress(); % return 1st MAC in string
% The format is like F0-4D-A2-DB-00-37 for Windows, f0:4d:a2:db:00:37 otherwise.
% 
%  macs = MACAddress(1); % return all MAC on the computer
% The output is cell even if only one MAC found.
% 
%  [macs, st] = MACAddress(1); % also return more info in st
% 
% To get numeric:
%  num = uint8(sscanf(MACAddress, '%2x%*c', 6))';

% 170510 Adapted this from RTBox code (Xiangrui.Li at gmail.com).
% 170525 Include mex for MS Windows.
% 171030 mex.c more robust. Include Octave 4 mex.
% 180525 use jsystem('ipconfig') for Windows; java method moved behind.
% 180626 implement 2nd optional output for both m and mex (almost rewritten).

if nargin<1 || isempty(allMac), allMac = false; end % default to first MAC

if ispc
    [~, str] = jsystem({'ipconfig.exe' '/all'});
	str = regexprep(str, '\r', '');
    mac_expr = 'Physical Address.*?:\s*((?:[0-9A-F]{2}-){5}[0-9A-F]{2})\s';
    nam_expr = '\nEthernet adapter\s+(.*?):?\n'; % nam/des/ip4/ip6 all in a block
    des_expr = 'Description.*?:\s*(.*?)\n';
    ip4_expr = 'IP(?:v4)? Address.*?:\s*((?:\d{1,3}\.){3}\d{1,3})';
    ip6_expr = 'IPv6 Address.*?:\s*((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    fmt = '%02X-%02X-%02X-%02X-%02X-%02X'; % adopt OS format preference
elseif ismac
    % [~, str] = jsystem({'networksetup' '-listallhardwareports'});
    [~, str] = jsystem({'ifconfig'});
    mac_expr = '\n\s+ether\s+((?:[0-9a-f]{2}:){5}[0-9a-f]{2})\s';
    des_expr = '\n(.*?):\s+';
    nam_expr = des_expr;
    ip4_expr = 'inet\s+((?:\d{1,3}\.){3}\d{1,3})';
    ip6_expr = 'inet6\s+((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    fmt = '%02x:%02x:%02x:%02x:%02x:%02x';
else % linux
    [err, str] = jsystem({'ip' 'address'}); % later Linux
    if ~err % almost always
        mac_expr = '\s+link/ether\s+((?:[0-9a-f]{2}:){5}[0-9a-f]{2})\s';
        des_expr = '\n\d+:\s+(.*?):\s+';
        ip4_expr = '\s+inet\s+((?:\d{1,3}\.){3}\d{1,3})';
        ip6_expr = '\s+inet6\s+((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    else % use ifconfig for old linux
        cmd = '/sbin/ifconfig';
        if ~exist(cmd, 'file'), cmd = 'ifconfig'; end
        [~, str] = jsystem({cmd});
        mac_expr = '\s+HWaddr\s+((?:[0-9a-f]{2}:){5}[0-9a-f]{2})\s';
        des_expr = '\n*(.*?)\s+';
        ip4_expr = 'inet addr:\s*((?:\d{1,3}\.){3}\d{1,3})';
        ip6_expr = 'inet6 addr:\s*((?:[0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})';
    end
    nam_expr = des_expr;
    fmt = '%02x:%02x:%02x:%02x:%02x:%02x';
end

if allMac, [mac, ind] = regexp(str, mac_expr, 'tokens', 'start');
else,      [mac, ind] = regexp(str, mac_expr, 'tokens', 'start', 'once');
end
mac = [mac{:}];
% if iscell(mac) && numel(ind)>1 % make mac unique
%     [~, ia] = unique(mac);
%     ia = sort(ia); % [mac, ia] = unique(mac, 'stable'); ind = ind(ia);
%     mac = mac(ia);
%     ind = ind(ia);
% end

if nargout>1 && ~isempty(mac)
    st = struct('FriendlyName', [], 'Description', [], 'MAC_address', mac, ...
                'IPv4_address', [], 'IPv6_address', []);
    i0 = [1 regexp(str, '\n\S') numel(str)]; % split str into blocks
    for i = 1:numel(ind)
        j = find(i0<ind(i), 1, 'last');
        a = str(i0(j) : i0(j+1)); % the block with current mac
        c = regexp(a, nam_expr, 'tokens', 'once');
        if ~isempty(c), st(i).FriendlyName = c{1}; end
        c = regexp(a, des_expr, 'tokens', 'once');
        if ~isempty(c), st(i).Description  = c{1}; end
        c = regexp(a, ip4_expr, 'tokens', 'once');
        if ~isempty(c), st(i).IPv4_address = c{1}; end
        c = regexp(a, ip6_expr, 'tokens', 'once');
        if ~isempty(c), st(i).IPv6_address = c{1}; end
    end
end

% java method is OS-independent, more reliable than regexp, but often slower and
% miss eth1 seen at least 1 Ubuntu machine
if isempty(mac)
try %#ok
    if nargout>1, st = []; end
    if allMac, mac = {}; end
    ni = java.net.NetworkInterface.getNetworkInterfaces;
    while ni.hasMoreElements
        aa = ni.nextElement;
        a = aa.getHardwareAddress;
        if numel(a)~=6 || all(a==0), continue; end % not valid mac
        a = typecast(a, 'uint8'); % from int8
        m = sprintf(fmt, a);
        if nargout>1
            st(end+1).FriendlyName = char(aa.getName); %#ok
            st(end).Description = char(aa.getDisplayName);
            st(end).MAC_address = m;
            aa = aa.getInetAddresses;
            while aa.hasMoreElements
                c = char(aa.nextElement);
                a = regexp(c, '(\d{1,3}\.){3}\d{1,3}', 'match', 'once');
                if ~isempty(a), st(end).IPv4_address = a; end
                a = regexp(c, '(([0-9a-f]{0,4}:){1,7}[0-9a-f]{0,4})', 'match', 'once');
                if ~isempty(a)
                    st(end).IPv6_address = strrep(a, 'fe80:0:0:0:', 'fe80::'); 
                end
            end
        end
        if allMac, mac{end+1} = m; %#ok
        else, mac = m; break; % done after finding 1st
        end
    end
end
end

% If all attemps fail, give warning and return a random MAC
if isempty(mac)
    warning('MACAddress:RandomMAC', 'Returned MAC are random numbers');
    a = randi(255, [1 6], 'uint8');
    a(5) = bitset(a(5), 1); % set 8th bit for random MAC, likely meaningless
    mac = sprintf(fmt, a);
    if nargout>1
        st = struct('FriendlyName', [], ...
            'Description', 'Failed to find network adapter', ...
            'MAC_address', mac, 'IPv4_address', [], 'IPv6_address', []);
    end
    if allMac, mac = {mac}; end
end

%% faster than system: based on https://github.com/avivrosenberg/matlab-jsystem
function [err, out] = jsystem(cmd)
% cmd is cell str, no quotation marks needed for file names with space.
try
    pb = java.lang.ProcessBuilder(cmd);
    pb.redirectErrorStream(true); % ErrorStream to InputStream
    process = pb.start();
    scanner = java.util.Scanner(process.getInputStream).useDelimiter('\A');
    if scanner.hasNext(), out = char(scanner.next()); else, out = ''; end
    err = process.exitValue; % err = process.waitFor() may hang
    if err, error('java.lang.ProcessBuilder error'); end
catch % fallback to system() if java fails like for Octave
    cmd = regexprep(cmd, '.+? .+', '"$0"'); % double quotes if with middle space
    [err, out] = system(sprintf('%s ', cmd{:}, '2>&1')); % Octave need 2>&1
end
