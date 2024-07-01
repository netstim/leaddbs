function MAC = ea_getMAC
% Get MAC address

if isunix
    % ifconfig should work for both Linux and macOS in most cases
    cmd = ['ifconfig | awk ''/^[a-z]/ ' ...
           '{if (iface && mac && ip) print iface, mac, ip; iface=$1; gsub(/:$/, "", iface); mac=""; ip=""} ' ...
           '/ether/ {mac=$2} /inet / && !/^127\./ {ip=$2} END' ...
           '{if (iface && mac && ip) print iface, mac, ip}'''];
    [status, cmdout] = system(cmd);
    if status
        % when ifconfig is not available, using ip instead (Linux only)
        cmd = ['ip -o link show | awk ''/ether/ {print substr($2, 1, length($2)-1), $(NF-2)}'' ' ...
               '| while read iface mac; do ip_address=$(ip -4 -o addr show "$iface" primary 2>/dev/null ' ...
               '| awk ''/inet/ {print $4}''); if [[ $ip_address ]]; then echo "$iface $mac ${ip_address%%/*}"; fi; done'];
        [~, cmdout] = system(cmd);
    end
    
    cmdout = strsplit(strip(cmdout), '\n')';
    cmdout(startsWith(cmdout, {'bridge', 'docker', 'tailscale', 'utun', 'vmnet'})) = [];
    cmdout = sort(cmdout);
    cmdout = split(cmdout{1}, ' ');
    MAC = cmdout{2};
elseif ispc
    [~, cmdout] = system('getmac /nh');
    cmdout = strip(cmdout);
    MAC = replace(lower(cmdout(1:17)), '-', ':');
end
