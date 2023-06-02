function clearXattr(folder, pattern, type)
% Clear extended attribute (most importantly 'com.apple.quarantine') of
% excutable or app to bypass Gatekeeper on macOS.

if nargin < 2
    pattern = '*maci64';
end

if nargin < 3
    type = 'f';
end

system(['find "', folder, '" -type ', type, ' -name "', pattern, '" -exec xattr -cr {} \;']);
