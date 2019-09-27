function ext = ea_prefsext()
% Returns the extension of the preferences file (.m or .json)

if isdeployed
    ext = '.json';
else
    ext = '.m';
end
