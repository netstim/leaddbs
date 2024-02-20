function prefs = ea_setprefs(key, value, whichPrefs)
% Set LeadDBS preference
% 'whichPrefs' can be either 'machine' (default) to set 'prefs.machine.*'
% or 'user' to set 'prefs.*'

% In case deployed, only support set user prefs saved in ea_prefs_user.json
if isdeployed
    prefPath = ea_prefspath('json');
    prefs = loadjson(prefPath);
    fields = strsplit(key, '.');
    prefs = setfield(prefs, fields{:}, value);
    savejson('', prefs, prefPath);
    return
end

if ~exist('whichPrefs', 'var') || isempty(whichPrefs)
    % Set prefs.machine.* by default
    whichPrefs = 'machine';
end

switch lower(whichPrefs)
    case {'m', 'machine'}
        prefs = ea_prefs;
        machine = prefs.machine;
        fields = strsplit(key, '.');
        machine = setfield(machine, fields{:}, value);
        try % may not have write permissions
            save(ea_prefspath('.mat'),'machine');
        catch
            warning('Could not save preferences to user home directory. Please check permission.');
        end
    case {'u', 'user'}
        if ~ischar(value)
            error('Input value should be a string!');
        end

        if isempty(str2num(value)) || ~isempty(regexp(value, '[a-zA-Z]+', 'once'))
            % value to be set is str
            value = ['''', value, ''''];
        end

        % Replace prefs in ea_prefs_user.m
        prefs = fileread(ea_prefspath);
        pattern = [strrep(['prefs.',key], '.', '\.'), ' *= *.*?;'];
        if ~isempty(regexp(prefs, pattern, 'once')) % Key exists
            prefs = regexprep(prefs, pattern, regexptranslate('escape',['prefs.',key,' = ',value,';']));
        else % New key added
            prefs = [prefs, sprintf('\nprefs.%s = %s;\n',key,value)];
        end

        try % may not have write permissions
            fid = fopen(ea_prefspath, 'w');
            fwrite(fid, prefs);
            fclose(fid);
        catch
            warning('Could not save preferences to user home directory. Please check permission.');
        end
end
