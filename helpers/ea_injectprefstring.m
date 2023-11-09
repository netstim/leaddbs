function ea_injectprefstring(varargin)
% Adds preference to ea_prefs with the following structure:
% prefs.varargin{1}.varargin{...}.varargin{n-1} = 'varargin{n}';

% check inputs
if nargin < 2
    error('Minimum number of inputs is 2')
end
if ~iscellstr(varargin)
    error('All inputs must be character vectors')
end

if ~isdeployed
    % generate output string
    str = ['prefs.', strjoin(varargin(1:end-1),'.'), '=''',varargin{end},''';'];
    % print
    fid=fopen(ea_prefspath,'a');
    fprintf(fid,'%s\n',str);
    fclose(fid);
else
    try
        % load existing prefs
        fid = fopen(ea_prefspath('json'),'rt');
        prefs = jsondecode(fread(fid,'*char')'); fclose(fid);
        % add new field
        prefs = setfield(prefs,varargin{1:end-1},varargin{end});
        % save
        fid = fopen(ea_prefspath(ea_prefsext),'wt');
        fwrite(fid, jsonencode(prefs), 'char'); fclose(fid);
    catch
        warning('User preferences file could not be read. Please set write permissions to Lead-DBS install directory accordingly.');
    end
end