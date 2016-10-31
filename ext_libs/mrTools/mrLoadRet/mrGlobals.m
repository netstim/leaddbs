% mrGlobals script
%
% Defines MLR as a global variable.
% If MLR is not yet initialized then do so.
% Runs as a script in the scope of the calling function.
%
% djh 6/2004
% Just testing merge erase this line later

global MLR
global mrDEFAULTS

% If MLR is not yet initialized then do so
if isempty(MLR) || (isfield(MLR,'session') && isempty(MLR.session))
  
    % make sure paths are fixed to not to conflict with vista
    mlrPath mrTools

    % read the preferences and figlocs
    mrDEFAULTS = loadMrDefaults;

    % Check Matlab version number
    [mlrVersion, expectedMatlabVersion, expectedToolboxNames] = mrLoadRetVersion;
    version = ver('Matlab');
    matlabVersion = str2num(version.Version(1:min(length(version.Version),4)));
    if ~ismember(matlabVersion, expectedMatlabVersion);
      oneTimeWarning('mrToolsMatlabVersionError',['(mrGlobals) mrTools has been tested on Matlab versions ',mlrnum2str(expectedMatlabVersion,'compact=1'),'. You are running Matlab ',version.Version]);
    end

    % Check for expected toolboxes
    versionAll = ver;
    for iToolbox = 1:length(expectedToolboxNames)
      if ~any(strcmp(expectedToolboxNames{iToolbox},{versionAll.Name}))
	oneTimeWarning(sprintf('mrToolsMissing%s',expectedToolboxNames{iToolbox}),sprintf('(mrGlobals) mrTools uses the %s which you do not have installed - functions that rely on this toolbox may fail to work.',expectedToolboxNames{iToolbox}));
      end
    end
    
    % Initialize MLR
    MLR.version = mlrVersion;
    MLR.homeDir = pwd;

    % Load session and groups structures from mrSESSION.mat
    [session, groups] = loadSession(MLR.homeDir);
    % check session
    if isempty(session)
      oneTimeWarning(sprintf('mrSession_%s',fixBadChars(MLR.homeDir)),sprintf('(mrGlobals) Could not find mrSession in %s',MLR.homeDir));
    end
    MLR.session = session;
    MLR.groups = groups;

    % Initialize MLR.views
    MLR.views = {};

    % Initialize graph window
    MLR.graphFigure = [];

    % setup caches
    MLR.caches = {};
    
    % Inform user that mrLoadRet has started up
    oneTimeWarning('mrLoadRetVersion',['(mrGlobals) mrLoadRet ',num2str(MLR.version),', Matlab ',num2str(matlabVersion)],1);

    % Clean up
    clear expectedMatlabVersion version matlabVersion session groups mlrVersion 
end

