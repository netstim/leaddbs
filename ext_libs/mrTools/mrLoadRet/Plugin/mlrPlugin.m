% mlrPlugin.m
%
%        $Id:$ 
%      usage: mlrPlugin(<v>)
%         by: justin gardner
%       date: 11/24/10
%    purpose: With no arguments, brings up plugin selection dialog. This function will
%             look in the directory mrLoadRet/Plugin for plugins (i.e. directories
%             that have a function names dirnamePlugin.m), you can also specify
%             a comma-delimited list of alternative plugin directories to also 
%             search for plugins with the preference 'pluginPaths' - i.e. if you had
%             a group plugin path you could do: mrSetPref('pluginPaths','/path/to/group/plugins');
% 
%             It saves in the preference 'selectedPlugins', a cell array with the names of 
%             all currently selected plugins.
%
%             With a view argument initializes all plugins that the user has selected for that view
%             If you call with altPlugins set to 1 it will check the directories
%             under mrLoadRet/PluginAlt for more plugin directories. These are 
%             advanced / undocumented / untested plugins that can be chosen 
%
%             mlrPlugin('altPlugins=1');
%
%             If you add a plugin with altPlugins set above then it will add that
%             alt plugin's parent directory to your list of pluginPaths and if
%             you remove all alt plugins then it will rekove that alt plugins parent
%             directory from your pluginPaths. This is so that when you call with
%             mlrPlugin it will show / not show relevant PluginAlt directories
%     
function v = mlrPlugin(varargin)

v = [];
altPlugins = 0;
mlrPath('mrTools');

% check arguments
if ~any(nargin == [0 1 2])
  help mlrPlugin
  return
end

if nargin > 0
  % if the first argument is a view than we are initializing
  if isview(varargin{1})
    v = varargin{1};
    % when initializing, set altPlugins to true so that it checks there
    % for the plugins
    getArgs({varargin{2:end}},'altPlugins=1');
  else
    % with no arguments, setting plugins default to not showing altPlugins
    getArgs(varargin,{'altPlugins=0'});
  end
end

% check Plugin directory for possible plugins
% first, find where this function lives. Under that
% should be directories for plugins that are distributed
% with MLR
mlrPluginPath = fileparts(which('mlrPlugin'));
if ~isdir(mlrPluginPath)
  disp(sprintf('(mlrPlugin) Could not find default plugin directory'));
  return
end
% get the plugins from the default location within the MLR distribution
plugins = getPlugins(mlrPluginPath);

% alternative plugin directory name
pluginAltDirname = sprintf('%sAlt',fileparts(which('mlrPlugin')));

% check any directories found in the alternate plugin path
pluginPaths = commaDelimitedToCell(mrGetPref('pluginPaths'));

for i = 1:length(pluginPaths)
  % check if the first part is equal to PluginAlt in which case we append
  % the correct path
  if strncmp('PluginAlt',pluginPaths{i},length('PluginAlt'))
    pluginPaths{i} = fullfile(pluginAltDirname,getLastDir(pluginPaths{i}));
  end
end

% add the plugin paths
pluginPaths = unique(pluginPaths);
for i = 1:length(pluginPaths)
  % add the plugin path
  plugins = getPlugins(pluginPaths{i},plugins);
end

% check alt plugins if asked for
if altPlugins
  % Check mrLoadRet/PluginAlt directory
  if ~isdir(pluginAltDirname)
    disp(sprintf('(mlrPlugin) Could not find %s directory',pluginAltDirname));
  else
    % chekc the top directory for plugins
    plugins = getPlugins(pluginAltDirname,plugins,false);
    % check each sub-directory
    pluginAltDir = dir(pluginAltDirname);
    for i = 1:length(pluginAltDir)
      if pluginAltDir(i).isdir 
	dirName = pluginAltDir(i).name;
	if (length(dirName) >= 1) && (dirName(1) ~= '.')
	  dirName = fullfile(pluginAltDirname,dirName);
	  if isdir(dirName)
	    % load the plugins
	    plugins = getPlugins(dirName,plugins);
	  end
	end
      end
    end
  end
end

% no plugins - display message and return
if isempty(plugins)
  disp(sprintf('(mlrPlugin) No plugin directories found. Plugins should be in directory mrLoadRet/Plugins'));
  return
end

% get which plugins user has selected on, this should be a list
% of plugin names
selectedPlugins = mrGetPref('selectedPlugins');
for i = 1:length(selectedPlugins)
  % get which number plugins are selected
  selected = find(strcmp(selectedPlugins{i},{plugins.name}));
  % and select them.
  if length(selected) == 1
    plugins(selected).selected = true;
  elseif length(selected)>1
    disp(sprintf('(mlrPlugin) Multiple (%i) plugins names %s were found. Only selecting the first one (%s).',length(selected),selectedPlugins{i},plugins(selected(1)).path));
    plugins(selected(1)).selected = true;
  end
end

% when run with no view...
if isempty(v)
  % keep the list of which plugins were selected
  previousSelectedPlugins = selectedPlugins;
  % put up a selection dialog box
  for i = 1:length(plugins)
    paramsInfo{i} = {plugins(i).name,plugins(i).selected,'type=checkbox',plugins(i).help};
  end
  params = mrParamsDialog(paramsInfo,'Choose plugins you wish to install');
  % if user didn't hit cancel, then save the new choices
  if ~isempty(params)
    % get user choices
    selectedPlugins = {};
    for i = 1:length(plugins)
      if params.(plugins(i).name)
	selectedPlugins{end+1} = plugins(i).name;
      end
    end
    if ~isempty(getMLRView) && ~isequal(previousSelectedPlugins,selectedPlugins)
      mrWarnDlg(sprintf('(mlrPlugin) Restart of MLR is required for change in plugin list to take effect'));
    end
    % save the choices
    mrSetPref('selectedPlugins',selectedPlugins);
    % if we have altPlugins selected and the user has selected a plugin from
    % one of those directories, then make that a default plugin path for the user
    % so that it comes up when they call mlrPlugin from the menu. Also remove
    % any default pluginPaths for a choice that has been deselected
    if altPlugins
      % get which have been deselected
      deselectedPlugins = setdiff(previousSelectedPlugins,selectedPlugins);
      deselectedAltPath = {};
      % find the path (if it is in Alt) for each deselectedPlugin
      for i = 1:length(deselectedPlugins)
	% get the path for the selected plugin
	selectedPluginNum = find(strcmp(deselectedPlugins{i},{plugins(:).name}));
	% check if it is an Alt
	if strcmp(fileparts(fileparts(plugins(selectedPluginNum).path)),pluginAltDirname)
	  deselectedAltPath{end+1} = fullfile('PluginAlt',getLastDir(fileparts(plugins(selectedPluginNum).path)));
	end
      end
      deselectedAltPath = unique(deselectedAltPath);
      % find the path (if it is in Alt) for each selectedPlugin
      selectedAltPath = {};
      for i = 1:length(selectedPlugins)
	% get the path for the selected plugin
	selectedPluginNum = find(strcmp(selectedPlugins{i},{plugins(:).name}));
	% check if it is an Alt
	if strcmp(fileparts(fileparts(plugins(selectedPluginNum).path)),pluginAltDirname)
	  selectedAltPath{end+1} = fullfile('PluginAlt',getLastDir(fileparts(plugins(selectedPluginNum).path)));
	end
      end
      selectedAltPath = unique(selectedAltPath);
      % deselect any alt paths that do not have any plugins in them anymore
      deselectedAltPath = setdiff(deselectedAltPath,selectedAltPath);
      pluginPaths = commaDelimitedToCell(mrGetPref('pluginPaths'));
      pluginPaths = setdiff(pluginPaths,deselectedAltPath);
      % and add and selectedAltPaths
      pluginPaths = union(pluginPaths,selectedAltPath);
      % and save
      mrSetPref('pluginPaths',cellToCommaDelimited(pluginPaths));
    end
  end
% with a view argument, then run plugins
else
  tic
  % get which plugins are selected
  selectedPlugins = find([plugins.selected]);

  % and run their install command
  for i = 1:length(selectedPlugins)
    disp(sprintf('(mlrPlugin) Installing plugin %s',plugins(selectedPlugins(i)).name));
    retval = eval(plugins(selectedPlugins(i)).installCommand);
    if isview(retval)
      v=retval; %for those plugins that need to modify the view
    end
  end
  pluginTime=toc;
  disp(['(mlrPlugin) Installing Plugins took ' num2str(pluginTime) ' sec']);

end
  

%%%%%%%%%%%%%%%%%%%%
%    getPlugins    %
%%%%%%%%%%%%%%%%%%%%
function plugins = getPlugins(pluginPath,plugins,warnMissingPluginFunction)

% empty plugins directory to start with
if nargin == 1,plugins = [];end
if nargin < 3,warnMissingPluginFunction = true;end

% check for plugin function in each directory
pluginDir = dir(pluginPath);
if isempty(pluginDir)
  mrWarnDlg(['(mlrPlugin) Plugin directory ' pluginPath ' does not exist.']);
end
for i = 1:length(pluginDir)
  if pluginDir(i).isdir && (length(pluginDir(i).name) > 1) && (pluginDir(i).name(1)~='.')
    % check for proper file, that is there should be a file
    % called directoryPlugin.m which can be run to get help
    % info and install the plugin
    pluginName = sprintf('%sPlugin',pluginDir(i).name);
    pluginFullName = fullfile(pluginPath,pluginDir(i).name,sprintf('%s.m',pluginName));
    pluginFullName = mlrReplaceTilde(pluginFullName);
    if isfile(pluginFullName)
      % Make sure plugin function exists
      if isempty(which(pluginName))
	if warnMissingPluginFunction
	  disp(sprintf('(mlrPlugin) Plugin function %s does not exist on path',pluginFullName));
	end
	continue
      end
      % make sure that calling the function returns the correct one on the path
      if ~isequal(pluginFullName,which(pluginName))
	% just give a warning
	mrWarnDlg(sprintf('(mlrPlugin) Plugin function %s found in %s but should be in %s',pluginName,which(pluginName),pluginFullName));
      end
      % only add the plugin if it does not already exist
      if isempty(plugins) || ~any(strcmp(pluginDir(i).name,{plugins(:).name}))
	% now keep info about the plugin
	plugins(end+1).name = pluginDir(i).name;
	plugins(end).path = fileparts(pluginFullName);
	plugins(end).help = eval(sprintf('%s(''help'')',pluginName));
        % make sure we got a string back for the help
	if ~isstr(plugins(end).help)
	  disp(sprintf('(mlrPlugin) Plugin %s did not return help string correctly',pluginName));
	  plugins(end).help = 'No help information for plugin';
	end
	plugins(end).installCommand = sprintf('%s(''install'',v);',pluginName);
	plugins(end).selected = 0;
      end
    else
      if warnMissingPluginFunction
	disp(sprintf('(mlrPlugin) Plugin function %s not found in %s',pluginName,fileparts(pluginFullName)));
      end
    end
  end
end
