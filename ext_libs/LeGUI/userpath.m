function p = userpath(inArg)
%USERPATH User environment path.
%   USERPATH returns a string specifying the first folder on the search path.
%   MATLAB adds the userpath to the search path upon startup.
%
%   The default userpath is platform-specific:
%     Windows: user's "Documents" folder appended with "MATLAB"
%     Mac:     user's "Documents" folder ($home/Documents) appended with "MATLAB"
%     Linux:   user's $home folder appended by "Documents" and "MATLAB"
%              (if there is no $home/Documents directory, userpath will not be set)
%
%   USERPATH(path) set the current value of userpath to the folder passed in.
%   It updates the current MATLAB path, and this new userpath will persist
%   across MATLAB sessions. 
%
%   USERPATH('reset') resets the userpath to the default.  It updates the current
%   MATLAB path, and this new userpath will persist across MATLAB sessions. 
%
%   USERPATH('clear') removes the userpath.  It updates the current MATLAB path,
%   and this new userpath will persist across MATLAB sessions.  
%
%   Note that additional folders can be added to the top of the search path upon
%   startup by specifying the path for the folders via the MATLABPATH environment
%   variable.
%
%   On Mac and Linux, if the directory $home/matlab exists, it will also be
%   added to the search path after any folders specified via the MATLABPATH
%   environment variable regardless of whether the MATLABPATH environment
%   variable is set. (This is done for compatibility with earlier releases.)
%
%   See also PATHDEF.

%   Copyright 1984-2017 The MathWorks, Inc.

narginchk(0, 1);

 if nargin > 0 
    validateattributes(inArg, {'string', 'char'}, {'scalartext'});
    inArg = convertStringsToChars(inArg);
 end

% If found, process argument and return.  
if nargin == 1 
    switch inArg
      case 'reset'
        resetUserPath;
      case 'clear'
        clearUserPath;
      otherwise
        setUserPath(inArg);
    end
    return
end

% Return the user work folder if it exists
p = system_dependent('getuserworkfolder');
if ~exist(p,'dir') || ~isAbsolute(p)
    p = '';
end

end % userpath

function resetUserPath
  oldUserPath = system_dependent('getuserworkfolder');
  rmpathWithoutWarning(oldUserPath);
  defaultUserPath = system_dependent('getuserworkfolder', 'default');
  addpath(defaultUserPath);
  s = matlab.settings.internal.settings;
  matlabNode = s.matlab;
  if matlabNode.UserPath.hasPersonalValue()
    matlabNode.UserPath.clearPersonalValue();
  end
end
 
function setUserPath(newPath)
  if exist(newPath, 'dir')
    % Insure that p is an absolute path
    if isAbsolute(newPath)
      oldUserPath = system_dependent('getuserworkfolder');
      rmpathWithoutWarning(oldUserPath);
      addpath(newPath);
      s = matlab.settings.internal.settings;
      s.matlab.UserPath.PersonalValue = newPath;
    else
      error(message('MATLAB:userpath:invalidInput'));
    end
  else
    error(message('MATLAB:userpath:invalidInput'));
  end
end

function clearUserPath
  oldUserPath = system_dependent('getuserworkfolder');
  rmpathWithoutWarning(oldUserPath);
  s = matlab.settings.internal.settings;
  s.matlab.UserPath.PersonalValue = '';
end


function rmpathWithoutWarning(pathToDelete)
  if ~isempty(pathToDelete)
    [lastWarnMsg, lastWarnId] = lastwarn;
    oldWarningState = warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(pathToDelete);
    warning(oldWarningState.state,'MATLAB:rmpath:DirNotFound')
    lastwarn(lastWarnMsg, lastWarnId);
  end
end

function status = isAbsolute(file)
  if strncmp(computer,'PC',2)   % ISPC is not available during MATLAB startup
    file = strrep(file,'\','/'); % Canonicalize file name
    status = ~isempty(regexp(file,'^[a-zA-Z]:\/','once')) || strncmp(file,'//',2);
  else
    status = strncmp(file,'/',1);
  end
end
