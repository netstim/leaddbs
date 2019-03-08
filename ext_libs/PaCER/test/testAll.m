% request explicitly from the user to launch test suite locally
if isempty(strfind(getenv('HOME'), 'jenkins'))
   reply = '';
   while isempty(reply)
       reply = input([' -> Do you want to launch the test suite locally? Time estimate: more than 60 minutes Y/N: '], 's');
   end

   if strcmpi(reply, 'y') || strcmpi(reply, 'yes')
       launchTestSuite = true;
   else
       launchTestSuite = false;
   end
else
   % on the CI, always reset the path to make absolutely sure, that we test
   % the current version
   restoredefaultpath()
   launchTestSuite = true;
end

% save the current folder
origDir = pwd;

% if the location of pacer is not yet known
if isempty(which('SETUP_PACER.m'))
   % move back to the root of the repository
   cd([fileparts(which('testAll.m')) filesep '..'])

   % assign the path
   PACERDIR = pwd;
else
   PACERDIR = fileparts(which('SETUP_PACER.m'));
   cd(PACERDIR);
end

% include the root folder and all subfolders.
addpath(genpath([pwd filesep 'test']));

try
   if launchTestSuite

      % define a success exit code
      exit_code = 0;

      % ensure that we ALWAYS call exit
      if ~isempty(strfind(getenv('HOME'), 'jenkins')) || ~isempty(strfind(getenv('USERPROFILE'), 'jenkins'))
         exit(exit_code);
      end
   end
catch ME
   if ~isempty(strfind(getenv('HOME'), 'jenkins')) || ~isempty(strfind(getenv('USERPROFILE'), 'jenkins'))
       % only exit on jenkins.
       exit(1);
   else
       % switch back to the folder we were in and rethrow the error
       cd(origDir);
       rethrow(ME);
   end
end