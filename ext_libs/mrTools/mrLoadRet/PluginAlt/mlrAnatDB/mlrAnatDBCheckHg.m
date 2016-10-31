% mlrAnatDBCheckHg.m
%
%        $Id:$ 
%      usage: tf = lrAnatDBCheckHg()
%         by: justin gardner
%       date: 06/22/15
%    purpose: Returns whether mercurial is correctly installed or not
%
function tf = mlrAnatDBCheckHg

tf = false;

% check for hg
[status,result] = system('which hg');
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) You do not have mercurial installed. You will need to install mercurial. Typicaly by going to the website: http://mercurial.selenic.com and following download instructions.'));
  return
end
% check here for config stuff
[status,result] = system('hg config');
if ~strfind(result,'extensions.largefiles')
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Your hg config needs to have the extension for largefiles enabled. This is done by running "hg config --edit" and adding the line "largefiles =" after the section header "[extensions]"'));
  return
end
[status,result] = system('hg config ui.username');
if status ~= 0
  mrWarnDlg(sprintf('(mlrAnatDBPlugin) Your hg config needs to have your username and email address for you to commit. Use "hg config --edit" to fix this.'));
  return
end
  
tf = true;

