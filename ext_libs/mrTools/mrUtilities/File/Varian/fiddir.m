% fiddir.m
%
%        $Id$ 
%      usage: fiddir()
%         by: justin gardner
%       date: 12/16/09
%    purpose: Lists info about all fids in the dir
%
function retval = fiddir(dirname)

% check arguments
if ~any(nargin == [0 1])
  help fiddir
  return
end

if nargin == 0, dirname = '.';end

d = dir(fullfile(dirname,'*.fid'));

for i= 1:length(d)
  [xform fidinfo] = fid2xform(fullfile(dirname,d(i).name));
  disp(sprintf('%s (%s): %s [%i %i %i %i]',fidinfo.startDatestr,fidinfo.elapsedTimeStr,d(i).name,fidinfo.dim(1),fidinfo.dim(2),fidinfo.dim(3),fidinfo.dim(4)));
end

  
