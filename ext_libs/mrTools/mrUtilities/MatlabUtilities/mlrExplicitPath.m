% mlrExplicitPath.m
%
%        $Id:$ 
%      usage: p = mlrExplicitPath(p)
%         by: justin gardner
%       date: 02/19/12
%    purpose: Makes a path explicit by removing ~ or .. and replacing
%
function p = mlrExplicitPath(p)

% check arguments
if ~any(nargin == [1])
  help mlrExplicitPath
  return
end

p = mlrReplaceTilde(p);
if (length(p) > 2) && isequal(p(1:2),'..')
  % replaece ..
  p = fullfile(fileparts(pwd),p(3:end));
end


