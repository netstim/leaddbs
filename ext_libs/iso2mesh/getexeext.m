function exesuff=getexeext()
%
% exesuff=getexeext()
%
% get meshing external tool extension names for the current platform
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% output:
%     exesuff: file extension for iso2mesh tool binaries
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if ispc
    exesuff = '.exe';
else
    exesuff = ['.', computer('arch')];
end
