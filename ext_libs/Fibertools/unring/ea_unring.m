%%%%%%%%
%
% unring - tool for removal of the Gibbs ringing artefact
% Usage: outvol = unring(invol,params)
% Options: invol - input volume 
%          params - 3x1 array with [minW maxW nsh]
%                     nsh discretization of subpixel spaceing (default 20)
%                     minW  left border of window used for TV computation (default 1)
%                     maxW  right border of window used for TV computation (default 3)


function v = unring(v,params)

if nargin == 1,
    params = [1 3 20];
end;

v = double(v);
v = ringRm(v,params);