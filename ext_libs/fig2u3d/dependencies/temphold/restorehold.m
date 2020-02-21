function [] = restorehold(ax, prevhold)
%RESTOREHOLD   restore hold status.
%   RESTOREHOLD(ax, prevhold) restores the hold status of axes object with
%   handle ax, to the state defined by prevhold. To save prevhold for use
%   with this function, call takehold or ishold.
%
% usage
%   RESTOREHOLD(ax, prevhold)
%
% input
%   ax = axes object handle
%   prevhold = hold status to restore
%            = 1 to set hold(ax, 'on') | 
%              0 to set hold(ax, 'off')
%
% See also TAKEHOLD, HOLD, ISHOLD.
%
% File:      restorehold.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.01.22 - 2012.05.10
% Language:  MATLAB R2012a
% Purpose:   restore hold status of axes object
% Copyright: Ioannis Filippidis, 2012-

if prevhold == 1
    hold(ax, 'on')
elseif prevhold == 0
    hold(ax, 'off')
else
    error('hold was neither on or off.')
end
