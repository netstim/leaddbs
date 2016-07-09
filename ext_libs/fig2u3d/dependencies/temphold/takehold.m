function [prevhold] = takehold(ax, mode)
%TAKEHOLD  hold axes and return previous hold status.
%   TAKEHOLD(ax) is equivalent to hold(ax, 'on'), where ax is the handle to
%   an axes object.
%
%   TAKEHOLD(ax, 'local') clears the axes ax if their previous hold state
%   is off and then holds the axes, i.e., hold(ax, 'on').
%   This is useful when a plotting function draws multiple
%   graphics objects, so it needs to hold the plot intermediately, but
%   maintains pre-existing graphics in the plot only if hold(ax 'on') has
%   been issued by another source. This is the same behavior that the plot
%   function exhibits.
%
%   TAKEHOLD(ax, 'add') or TAKEHOLD(ax) does not clear the axes ax, but
%   instead holds it. It adds to the plot, without erasing previous
%   graphics.
%
%   prevhold = TAKEHOLD(ax, ...) returns the previous hold status of the
%   axes object ax. This can then be restored by calling 
%
% usage
%   prevhold = TAKEHOLD
%   prevhold = TAKEHOLD(ax)
%   prevhold = TAKEHOLD(ax, 'add')
%   prevhold = TAKEHOLD(ax, 'local')
%
% input
%   ax = axes object handle
%   mode = 'add' | 'local'
%
% output
%   held = 1 if ax already held
%          0 if ax not already held
%
% See also RESTOREHOLD, HOLD, ISHOLD.
%
% File:      takehold.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.01.22 - 2012.05.10
% Language:  MATLAB R2012a
% Purpose:   hold(ax, 'on') and return current hold status
% Copyright: Ioannis Filippidis, 2012-

% which axes' handle ?
if nargin < 1
    ax = gca;
end

% which mode ?
if nargin < 2
    mode = 'add';
end

prevhold = ishold(ax);

if prevhold == 1
    return
end

if strcmp(mode, 'local')
    disp('heh')
    cla(ax)
end

hold(ax, 'on')
