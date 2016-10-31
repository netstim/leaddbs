function pp_ = dsErrorsurface(x, y, e, surfaceColor, surfaceAlpha)
%DSERRORSURFACE - Error surface plot.
%
%   pp_ = dsErrorsurface(x, y, e, [surfaceColor], [surfaceAlpha])
%
%   plots color surface instead of errorbars. looks more 
%   attractive in ds's opinion. the surface covers [y-e] to [y+e], 
%   i.e. these are symmetrical errorbars. 
%
%   x, y            - data vectors
%   e               - std/se [error]
%                       [if e is 2 x length(y), it is interpreted as
%                       [lower; upper] and the asymmetric errors
%                       are plotted
%
%   surfaceColor    - 1x3 color-vector (e.g. gray is [0.5, 0.5, 0.5])
%                       (default value is gray)
%   surfaceAlpha    - scalar alpha between 0 and 1 (sets translucence)
%                       (default value is 1)
%                     [on Apple X11, anything but 1 causes printing
%                     problems due to a ML/X11 bug]
%
%         For example,
%            x = 0:0.2:(2*pi);
%            y = sin(x);
%            e = std(y)*ones(size(x));
%            figure, plot(x,y) 
%            hold on
%            dsErrorsurface(x,y,e)
%            hold off            
%         draws symmetric error surface of unit standard deviation.
%
% ds 10/26/03 - written
% ds 01/25/04 - wrote some more comments
% ds 02-11-04 - added asymmetric errorsurface [e needs to be Mx2]

% default surfaceColor
if (~exist('surfaceColor','var') )
    surfaceColor = [1 1 1]*0.6; 
end

% default surfaceAlpha
if (~exist('surfaceAlpha','var') )
    surfaceAlpha = 0.6; % default
end

% add some logic to deal with asymmetric errors
% x, y, e <-- this is assumed
% x, y, l,u <-- add stuff to deal this one

% check sizes of input
if (prod(size(x)) ~= prod(size(y))) 
 error('check the sizes of x and y - they look off! echume provlima')
end

%  prod(size(e))/2
% prod(size(y))

if  ( prod(size(y)) ~= prod(size(e)) ) & ( prod(size(y)) ~= prod(size(e))/2 ) 
 error('e must be length(Y) or length(Y)x2');
end

% check that the x/y/e are [.....] rather than [....]'
% so that fliplr works below.
% vectors needs to be 1xN 

if size(x,1) > size(x,2)
    x = x';
end

if size(y,1) > size(y,2)
    y = y';
end

if size(e,1) > size(e,2)
    e = e';
end

% make one big polygon...
%
%   <-- fliplr[x]--| y+e 
%                  |    
%   [x] ---------->| y-e

X = [ x fliplr(x) ];
if size(e,1) == 1
    Y = [y-e fliplr(y+e)];
elseif size(e,1) == 2 % asymmetric error bars
    % Y = [y-e(2,:) fliplr(y+e(1,:))]
    Y = [y-abs(e(1,:)) fliplr( y+abs(e(2,:)) )];
end  
% now use patch to make the plot...
% pp_  the handles get returned to they can be modified outside the
% function.

pp_ = patch(X,Y,surfaceColor);
set(pp_,'FaceAlpha', surfaceAlpha, ...
    'edgecolor', [0, 0, 0], 'edgealpha', 0)
% edge color and alpha are 0 by default...

return;

%%%%%%%%%%%%%%%%%%%%%%%%
% errorsurface example %
%%%%%%%%%%%%%%%%%%%%%%%%

% make some data.
x = linspace(0, 3*pi, 30);
y1 = sin(x) + rand(size(x));
y2 = cos(x) + 0.1*rand(size(x));
e1 = 0.2+rand(size(x))./2;
e2 = 0.2+rand(size(x))./2;


figure(1)
% pass handles into variables that you keep:
h1_ = plot(x, y1, 'ko-');
set(h1_,'markerfacecolor',[0 0 0], 'linewidth',1.5);
hold on
% x, y, e, color, alpha
dsErrorsurface(x, y1, e1, [0.4 0.4 0.4], 0.6);
h2_ = plot(x, y2, 'ko--');
dsErrorsurface(x, y2, e2, [0.6 0.6 0.6], 0.6);
set(h2_,'linewidth',1.5);
hold off
% shape axes 
pbaspect([2 1 1]) 
axis square 
ylabel('% \Delta MR activity')
xlabel('time [s]')
title('plot')
% now for the legend... you can specify what you want for each 
legend([h1_, h2_],'first plot', 'second plot')

figure(2)
% asymmetric errorbars - if needed when using, e.g. 5/95 percentiles
eAsym = [0.1+rand(size(x))./2;
    0.5+rand(size(x))./2];
% plot same data as in Figure 1, but with asym errorbars
h1asym_ = plot(x, y1, 'ko-');
set(h1_,'markerfacecolor',[0 0 0]);
hold on
% x, y, e, color, alpha
dsErrorsurface(x, y1, eAsym, [0.6 0.6 0.6], 1 );
title('asymmetric errorbars - like matlab "errorbar" ')
hold off
