function [h,myTitle] =  myboxplot2(x, y, stdA, stdE, dx, Col, Width, noWhisker)
%' myboxplot(x, y, dx, Col, Width, noWhisker)' Produces a single box plot.
%   PRE:  x is the x-value (scalar!) over which to plot
%         y = scalar/(column!!!) vector w/ the data points to plot
%         dx = width of the box in default units (of x)
%         Col = color to plot in, see helpMatlab-file...
%         Width = line thickness for plots
%         noWhisker: pass any six'th argument if you do not want 
%                    a whisker to be plotted
%
%   POST: plots a 
%         blue square  = mean
%         blue box     = mean +- standard error := sigma/sqrt(n)
%         blue whisker = mean +- standard deviation = sigma


% Make sure y is a vector.
if min(size(y)) ~= 1, 
    error('Requires a vector first argument.'); 
end

if ~exist('noWhisker','var')
  noWhisker = 0;
end

if nargin <5
    error('Requires five input arguments.  x is the x-value (scalar!) over which to plot,  y = scalar/(coulmn) vector w/ the data points to plot,  dx = width opf the box in default units (of x),  Col = color to plot in, see helpMatlab-file...');
end

dxWhisker = 0.6*dx;
% define the mean instead of the data and the 
% +- sigma = stdDeviation instead

  %%if(max(size(y))>1)  %% if y is really a vector and not a scalar...
    median_y    = y; %= median(y);
    mean_y      = y; %mean(y);
    sigma_y     = stdA; %std( y );
    stdError_y  = stdE; %sigma_y / sqrt(length(y));
  %%end
  %ymin = min(y);
  %ymax = max(y);

  % dx = (ymax-ymin)/20;
  % dx   = sigma_y/20;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot the \sigma interval as a whisker
  %% errorbar(...,LineSpec) draws the error bars using the line type, marker symbol, and color specified by LineSpec.
  %%  h_sigma_yInterval = errorbar(x, mean_y ,sigma_y  );
  %%  set(h_sigma_yInterval,'Color','Col');
  if ~noWhisker %% else don't plot the whisker
    whiskerLine_x  = [x;              x];
    whiskerLine_y  = [mean_y-sigma_y; mean_y+sigma_y];
    h_whisUp   = plot(whiskerLine_x, whiskerLine_y, 'Color',Col, 'LineWidth',Width);

    S_sigma_yInterval  = ['whisker: +-1\sigma interval = stdDev = ', num2str(sigma_y) ];

    hold on;
    whisker_x  = [x-dxWhisker;           x+dxWhisker];
    whisker_y  = [mean_y+sigma_y; mean_y+sigma_y];
    h_whisUp   = plot(whisker_x, whisker_y, 'Color',Col, 'LineWidth',Width);
    
    whisker_x  = [x-dxWhisker;           x+dxWhisker];
    whisker_y  = [mean_y-sigma_y; mean_y-sigma_y];
    h_whisDown = plot(whisker_x, whisker_y, 'Color',Col, 'LineWidth',Width);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot the box for the standard error
  boxXmin = x-dx;
  % boxXmax = x+dx;
  boxYmin = mean_y - stdError_y;
  % boxYmax = mean_y + stdError_y;
  % boxX    = [boxXmin; boxXmax; boxXmax; boxXmin; boxXmin];
  % boxY    = [boxYmin; boxYmin; boxYmax; boxYmax; boxYmin];
  % h_box   = plot(boxX, boxY, 'Color',Col, 'LineWidth',Width);
  % rectangle('Position', [boxXmin boxYmin 2*dx 2*stdError_y],'FaceColor',Col,'EdgeColor','none');
  rectangle('Position', [boxXmin boxYmin 2*dx 2*stdError_y],'FaceColor',Col);
  S_box   = ['box: stdError of mean = \sigma/sqrt(n)'];

  hold on;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot mean and data points
  % h_mean_y   = plot(x, mean_y , '  square', 'Color',Col, 'LineWidth',Width);

  h_mean_y   = plot(x, mean_y , 'w.','MarkerSize',6);
  h_mean_y   = plot(x, mean_y , 'k.','MarkerSize',4);

  S_mean_y   = ['mean: ', num2str(x), num2str(mean_y, 4)];  %% create String for labelling the plot

%  h_median_y   = plot(x, median_y , ' *', 'Color',Col, 'LineWidth',Width);
%  S_median_y   = ['median: ', num2str(x), num2str(median_y, 4)];  %% create String for labelling the plot


%%  h_y        = plot( x*ones(size(y,1)),  y , 'r+' );
%%  S_y        = ['y(s)'];  %% create String for labelling the plot

h=h_whisUp; %h_mean_y;  %% handle for the mean plot (for legend)

%myTitleCell={[S_mean_y];...
%              [S_box]; ...
%              [S_sigma_yInterval]};
%             %% convert it to a padded character array (equal length, rest is zero), see intro, p. 64
%myTitle = char(myTitleCell);
if nargin < 6  %% else don't plot the whisker
  myTitle=[S_mean_y,'  ',S_box,'  ',S_sigma_yInterval];
else
  myTitle=[S_mean_y,'  ',S_box];
end
%  h      = [h_mean_y, h_y, h_box]; %% define a handle for all the graphics
%  legend(h, S_mean_y, S_y, S_box, 0); %%0 = pos);


