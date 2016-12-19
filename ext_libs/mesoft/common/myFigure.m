function  h = myFigure()
%      'h = myFigure()' 
%       creates a standard plottable figure
% PRE:  
%        
%      

%%  rect = [0.25, 0.25, 11.193, 7.7677];  %%% for normalized PaperUnits
  rect = [0.1, 0.1, 11.59, 8.16];  %%% 0nly 2.54 mm border around figure! --> for psnup! for normalized PaperUnits
  h=figure('PaperOrientation','landscape',  ...
           'PaperType'       ,'a4', ...
           'PaperUnits'      ,'inches', ...  %% default value.
           'PaperPosition'   ,rect);             %%  rect = [left, bottom, width, height]

% a4: width = 210, height = 297 mm
% ==       8.2677 * 11.692 inches
