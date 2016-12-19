function h = myFigPortrait()
%      'h = myFigPortrait()' 
%       creates a standard plottable figure
% PRE:  
%        
%      
Position = [10 180  420 560];

%%  rect = [0.25, 0.25, 11.193, 7.7677];  %%% for normalized PaperUnits
% rect = [0.2, 0.2, 8.06, 11.49];  %%% 0nly 2.54 mm border around figure! --> for psnup! for normalized PaperUnits
  rect = [0.1, 0.1, 8.16, 11.59];  %%% 0nly 2.54 mm border around figure! --> for psnup! for normalized PaperUnits
  h=figure('PaperOrientation','portrait',  ...
           'PaperType'       ,'a4', ...
           'PaperUnits'      ,'inches', ...  %% default value.
           'PaperPosition'   ,rect, ...             %%  rect = [left, bottom, width, height]
	   'Position', Position);            

% a4: width = 210, height = 297 mm
% ==       8.2677 * 11.692 inches
