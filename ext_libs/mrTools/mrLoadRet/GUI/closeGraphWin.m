function closeGraphWin(h)
%
% closeGraphWin([h])
%
%  $Id$
%
% Closes the current graphWin, or graphWin with handle h, that's part of MLR graph figures.  Sets global variable MLR
%
% djh, 3/3/98

mrGlobals

if ~isempty(MLR.graphFigure)
  if ieNotDefined('h')
    h = selectGraphWin(1);
  end
  % save the position
  if ismember(h,MLR.graphFigure)
    mrSetFigLoc(['graphFigure' int2str(mlrGetFignum(h))],get(h,'Position'));
    MLR.graphFigure = setdiff(MLR.graphFigure,h); 
    delete(h)
    return;
  end
end

% Don't do that unless the handle to the main GUI handle is made invisible:
 
% %if there are no graph figures in MLR, just delete the figure
% if isempty(h)
%   h = gcf;               
% end
%   delete(h)
% end


