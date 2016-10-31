%mrLoadRet color scale function, for use in the overlay field 'userDefinedCMap'
%creates a color map that is dark red values above a threshold and goes to yellow for values below
%the threshold is set to threshold % of the scale, but can be changed here
% jb, 13/10/2010
%        $Id$

%not in use right now

function colormap = statsColorMap(numberColors)

threshold = .05;
% colormap is made with a little bit less on the dark end
belowThreshold = floor(numberColors*threshold);
hotRange = ceil(numberColors*threshold);

reducedHot = flipud(hot(hotRange));
colormap = repmat(reducedHot(belowThreshold+1,:),numberColors,1);
colormap(1:belowThreshold,:) = reducedHot(1:belowThreshold,:);
