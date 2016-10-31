function cmap = grayCirc(numGrays)
  if ~iseven(numGrays);
    numGrays = numGrays + 1;
  end
  
  cmapA = gray(numGrays/2 + 1);
  cmapB = flipud(cmapA);
  
  cmap = cat(1,cmapA(1:numGrays/2,:), cmapB(1:numGrays/2,:));
