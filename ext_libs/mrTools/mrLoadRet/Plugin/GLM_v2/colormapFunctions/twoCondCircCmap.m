function cmap = twoCondCircCmap(numGrays)
  if isodd(numGrays)
    numGrays = numGrays +1;
  end
  cmapA = twoCondCmap(numGrays/2);
  cmapB = flipud(cmapA);
  cmap = cat(1,cmapA, cmapB);
return
  
