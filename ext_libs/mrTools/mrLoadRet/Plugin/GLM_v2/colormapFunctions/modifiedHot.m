function cmap = modifiedHot(numColors)

cmap = hot(round(numColors*312/256));
cmap = cmap(end-numColors+1:end,:);

  
