function cmap = hotColdCmap(n)


  h = hot(floor(n/2));
  c(:,1) = h(:,3);
  c(:,2) = h(:,2);
  c(:,3) = h(:,1);
  
  if iseven(h)
    cmap = [flipud(c);h];
  else
    cmap = [flipud(c);[0 0 0];h];
  end

return

