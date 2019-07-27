function [i] = ea_searchclosest(x,v)
z=abs(x-v); 
[~,ix]=sort(z,'ascend');
i=ix(1);