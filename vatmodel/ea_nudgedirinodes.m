function nodes=ea_nudgedirinodes(nodes,centroid)
% get to ~1000 comparison points
div=round(length(nodes)/1000);
if div<1
    div=1;
end
nodes=nodes(1:div:end,:);
nodes=nodes.*9;
nodes=nodes+repmat(centroid,length(nodes),1);
nodes=nodes./10;