
load data
h=figure;
     plotmesh(allnode,allface,'linestyle','-','facealpha',0.1);

[node,~,face]=s2m(allnode,{allface{:}},electrodetrisize,100,'tetgen',[],[]); % generate a tetrahedral mesh of the cylinders
