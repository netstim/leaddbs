function ea_genvol_abbott(meshel,elspec,vizz)

electrodetrisize=0.1;  % the maximum triangle size of the electrode mesh

%% loading the electrode surface model
ncyl=[];
fcyl=[];
scyl=[];
seeds=[];

for i=1:length(meshel.ins)
    fcyl=[fcyl; meshel.ins{i}.faces+size(ncyl,1)];

        scyl=[scyl; meshel.ins{i}.endplates+size(ncyl,1)]; % had to rebuild the endplates
    ncyl=[ncyl; meshel.ins{i}.vertices];
    seeds=[seeds; mean(meshel.ins{i}.vertices)];
end
for i=1:length(meshel.con)
    fcyl=[fcyl; meshel.con{i}.faces+size(ncyl,1)];
    if i>1 % not for active tip
        scyl=[scyl; meshel.con{i}.endplates+size(ncyl,1)];
    end
    ncyl=[ncyl; meshel.con{i}.vertices];
    seeds=[seeds; mean(meshel.con{i}.vertices)];
end

[unique_ncyl, I, J]=unique(ncyl, 'rows');
unique_fcyl=unique(round(J(fcyl)),'rows');
unique_scyl=unique(round(J(scyl)),'rows');
if vizz
    figure
    patch('faces',fcyl,'vertices',ncyl,'edgecolor','k','facecolor','none');

    figure
    patch('faces',unique_fcyl,'vertices',unique_ncyl,'edgecolor','b','facecolor','none');
end
fcyl=num2cell(unique_fcyl,2);
scyl=num2cell(unique_scyl,2);

% clean from duplicate indices:
for ff=1:length(fcyl)
    [has,which]=ea_hasduplicates(fcyl{ff});
    if has
        doubles=find(fcyl{ff}==which);
        fcyl{ff}(doubles(2:end))=[];
    end
end

for ff=1:length(scyl)
    [has,which]=ea_hasduplicates(scyl{ff});
    if has
        doubles=find(scyl{ff}==which);
        scyl{ff}(doubles(2:end))=[];
    end
end

%% convert to obtain the electrode surface mesh model
[node,~,face]=s2m(unique_ncyl,{fcyl{:}, scyl{:}},electrodetrisize,100,'tetgen',seeds,[]); % generate a tetrahedral mesh of the cylinders
save([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'],'node','face');
