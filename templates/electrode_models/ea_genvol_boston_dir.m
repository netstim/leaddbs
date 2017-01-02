function ea_genvol_boston_dir(meshel,elspec,vizz)

electrodetrisize=0.1;  % the maximum triangle size of the electrode mesh

%







%% loading the electrode surface model

ncyl=[];
fcyl=[];
scyl=[];
seeds=[];


for i=1:length(meshel.con)
    fcyl=[fcyl; meshel.con{i}.faces+size(ncyl,1)];
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

% clean from duplicate indices:

for ff=1:length(fcyl)
    [has,which]=ea_hasduplicates(fcyl{ff});
    if has
        doubles=find(fcyl{ff}==which);
        fcyl{ff}(doubles(2:end))=[];
    end
end


%% convert to obtain the electrode surface mesh model

[node,~,face]=s2m(unique_ncyl,fcyl,electrodetrisize,100,'tetgen',seeds,[]); % generate a tetrahedral mesh of the cylinders

keyboard







% now add insulation

ncyl=[];
fcyl=[];
scyl=[];
seeds=[];

for i=1:length(meshel.ins)
    fcyl=[fcyl; meshel.ins{i}.faces+size(ncyl,1)];
    
    ncyl=[ncyl; meshel.ins{i}.vertices];
    seeds=[seeds; mean(meshel.ins{i}.vertices)];
end

[ncyl,fcyl]=meshcheckrepair(ncyl,fcyl,'dup');
[ncyl,fcyl]=meshcheckrepair(ncyl,fcyl,'deep');
if vizz
   figure
   patch('vertices',ncyl,'faces',fcyl,'facecolor','none')
    axis equal
    zlim([0,10]);
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

% clean from duplicate indices:

for ff=1:length(fcyl)
    [has,which]=ea_hasduplicates(fcyl{ff});
    if has
        doubles=find(fcyl{ff}==which);
        fcyl{ff}(doubles(2:end))=[];
    end
end
[ncyl,fcyl]=meshcheckrepair(ncyl,fcyl,'deep');



[nodeins,~,faceins]=s2m(unique_ncyl,fcyl,electrodetrisize,100,'tetgen',seeds,[]); % generate a tetrahedral mesh of the cylinders












save([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'],'node','face');


