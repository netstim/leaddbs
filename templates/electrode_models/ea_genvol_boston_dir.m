function ea_genvol_boston_dir(meshel,elspec,z3ix,vizz)

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

% clean from duplicate indices:

% for ff=1:length(fcyl)
%     [has,which]=ea_hasduplicates(fcyl{ff});
%     if has
%         doubles=find(fcyl{ff}==which);
%         fcyl{ff}(doubles(2:end))=[];
%     end
% end


%% convert to obtain the electrode surface mesh model


nodecon=unique_ncyl;
facecon=unique_fcyl;


%[nodecon,~,facecon]=s2m(nodecon,facecon,electrodetrisize,100,'tetgen',seeds,[]); % generate a tetrahedral mesh of the cylinders


if vizz
    figure
    patch('faces',facecon(:,1:3),'vertices',nodecon,'edgecolor','k','facecolor','none');

end






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
%
% for ff=1:length(fcyl)
%     [has,which]=ea_hasduplicates(fcyl{ff});
%     if has
%         doubles=find(fcyl{ff}==which);
%         fcyl{ff}(doubles(2:end))=[];
%     end
% end

keyboard
nodeins=unique_ncyl;
faceins=unique_fcyl;

if vizz
    figure
    patch('faces',faceins(:,1:3),'vertices',nodeins,'edgecolor','k','facecolor','none');

end

%[ncyl,fcyl]=meshcheckrepair(unique_ncyl,fcyl,'deep');
ISO2MESH_SURFBOOLEAN='cork';   % now intersect the electrode to the nucleus
[outer_node,outer_face]=surfboolean(nodecon,facecon(:,1:3),'all',nodeins,faceins);

convhull=convhulln(outer_node);

[ins_node,ins_face]=surfboolean(outer_node,convhull,'diff',nodecon,facecon(:,1:3));

[comb_node,comb_face]=surfboolean(ins_node,ins_face,'xor',nodecon,facecon(:,1:3));


clear ISO2MESH_SURFBOOLEAN;
if vizz
    figure

    patch('faces',comb_face(:,1:3),'vertices',comb_node(:,1:3),'edgecolor','k','facecolor','none');

    figure

    patch('faces',outer_face(:,1:3),'vertices',outer_node(:,1:3),'edgecolor','k','facecolor','none');

  figure

    patch('faces',ins_face(:,1:3),'vertices',ins_node(:,1:3),'edgecolor','k','facecolor','none');

end


%comb_face=unique(comb_face,'rows');

keyboard
[ncomb_node,~,ncomb_face]=s2m(comb_node,comb_face,electrodetrisize,100,'tetgen',seeds,[]); % generate a tetrahedral mesh of the cylinders












save([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'],'node','face');


