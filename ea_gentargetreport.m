function ea_gentargetreport(M)



thresh = inputdlg('Please enter threshold in mm:','Enter threshold...',1,{'0.5'});

thresh=str2double(thresh{1});


for thr=0:1

for target=1:length(M.ui.volumeintersections)
    if thr
        rf(target,thr+1)=figure('color','w','Numbertitle','off','name',['Electrode centers residing in ',M.vilist{M.ui.volumeintersections(target)}]);
    else
rf(target,thr+1)=figure('color','w','Numbertitle','off','name',['Distances of electrode centers to nearest voxel in ',M.vilist{M.ui.volumeintersections(target)}]);
    end
distances=zeros(8,length(M.ui.listselect));
    
    for pt=1:length(M.ui.listselect)
        
        distances(:,pt)=[M.stats(M.ui.listselect(pt)).ea_stats.conmat{1}(:,M.ui.volumeintersections(target));... % right side
            M.stats(M.ui.listselect(pt)).ea_stats.conmat{2}(:,M.ui.volumeintersections(target))];
        
        
        
    end
    
    if thr
        distances=distances<thresh;
    end
        
    
   R{target,thr+1}=distances; 

   cnames=M.patient.list(M.ui.listselect');
   [~,cnames]=cellfun(@fileparts,cnames,'UniformOutput',0);
   rnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
   
   t(target,thr+1)=uitable(rf(target,thr+1),'Data',R{target,thr+1},'ColumnName',cnames,'RowName',rnames);
   
      figdims=get(rf(target,thr+1),'Position');
      textend=get(t(target,thr+1),'Extent');
   set(rf(target,thr+1),'Position',[figdims(1:2),textend(3:4)]);
   figdims=get(rf(target,thr+1),'Position');
   set(t(target,thr+1),'Position',[0,0,figdims(3),figdims(4)])
   
   
  
end
end

assignin('base','R',R);
