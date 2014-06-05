function fib_plots=showfibres(fig,coords,options)
ht=uitoolbar(fig);
    colornames='rbgcmywkrbgcmywkrbgcmywk';

load(fullfile('fibres','gibbsconnectome5.mat'));
cnt=1;
for fib=1:100:length(gibbsconnectome)    
    
    for contact=1:8
        [IDX,D]=rangesearch(coords(contact,:),2,gibbsconnectome{fib},2);
        %plot3(gibbsconnectome{fib}(:,1),gibbsconnectome{fib}(:,2),gibbsconnectome{fib}(:,3),'k-');
        if ~isempty(IDX)
        connectingfibs{cnt}=gibbsconnectome{fib}';
        cnt=cnt+1;
        end
    end
end
if ~isempty(connectingfibs)

    for fib=1:length(connectingfibs)
        %for segment=1:length(connectingfibs{fib})-1;

        connectingfibs{fib}(4,:)=detcolor(connectingfibs{fib}); % add coloring information to the 4th column.
       
        for dim=1:4
        thisfib(dim,:)=double(interp1q([1:size(connectingfibs{fib},2)]',connectingfibs{fib}(dim,:)',[1:0.1:size(connectingfibs{fib},2)]')');
        end
        fib_plots(fib)=surface([thisfib(1,:);thisfib(1,:)],...
            [thisfib(2,:);thisfib(2,:)],...
            [thisfib(3,:);thisfib(3,:)],...
            [thisfib(4,:);thisfib(4,:)],'facecol','no','edgecol','interp','linew',1.5);
        clear thisfib
        
        %plot3(connectingfibs{fib}(segment:segment+1,1),connectingfibs{fib}(segment:segment+1,2),connectingfibs{fib}(segment:segment+1,3),'Color',segclr);

        %end
    end
    
    set(fib_plots(:),'EdgeAlpha',0.05);
    
    
    
    
    fiberbutton=uitoggletool(ht,'CData',get_icn('fibers'),'TooltipString','Fibers','OnCallback',{@fibresvisible,fib_plots},'OffCallback',{@fibresinvisible,fib_plots},'State','on');

end




if options.showconnectivities
    
    %% extract areas connected by fibres.
   atlas=load_nii(fullfile('templates','aal','aal.nii'));
   V=spm_vol(fullfile('templates','aal','aal.nii'));
   aID = fopen(fullfile('templates','aal','aal.txt'));
   atlas_lgnd=textscan(aID,'%d %s');
   allcareas=[];
   for fib=1:length(connectingfibs)
       
       thisfibendpoints=[connectingfibs{fib}(1:3,1),connectingfibs{fib}(1:3,end)];
       [dummy,thisfibendpoints]=map_coords(thisfibendpoints,fullfile('templates','aal','aal.nii'));
       thisfibendpoints=double(thisfibendpoints); clear dummy;
       conareas=spm_sample_vol(V,thisfibendpoints(1,:),thisfibendpoints(2,:),thisfibendpoints(3,:),0);
      allcareas=[allcareas,conareas]; 
   end
   allcareas=unique(allcareas(allcareas>0));
   
   %% now show areas
   atlas.img=round(atlas.img);
   for anatarea=1:length(allcareas)
       [xx,yy,zz]=ind2sub(size(atlas.img),find(atlas.img==allcareas(anatarea)));
       XYZ=[xx,yy,zz];
       XYZ=map_coords(XYZ',fullfile('templates','aal','aal.nii'));
       XYZ=XYZ';
       k=convhulln(XYZ);
       [dummy,centroid]=kmeans(XYZ,1);
       centroid=centroid(1,:); clear dummy;
       regionsurfs(anatarea)=trisurf(k,XYZ(:,1),XYZ(:,2),XYZ(:,3),...
                abs(repmat(anatarea*(64/length(allcareas)),length(XYZ),1)...
                +randn(length(XYZ),1)*0.1*length(allcareas))');
            
            %% shading etc.
        colorc=colornames(anatarea);
        colorc=rgb(colorc);
        spec_atlas(regionsurfs(anatarea),colorc,'aal',1);
        
        
        %% put a label to it
        thislabel=atlas_lgnd{2}(atlas_lgnd{1}==allcareas(anatarea));
        labels(anatarea)=text(centroid(1),centroid(2),centroid(3),thislabel);
   end
   
       regionbutton=uitoggletool(ht,'CData',get_icn('regions'),'TooltipString','Fibers','OnCallback',{@regionsvisible,regionsurfs},'OffCallback',{@regionsinvisible,regionsurfs},'State','on');
       captionbutton=uitoggletool(ht,'CData',get_icn('captions'),'TooltipString','Fibers','OnCallback',{@regionsvisible,labels},'OffCallback',{@regionsinvisible,labels},'State','on');
   
   
end


function indcol=detcolor(mat) % determine color based on traversing direction.

xyz=abs(diff(mat,1,2));
rgb=xyz/max(xyz(:));
rgb=[rgb,rgb(:,end)];
rgbim=zeros(1,size(rgb,2),3);
rgbim(1,:,:)=rgb';
indcol=double(rgb2ind(rgbim,jet));

function fibresvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);

function fibresinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');



function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);


function regionsvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);

function regionsinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);

