function [atlassurfs]=showatlas(fig)
atlascnt=1;
set(0,'CurrentFigure',fig)

atlases=genatlastable();




setinterpol=1;

ht=uitoolbar(fig);

for atlas=1:length(atlases.names)
    
    switch atlases.types(atlas)
        case 1 % left hemispheric atlas.
            nii=load_nii(['atlases',filesep,'lh',filesep,atlases.names{atlas}]);
            nii.img=round(nii.img);
            
        case 2 % right hemispheric atlas.
            nii=load_nii(['atlases',filesep,'rh',filesep,atlases.names{atlas}]);
            nii.img=round(nii.img);
        case 3 % both-sides atlas composed of 2 files.
            lnii=load_nii(['atlases',filesep,'lh',filesep,atlases.names{atlas}]);
            rnii=load_nii(['atlases',filesep,'rh',filesep,atlases.names{atlas}]);
            lnii.img=round(lnii.img);
            rnii.img=round(rnii.img);
        case 4 % mixed atlas (one file with both sides information.
            nii=load_nii(['atlases',filesep,'mixed',filesep,atlases.names{atlas}]);
            nii.img=round(nii.img);
    end
    
    for sides=detsides(atlases.types(atlas));
        if atlases.types(atlas)==3 % both-sides atlas composed of 2 files.
            if sides==1
                nii=lnii;
            elseif sides==2
                nii=rnii;
            end
        end
        
        colornames='rbgcmywkrbgcmywkrbgcmywk';
        for color=1:max(nii.img(:))
            
            colorc=colornames(color);
            colorc=rgb(colorc);
            
            [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img==color)); % find 3D-points that have correct value.
            if ~isempty(xx)
                disp(['color= ',num2str(color),', length= ',num2str(length(xx)),'.']);
                
                XYZ=[xx,yy,zz]; % concatenate points to one matrix.
                
                XYZ=map_coords_proxy(XYZ',nii.fileprefix)'; % map to mm-space
                
                
                if atlases.types(atlas)==4 % mixed atlas, divide
                    if sides==1
                        XYZ=XYZ(XYZ(:,1)<0,:,:);
                    elseif sides==2
                        XYZ=XYZ(XYZ(:,1)>0,:,:);
                    end
                end
                %surface(xx(1:10)',yy(1:10)',zz(1:10)',ones(10,1)');
                hold on
                
                %trisurf(delaunay3(xx,yy,zz),xx,yy,zz,repmat(color,length(xx),1),'edgecolor',[0.8 0.4 0]);
                
                
                k=convhulln(XYZ);
                
                %CMatrix=ones(length(XYZ));          %(ones(length(XYZ),length(XYZ),3));
                atlassurfs(atlascnt)=trisurf(k,XYZ(:,1),XYZ(:,2),XYZ(:,3),...
                    abs(repmat(atlas*(64/length(atlases.names)),length(XYZ),1)...
                    +randn(length(XYZ),1)*(length(atlases.names)/1))');
                caxis([1 64]);
                
                
                jetlist=jet;
                atlasc=squeeze(jetlist(round(atlas*(64/length(atlases.names))),:));
                
                
                
                colorbuttons(atlascnt)=uitoggletool(ht,'CData',get_icn('atlas',atlasc),'TooltipString',atlases.names{atlas},'OnCallback',{@atlasvisible,atlassurfs(atlascnt)},'OffCallback',{@atlasinvisible,atlassurfs(atlascnt)},'State','on');
                
                
                
                
                spec_atlas(atlassurfs(atlascnt),colorc,atlases.names{atlas},setinterpol);
                
                atlascnt=atlascnt+1;
                
                
                drawnow
                
            end
        end
    end
end




end






function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);
end

function atlasvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);
end
function atlasinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);
end
function sides=detsides(opt)

switch opt
    case 1 % left hemispheric atlas
        sides=1;
    case 2 % right hemispheric atlas
        sides=2;
    case 3
        sides=1:2;
    case 4
        sides=1:2;
        
end
end



function coords=map_coords_proxy(XYZ,fname)

if exist([fname,'.nii'],'file')

    coords=map_coords(XYZ,[fname,'.nii']);
elseif exist([fname,'.nii.gz'],'file')
    copyfile([fname,'.nii.gz'],[fname,'_temp.nii.gz'])
    system(['gunzip ',fname,'_temp.nii.gz']);
    coords=map_coords(XYZ,[fname,'_temp.nii']);
    delete([fname,'_temp.nii']);
end
end

