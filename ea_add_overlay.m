function cuts=ea_add_overlay(boundbox,cuts,side,tracor,patientname,options)
% This function overlays atlas data over 2d-slice views to be exported by
% eAuto-DBS. The function is called from ea_writeplanes.m
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn
useatlases=~strcmp(options.atlasset,'Use none');

needtodelete=0; % small flag to cleanup files involved in .nii.gz support.
try
IV=spm_vol([options.root,patientname,filesep,options.prefs.gtranii]);
catch
IV=spm_vol([options.root,patientname,filesep,options.prefs.tranii]);
end
if useatlases
    
    if ~exist([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'file')
        atlases=ea_genatlastable(options.earoot,options);
    else
        load([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat']);
    end
    
    
    
    for atlas=1:length(atlases.names)
        
        switch atlases.types(atlas)
            case 1 % left hemispheric atlas.
                Vn=[options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
                
            case 2 % right hemispheric atlas.
                Vn=[options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
            case 3 % both-sides atlas composed of 2 files.
                switch side
                    case 1
                        Vn=[options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep,atlases.names{atlas}];
                        
                    case 2
                        Vn=[options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep,atlases.names{atlas}];
                end
            case 4 % mixed atlas (one file with both sides information but two clusters).
                Vn=[options.earoot,'atlases',filesep,options.atlasset,filesep,'mixed',filesep,atlases.names{atlas}];
            case 5 % midline atlas (one file with both sides information but only one cluster).
                Vn=[options.earoot,'atlases',filesep,options.atlasset,filesep,'midline',filesep,atlases.names{atlas}];
        end
        
        [p,f,ext]=fileparts(Vn);
        switch ext
            case '.gz'
                
                copyfile(fullfile(p,[f,'.gz']),fullfile(p,[f,'_temp.nii.gz']))
                try
                    gunzip(fullfile(p,[f,'_temp.nii.gz']));
                    
                catch
                    
                    system(['gunzip ',fullfile(p,[f,'_temp.nii.gz'])]);
                end
                V=spm_vol(fullfile(p,[f,'_temp.nii']));
                needtodelete=1;
            otherwise
                V=spm_vol(Vn);
        end
        
        
        xyzvox=[boundbox{1}(:),boundbox{2}(:),boundbox{3}(:),ones(size(boundbox{1},2),1)]';
        xyzmm=IV.mat*xyzvox;
        xyzatlvox=V.mat\xyzmm;
        switch tracor
            case 1 % transversal images
                [cmesh.X,cmesh.Y]=meshgrid(xyzatlvox(1,:),xyzatlvox(2,:));
                cmesh.Z=repmat(xyzatlvox(3,1),numel(cmesh.X),1);
                ima=spm_sample_vol(V,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),-1);
            case 2 % coronar images
                [cmesh.X,cmesh.Z]=meshgrid(xyzatlvox(1,:),xyzatlvox(3,:));
                cmesh.Y=repmat(xyzatlvox(2,1),numel(cmesh.X),1);
                ima=spm_sample_vol(V,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),-1);
            case 3 % saggital images
                [cmesh.Y,cmesh.Z]=meshgrid(xyzatlvox(2,:),xyzatlvox(3,:));
                cmesh.X=repmat(xyzatlvox(1,1),numel(cmesh.Y),1);
                ima=spm_sample_vol(V,cmesh.X(:),cmesh.Y(:),cmesh.Z(:),-1);
        end
        
        slice=reshape(ima,length(xyzvox),length(xyzvox));
        slice=flipud(slice);
        
        
        maxcolor=64;
        try
            jetlist=options.colormap;
        catch
            try
                jetlist=atlases.colormap;
            catch
                try
                    jetlist=eval(atlases.colormap);
                catch
                    jetlist=jet;
                end
            end
        end
        atlases.colormap=jetlist;
        
        try % see if there is explicit color information for this atlas
            atlasc=squeeze(jetlist(round(atlases.colors(atlas)),:)); % color for toggle button icon
        catch
            atlases.colors(atlas)=atlas*(maxcolor/length(atlases.names));
            atlasc=squeeze(jetlist(round(atlases.colors(atlas)),:)); % color for toggle button icon
            
        end
        
        colorf=zeros(length(boundbox{1}),length(boundbox{2}),3);
        colorf(:,:,1)=atlasc(1);
        colorf(:,:,2)=atlasc(2);
        colorf(:,:,3)=atlasc(3);
        
        
        %% color_overlay:
        if options.d2.col_overlay
            cof = imshow(colorf);
            
            
            slice=slice*options.d2.atlasopacity;
            
            set(cof, 'AlphaData', slice*0.3)
        end
        %% isoval overlay:
        if options.d2.con_overlay
            warning('off','all');
            con = contour(slice,1,'Edgecolor',options.d2.con_color);
        end
        
        if options.d2.lab_overlay
            
            centr=ea_regionprops(logical(slice),'Centroid');
            
            an=atlases.names{atlas}(1:find(atlases.names{atlas}=='.')-1);
           
           try
                text(centr.Centroid(1),centr.Centroid(2),an,'color',options.d2.con_color);
           catch
               
           end
        end
        
      
        if needtodelete
            
                    delete(fullfile(p,[f,'_temp.nii']));
                    delete(fullfile(p,[f,'_temp.nii.gz']));

                    needtodelete=0;
        end
        
    end
    
    axis off
    
    % save table information
    save([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'atlases');
end
