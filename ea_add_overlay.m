function cuts=ea_add_overlay(boundboxmm,cuts,side,tracor,patientname,options)
% This function overlays atlas data over 2d-slice views to be exported by
% LEAD-DBS. The function is called from ea_writeplanes.m
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
    % load/generate atlas_index.mat
    
    if ~exist([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'file')
        atlases=ea_genatlastable([],options.earoot,options);
    else
        load([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat']);
        atlases=ea_genatlastable(atlases,options.earoot,options);
    end
    
    
    
    for atlas=1:length(atlases.names)
        for side=options.sides
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

            
            
            planemm=[boundboxmm{1}(:),boundboxmm{2}(:),boundboxmm{3}(:)];
            planemm=round(planemm);
            
            thresh=ea_detthresh(atlases,atlas,atlases.XYZ{atlas,side}.val)*1.5;
            atlmm=atlases.XYZ{atlas,side}.mm;
            ratlmm=round(atlmm);
            
            atlvx=atlases.XYZ{atlas,side}.vx;
            atlval=atlases.XYZ{atlas,side}.val;
            
            for d=ea_planesdim(tracor)
                [~,minix]=min(atlvx(:,d));
                [~,maxix]=max(atlvx(:,d));
                atlbb(d,:)=[atlmm(minix,d);atlmm(maxix,d)];
                
            end
            
                    
            if any(ratlmm(:,ea_intersecdim(tracor))==planemm(1,ea_intersecdim(tracor))) % only if intersection exists plot the atlas.
                keyboard
                xyatl=atlvx(ratlmm(:,ea_intersecdim(tracor))==planemm(1,ea_intersecdim(tracor)),ea_planesdim(tracor));
                valatl=atlval(ratlmm(:,ea_intersecdim(tracor))==planemm(1,ea_intersecdim(tracor)));
                
                
                
                
                xyatl(:,1)=xyatl(:,1)-min(xyatl(:,1))+1;
                xyatl(:,2)=xyatl(:,2)-min(xyatl(:,2))+1;
                slice=zeros(max(xyatl(:,1)),max(xyatl(:,2)));
                slice(sub2ind(size(slice),xyatl(:,1),xyatl(:,2)))=valatl;
                slice=interp2(slice,3);
                slice(slice<thresh)=0;
                slice=slice';
                
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
                
                colorf=zeros(size(slice,1),size(slice,2),3);
                colorf(:,:,1)=atlasc(1);
                colorf(:,:,2)=atlasc(2);
                colorf(:,:,3)=atlasc(3);
                
                
                %% color_overlay:
                if options.d2.col_overlay
                    cof = imshow(colorf);

                            set(cof,'XData',linspace(atlbb(1,1),atlbb(1,2)),'YData',linspace(atlbb(2,1),atlbb(2,2)));
                    slice=slice*options.d2.atlasopacity;
                    
                    set(cof, 'AlphaData', slice*0.3)
                    
                    
                end
                %% isoval overlay:
                if  options.d2.con_overlay
                    if any(slice(:));
                        warning('off','all');
                        
                        %slice=logical(slice);
                        
                        %[~,con] = contour(slice,1);
                        
                        bw=bwconncomp(slice);
                        
                        for cp=1:length(bw.PixelIdxList)
                            slice(:)=0;
                            slice(bw.PixelIdxList{cp})=1;
                            [c] = contourc(slice,1);
                            
                            [~,yy]=find(~floor(c));
                            c(:,yy)=[];
                            
                            
                            cscale=c;
                            for dim=1:2
                                cscale(dim,:)=((c(dim,:)/size(slice',dim))*(atlbb(dim,2)-atlbb(dim,1)))+atlbb(dim,1);

                            end
                            plot(cscale(1,:),cscale(2,:),'color',options.d2.con_color,'LineSmoothing','on');
                        end
                    end
                    %set(con,'XData',linspace(atlbb(1,2),atlbb(1,1),20),'YData',linspace(atlbb(2,2),atlbb(2,1),20)); %,'LineWidth',1,'Edgecolor',options.d2.con_color);
                    
                end
                
                
                if options.d2.lab_overlay
                    
                    
                    
                    centr=mean(atlmm(atlmm(:,2)<0,:));%ea_centroid(logical(slice));
                    an=atlases.names{atlas}(1:find(atlases.names{atlas}=='.')-1);
                    
                    try
                        text(centr(1),centr(2),an,'color',options.d2.con_color);
                    catch
                        
                    end
                    
                end
                
            end
        end
        
    end
    
    
    
    
    %axis off
    
    % save table information
    save([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'atlases');
end



function in=ea_intersecdim(tracor)

switch tracor
    case 1
        in=3;
    case 2
        in=2;
    case 3
        in=1;
end

function pl=ea_planesdim(tracor)

switch tracor
    case 1
        pl=[1,2];
    case 2
        pl=[1,3];
    case 3
        pl=[2,3];
end