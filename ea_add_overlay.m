function cuts=ea_add_overlay(boundboxmm,cuts,tracor,options)
% This function overlays atlas data over 2d-slice views to be exported by
% LEAD-DBS. The function is called from ea_writeplanes.m
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

useatlases=~strcmp(options.atlasset,'Use none');
    set(0,'CurrentFigure',cuts)
needtodelete=0; % small flag to cleanup files involved in .nii.gz support.
try
    IV=spm_vol([options.root,options.patientname,filesep,options.prefs.gtranii]);
catch
    IV=spm_vol([options.root,options.patientname,filesep,options.prefs.tranii]);
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
        for side=detsides(atlases.types(atlas))
            planemm=[boundboxmm{1}(:),boundboxmm{2}(:),boundboxmm{3}(:)];
            %planemm=round(planemm);
            
            thresh=ea_detthresh(atlases,atlas,atlases.XYZ{atlas,side}.val);

            atlmm=atlases.XYZ{atlas,side}.mm;
            %ratlmm=round(atlmm);
            
            atlvx=atlases.XYZ{atlas,side}.vx;
            atlval=atlases.XYZ{atlas,side}.val;
            pdcnt=1;
            for d=ea_planesdim(tracor)
                [~,minix]=min(atlvx(:,d));
                [~,maxix]=max(atlvx(:,d));
                atlbb(pdcnt,:)=[atlmm(minix,d);atlmm(maxix,d)];
                pdcnt=pdcnt+1;                
            end
            
                  
               
               atlhts=atlmm(:,ea_intersecdim(tracor));
               planehts=planemm(:,ea_intersecdim(tracor));
               dists=abs(atlhts-planehts(1));
               
               dists=dists<abs(atlases.XYZ{atlas,side}.dims(ea_intersecdim(tracor)))/2;
            if any(dists) % only if intersection exists plot the atlas.
                
                xyatl=atlvx(dists,ea_planesdim(tracor));
                valatl=atlval(dists);
                
                
                
                
                xyatl(:,1)=xyatl(:,1)-min(xyatl(:,1))+1;
                xyatl(:,2)=xyatl(:,2)-min(xyatl(:,2))+1;
                slice=zeros(max(xyatl(:,1)),max(xyatl(:,2)));
                slice(sub2ind(size(slice),xyatl(:,1),xyatl(:,2)))=valatl;
                if ~any(size(slice)==1) % exception for problems with onedimensional slices
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
                        set(0,'CurrentFigure',cuts)
                        cof = imshow(colorf);
                        
                        
                        set(cof,'XData',linspace(atlbb(1,1),atlbb(1,2)),'YData',linspace(atlbb(2,1),atlbb(2,2)));
                        slice=slice*options.d2.atlasopacity;
                        
                        set(cof, 'AlphaData', slice)
                        
                        
                    end
                    %% isoval overlay:
                    if  options.d2.con_overlay
                        if any(slice(:));
                            bw=bwconncomp(slice);
                            for cp=1:length(bw.PixelIdxList)
                                slice(:)=0;
                                slice(bw.PixelIdxList{cp})=1;
                                [c] = contourc(slice,1);
                                c=c(:,2:end);
                                [~,yy]=find(c<1);
                                c(:,yy)=[];
                                
                                
                                cscale=c;
                                for dim=1:2
                                    cscale(dim,:)=((c(dim,:)/size(slice',dim))*(atlbb(dim,2)-atlbb(dim,1)))+atlbb(dim,1);
                                    
                                end
                                set(0,'CurrentFigure',cuts)
                                plot(cscale(1,:),cscale(2,:),'color',options.d2.con_color,'LineSmoothing','on');
                            end
                        end
                    end
                    
                    
                    if options.d2.lab_overlay
                        
                        
                        
                        centr=mean(atlmm(:,ea_planesdim(tracor)));%ea_centroid(logical(slice));
                        an=atlases.names{atlas}(1:find(atlases.names{atlas}=='.')-1);
                        
                        try
                            set(0,'CurrentFigure',cuts)
                            text(centr(1),centr(2),an,'color',options.d2.con_color);
                        catch
                            
                        end
                        
                    end
                end
            end
        end
        
    end
    
    
    
    
    %axis off
    
    % save table information
    save([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'atlases');
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
    case 5
        sides=1; % midline
        
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