function cuts=ea_add_overlay(boundboxmm,cuts,tracor,options)
% This function overlays atlas data over 2d-slice views to be exported by
% LEAD-DBS. The function is called from ea_writeplanes.m
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

    set(0,'CurrentFigure',cuts)
    try
    set(cuts,'GraphicsSmoothing','on')
    end
    % load/generate atlas_index.mat
    if ~isfield(options,'atlases') % atlases structure can be handed down directly within options struct.
        if ~exist([ea_space(options,'atlases'),options.atlasset,filesep,'atlas_index.mat'],'file')
            atlases=ea_genatlastable([],ea_space(options,'atlases'),options);
        else
            load([ea_space(options,'atlases'),options.atlasset,filesep,'atlas_index.mat']);
            atlases=ea_genatlastable(atlases,ea_space(options,'atlases'),options);
        end
    else
        atlases=options.atlases;
    end



    for atlas=1:length(atlases.names)
        for side=detsides(atlases.types(atlas))
            planemm=[boundboxmm{1}(:),boundboxmm{2}(:),boundboxmm{3}(:)];
            %planemm=round(planemm);


            try
                thresh=ea_detthresh(atlases,atlas,atlases.XYZ{atlas,side}.val);
            catch % fibertracts
                thresh=0.5;
            end
            atlmm=atlases.XYZ{atlas,side}.mm(:,1:3);
            %ratlmm=round(atlmm);
            try
                atlvx=atlases.XYZ{atlas,side}.vx;
                atlval=atlases.XYZ{atlas,side}.val;
            catch % fibertracts
                atlvx=round(atlmm);
                for dim=1:3
                atlvx(:,dim)=atlvx(:,dim)-min(atlvx(:,dim));
                end
                atlval=ones(size(atlmm,1),1);
            end

               atlhts=atlmm(:,ea_intersecdim(tracor));
               planehts=planemm(:,ea_intersecdim(tracor));
               dists=abs(atlhts-planehts(1));
               try
               dists=dists<abs(atlases.XYZ{atlas,side}.dims(ea_intersecdim(tracor)))*3;
               catch % fibertracts
                   dists=dists<1*2.5;
               end
            if any(dists) % only if intersection exists plot the atlas.

                xyatl=atlvx(dists,ea_planesdim(tracor));
                xyatlmm=atlmm(dists,ea_planesdim(tracor));
                valatl=atlval(dists);

                for d=1:2
                    [~,minix]=min(xyatl(:,d));
                    [~,maxix]=max(xyatl(:,d));
                    try
                        atlbb(d,:)=[xyatlmm(minix,d);xyatlmm(maxix,d)];
                    catch
                        keyboard
                    end
                end

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
                        if ischar(options.colormap)
                            jetlist = eval(options.colormap);
                        else
                            jetlist = options.colormap;
                        end
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
                        cof = image(colorf);

                        set(cof,'XData',linspace(atlbb(1,1),atlbb(1,2)),'YData',linspace(atlbb(2,1),atlbb(2,2)));
                        slice=slice./max(slice(:));
                        slice=slice*options.d2.atlasopacity;

                        set(cof, 'AlphaData', slice)
                    end
                    %% contour overlay:
                    if  options.d2.con_overlay
                        if any(slice(:));

                            try
                                bw=bwconncomp(slice);
                            catch % no image toolbox available.
                                bw=ea_conncomp(logical(slice));
                            end
                            for cp=1:length(bw.PixelIdxList)
                                slice(:)=0;
                                slice(bw.PixelIdxList{cp})=1;
                                [c] = contourc(ea_zeroframe(slice),1);


                                c=c(:,2:end);
                                [~,yy]=find(c<1);
                                c(:,yy)=[];

                                   dd=sum(abs(diff(c'))');
                                   ix=[];
                                if any(dd>3) % detect jumps in contour (=inner holes)
                                    ix=find(dd>3);
                                    ix=[ix,size(c,2)];
                                end

                                cscale=c;
                                for dim=1:2
                                    cscale(dim,:)=((c(dim,:)/size(slice',dim))*(atlbb(dim,2)-atlbb(dim,1)))+atlbb(dim,1);

                                end
                                set(0,'CurrentFigure',cuts)
                                if isempty(ix)
                                    warning('off')
                                    plot(cscale(1,:),cscale(2,:),'color',options.d2.con_color);
                                warning('on')
                                else
                                    startPoint=1;
                                    for plots=1:length(ix) % this happens if contour has an "inner hole"
                                    plot(cscale(1,startPoint:ix(plots)),cscale(2,startPoint:ix(plots)),'color',options.d2.con_color);
                                    startPoint=ix(plots)+1;
                                    end
                                end
                            end
                        end
                    end

                    if options.d2.lab_overlay

                        if any(slice(:));
                            centr=mean(xyatlmm(valatl>thresh,:));%ea_centroid(logical(slice));
                            an=ea_underscore2space(atlases.names{atlas}(1:find(atlases.names{atlas}=='.')-1));

                            try
                                set(0,'CurrentFigure',cuts)
                                text(centr(1),centr(2),an,'color',options.d2.con_color,'VerticalAlignment','middle','HorizontalAlignment','center');
                            catch

                            end
                        end
                    end
                end
            end
        end

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


function fslice=ea_zeroframe(slice)

fslice=zeros(size(slice)+[2,2]);
fslice(2:end-1,2:end-1)=slice;
