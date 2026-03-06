function ea_unified_draw_volumetric(obj,vals,group,side,gradientLevel,voxcmap)
if strcmp(obj.drawTool,'sweetspotmapping')
    % Plot voxels if any survived
    res=obj.results.sweetspotmapping.space{side};
    res.dt(1) = 16;
    res.img(:)=nan;
elseif strcmp(obj.drawTool,'networkmapping')
    if contains(obj.calcsettings.netmap_connectome, ' > ')
        % For functional connectome, should use the spacedef provided by the connectome itself
        res = obj.results.networkmapping.(ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome)).space;
    else
        switch obj.calcsettings.structuralresolution  % or functionalresolution if relevant
            case '2 mm'
                outspace = '222';
            case '1 mm'
                outspace = '111';
            case '0.5 mm'
                outspace = '555';
            otherwise
                error('Unsupported resolution: %s', obj.calcsettings.structuralresolution);
        end
        res = ea_load_nii([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,outspace,'.nii.gz']);
    end
    res.dt(1) = 16;


end
%plot based on vizmode
if strcmp(obj.drawTool,'sweetspotmapping') || strcmp(obj.vizmode,'Regions')
    if strcmp(obj.drawTool,'sweetspotmapping')
        subfield = 'sweetspotmapping';
    else
        subfield = 'networkmapping';
    end

    % Plot voxels if any survived
    if obj.posvisible
        % plot positives:
        posvox=res;
        posvox.img(:)=0;
        posvox.img(vals{group,side}>0)=vals{group,side}(vals{group,side}>0);

        pobj.nii=posvox;
        pobj.name='Positive';
        pobj.niftiFilename='Positive.nii';
        pobj.binary=0;
        pobj.usesolidcolor=0;
        pobj.color=obj.poscolor;
        pobj.colormap=ea_colorgradient(gradientLevel, obj.posBaseColor, obj.poscolor);
        pobj.smooth=10;
        pobj.hullsimplify=0.5;
        pobj.threshold=0;
        obj.drawobject.(subfield){group,side}{1}=ea_roi('Positive.nii',pobj);

        res=posvox; % keep copy for export
    end

    if obj.negvisible
        % plot negatives:
        negvox=res;
        negvox.img(:)=0;
        negvox.img(vals{group,side}<0)=-vals{group,side}(vals{group,side}<0);

        pobj.nii=negvox;
        pobj.name='Negative';
        pobj.niftiFilename='Negative.nii';
        pobj.binary=0;
        pobj.usesolidcolor=0;
        pobj.color=obj.negcolor;
        pobj.colormap=ea_colorgradient(gradientLevel, obj.negcolor, obj.negBaseColor);
        pobj.smooth=10;
        pobj.hullsimplify=0.5;
        pobj.threshold=0;
        obj.drawobject.(subfield){group,side}{2}=ea_roi('Negative.nii',pobj);

        res.img(:)=nansum([res.img(:),-negvox.img(:)],2); % keep copy for export.
    end
    if strcmp(obj.drawTool,'sweetspotmapping')
        if side == 1
            sideCode = 'lh';
        else
            sideCode = 'rh';
        end
        res.img(res.img==0)=nan;
        obj.spotdrawn.(sideCode) = res;
    end

elseif strcmp(obj.vizmode,'Surface (Elvis)')

    % sides=1:2;
    % keep=sides([obj.NMviz.modelRH,obj.NMviz.modelLH]);
    % sides = keep.*sides;
    sides = find([obj.NMviz.modelRH, obj.NMviz.modelLH]);  % gives indices of sides to keep

    % Check cmap
    if exist('voxcmap','var') && ~isempty(voxcmap{group})
        % defaultColor = [1 1 1]; % Default color for nan values
        % cmap = [voxcmap{group}; defaultColor];
        cmap = prep_voxel_colormap(voxcmap{group});
        if isempty(cmap)
            return
        end
    else
        warning('Colormap not defined!')
        return
    end
    res.img(:)=vals{group};
    if isfield(res,'cifti')
        res.cifti.cdata(res.inidx)=vals{group}(res.outidx);
    end
    h=ea_heatmap2surface(res,obj.model,sides,cmap,obj);
    % if ismember(sides,1)
    %     obj.drawobject.networkmapping{group}{1} = h{1};
    % end
    % if ismember(sides,2)
    %     obj.drawobject.networkmapping{group}{2} = h{2};
    % end
    if obj.NMviz.modelRH && obj.NMviz.modelLH
        obj.drawobject.networkmapping{group}{1} = h{1};
        obj.drawobject.networkmapping{group}{2} = h{2};
    elseif obj.NMviz.modelRH
        obj.drawobject.networkmapping{group}{1} = h{1};
    elseif obj.NMviz.modelLH
        obj.drawobject.networkmapping{group}{2} = h{2};
    end

elseif strcmp(obj.vizmode,'Surface (Surfice)')
    res.img(:)=vals{group};
    if isfield(res,'cifti')
        res.cifti.cdata(res.inidx)=vals{group}(res.outidx);
    end
    res.fname=[fileparts(obj.leadgroup),filesep,'model.nii'];

    if ~obj.posvisible
        res.img(res.img>0)=0;
    end

    if ~obj.negvisible
        res.img(res.img<0)=0;
    end
    ea_write_nii(res);
    % det mesh to plot:
    switch obj.model
        case 'Smoothed'
            mesh=([ea_space,'surf_smoothed.mz3']);
            azimuth = '90'; % Right lateral side
            hemiCode = '0'; % Show both hemishperes
        case 'Full'
            mesh=([ea_space,'surf.mz3']);
            azimuth = '90'; % Right lateral side
            hemiCode = '0'; % Show both hemishperes

    end

    threshs=ea_sfc_getautothresh({res.fname});

    script=['BEGIN;',...
        ' RESETDEFAULTS;'...
        ' ORIENTCUBEVISIBLE(FALSE);'];

    script=[script,...
        ' MESHLOAD(''',mesh,''');',...
        ' MESHCOLOR(255,255,255);'];

    cnt=1;

    if ~any(isnan(threshs(1,1:2)))
        script=[script,...
            ' OVERLAYLOAD(''',ea_path_helper(res.fname),''');',...
            ' OVERLAYCOLORNAME(',num2str(cnt),', ''Red-Yellow'');',...
            ' OVERLAYMINMAX(',num2str(cnt),',',num2str(threshs(1,1)),',',num2str(threshs(1,2)),');'];
        cnt=cnt+1;
    end

    if ~any(isnan(threshs(1,3:4)))
        script=[script,...
            ' OVERLAYLOAD(''',ea_path_helper(res.fname),''');',...
            ' OVERLAYCOLORNAME(',num2str(cnt),', ''Blue-Green'');',...
            ' OVERLAYMINMAX(',num2str(cnt),',',num2str(threshs(1,3)),',',num2str(threshs(1,4)),');'];
    end

    script=[script,...
        ' COLORBARVISIBLE(','false',');',...
        ' AZIMUTHELEVATION(',azimuth,', 0);',...
        ' MESHHEMISPHERE(',hemiCode,');'];

    script=[script,...
        ' END.'];

    ea_surfice(script,0);
end
if strcmp(obj.drawTool,'networkmapping')
    res.img(:)=vals{group};
    obj.networkdrawn = res;
end

    function cmap = prep_voxel_colormap(voxcmap_group)
        % Always prepend white as neutral
        if isempty(voxcmap_group)
            warning('Colormap not defined!');
            cmap = [];
            return;
        end
        neutralColor = [1 1 1];
        cmap = [neutralColor; voxcmap_group];
    end

end