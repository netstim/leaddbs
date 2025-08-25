function ea_unified_draw_volumetric(obj,vals,group,side,gradientLevel,voxcmap)
if strcmp(obj.drawTool,'sweetspotmapping')
    % Plot voxels if any survived
    res=obj.results.sweetspotmapping.space{side};
    res.dt(1) = 16;
    res.img(:)=nan;
    obj.vizmode = 'sweetspot';
elseif strcmp(obj.drawTool,'networkmapping')
    if contains(obj.calcsettings.netmap_connectome, ' > ')
        % For functional connectome, should use the spacedef provided by the connectome itself
        res = obj.results.networkmapping.(ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome)).space;
    else
        res = ea_load_nii([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,obj.outputspace,'.nii.gz']);
    end
    res.dt(1) = 16;


end
%plot based on vizmode
if strcmp(obj.vizmode,'sweetspot') || strcmp(obj.vizmode,'Regions')
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
        obj.drawobject.sweetspotmapping{group,side}{1}=ea_roi('Positive.nii',pobj);

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
        obj.drawobject.sweetspotmapping{group,side}{2}=ea_roi('Negative.nii',pobj);

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

    sides=1:2;
    % Check cmap
    if exist('voxcmap','var') && ~isempty(voxcmap{group})
        defaultColor = [1 1 1]; % Default color for nan values
        cmap = [voxcmap{group}; defaultColor];
    else
        warning('Colormap not defined!')
        return
    end
    res.img(:)=vals{group};
    if isfield(res,'cifti')
        res.cifti.cdata(res.inidx)=vals{group}(res.outidx);
    end
    h=ea_heatmap2surface(res,obj.model,sides,cmap,obj);
    obj.drawobject.networkmapping{group}{1} = h{1};
    obj.drawobject.networkmapping{group}{2} = h{2};


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
end