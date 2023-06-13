function [cuts,expslice,boundboxmm,allcontour]=ea_writeplanes(varargin)

% This function exports slice views of all electrode contacts reconstructed
% priorly. Images are written as .png image files. Bot transversal and
% coronal views are being exported. Additionally, overlays from atlas-data
% can be visualized via the function ea_add_overlay which uses all atlas
% files that are found in the lead_dbs atlas directory.
% inputs: options (struct using standard lead-dbs fields), optional:
% elstruct (for group visualization).

% usage: [cuts,expslice]=ea_writeplanes(options, elstruct or manualheight (in mm), manualtracor,manualV, visible ['on'/'off'], savefig [1/0])

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

options=varargin{1};
manualtracor=0;
svfig=1;
figvisible='on';

if ~isfield(options,'shifthalfup')
    options.shifthalfup=0;
end

if nargin==1
    % load prior results
    coords_mm=ea_load_reconstruction(options);
    % If there is only one patient to show, ave_coords_mm are the same as the single entry in elstruct(1).coords_mm.
    elstruct(1).coords_mm = coords_mm;
    ave_coords_mm = coords_mm;
elseif nargin>1 % elstruct has been supplied, this is a group visualization
    if isstruct(varargin{2})
        elstruct=varargin{2};
        ave_coords_mm=ea_ave_elstruct(elstruct,options);
    else % concrete height is being supplied (without electrode star plotting).
        elstruct=varargin{2};
        options.elspec.numel=1; % only iterate once below.
    end
end

if nargin>5 % also has flags to hide and not to save the result (as it will be used by 3D-figure).
    manualtracor=varargin{3}; % manually specify if to export tra, cor or sag image.
    manualV=varargin{4};
    figvisible=varargin{5};
    svfig=varargin{6};
    try
        atlases=varargin{7};
    end
end

if svfig
    disp('Exporting 2D slice output...');
end

scrsz = get(0,'ScreenSize');

titilePrefix = erase(options.patientname, 'sub-');

cuts=figure('name',[titilePrefix,': 2D cut views...'],'numbertitle','off','Position',[1 scrsz(4)/1.2 scrsz(3)/1.2 scrsz(4)/1.2],'Visible',figvisible,'MenuBar', 'none', 'ToolBar', 'none');
axis off
set(cuts,'position',[100, 100, 800 ,800]);
set(cuts,'color','w');
tracorpresent=zeros(3,1); % check if files are present.

if ~manualtracor
    tracorpresent=ones(3,1);
    [Vtra,Vcor,Vsag]=ea_assignbackdrop(options.d2.backdrop,options,'Patient');
else
    tracorpresent(manualtracor)=1; % only export specified orientation.
end

if isstruct(elstruct)
    for side=1:length(ave_coords_mm)
        coords{side}=ea_mm2vox(ave_coords_mm{side}, Vtra.mat);
    end
else
    coordsmm = elstruct;
    elstruct = ea_mm2vox(elstruct,manualV.mat);
    coords = elstruct;
    planedim=ea_getdims(manualtracor,1);
end

if startsWith(options.patientname, 'gs_')
    export_2D_folder = fullfile(options.root, options.patientname, '2D');
else
    export_2D_folder = fullfile(options.root, options.patientname, 'export', '2D');
end
ea_mkdir(export_2D_folder);
fileBaseName = fullfile(export_2D_folder, [options.patientname, '_desc-2D_']);

fid = fopen([fileBaseName 'viewplane.txt'],'w');
for iside=1:length(options.sides)
    side=options.sides(iside);
    % write out axial/coronal/sagittal images
    for tracor=find(tracorpresent)'
        for elcnt=1:(options.elspec.numel-options.shifthalfup)
            el=elcnt+options.elspec.numel*(side-1);

            % Show MR-volume
            set(0,'CurrentFigure',cuts)
            colormap(gray)
            try
                custom_cmap=evalin('base','custom_cmap');
                colormap(custom_cmap);
            end

            switch tracor
                case 1 % transversal images
                    if manualtracor
                        V=manualV;
                    else
                        V=Vtra;
                    end
                    if isstruct(elstruct) && length(options.sides)==2
                        showBothSide = 1;
                        el = [el, elcnt+options.elspec.numel];
                        if side==2
                            continue;
                        end
                    end
                case 2 % coronal images
                    if manualtracor
                        V=manualV;
                    else
                        V=Vcor;
                    end
                    if isstruct(elstruct) && length(options.sides)==2
                        showBothSide = 1;
                        el = [el, elcnt+options.elspec.numel];
                        if side==2
                            continue;
                        end
                    end
                case 3 % sagittal images
                    if manualtracor
                        V=manualV;
                    else
                        V=Vsag;
                    end
                    showBothSide = 0;
            end

            [planedim,onedim, secdim , dstring, lstring, Ltxt, Rtxt,plusminusc,plusminusr,plusminusl]=ea_getdims(tracor,side);

            [slice,~,boundboxmm,sampleheight]=ea_sample_slice(V,dstring,options.d2.bbsize,'mm',coords,el);
            if length(slice)==1 && isnan(slice)
                %there was no electrode here, skip this slice
                continue;
            end

            cont=1;
            try
                cont=evalin('base','custom_cont');
            end

            offs=1;
            try
                offs=evalin('base','custom_offs');
            end

            slice=ea_contrast(slice,cont,offs);
            level='';

            try
                level=evalin('base','level_offset');
            end

            if ~isempty(level)
                ms=ea_robustmean(slice(:));
                slice=slice-ms;
                maxmin=[max(slice(:)),min(slice(:))];
                mm=max((maxmin));
                slice(slice>mm)=mm;
                slice(slice<-mm)=-mm;
            end

            try
                level=evalin('base','level_offset');
            end

            try
                printstr_el_stat=['Electrode(s) ', strjoin(options.elspec.contactnames(el)) ', ',dstring,' view: ',lstring,'',num2str(sampleheight),' mm.'];
                disp(printstr_el_stat);
                
                if fid>0 % only if file exists (does sometimes not exist if called from lead anatomy or the slice-cuts feature of elvis)
                    fprintf(fid,'%s\n',printstr_el_stat);
                end
            end

            set(0,'CurrentFigure',cuts)

            hi=imagesc(slice);
            if ~isempty(level)
                clims=1*[-mm,mm];
                caxis(clims);
            end

            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])

            hold on

            if manualtracor
                aspectratio=V.dim(onedim)/V.dim(secdim);
            else
                aspectratio=1;
            end

            % Show overlays
            if options.d2.writeatlases
                if exist('atlases', 'var')
                    options.atlases = atlases;
                end
                if isfield(options, 'atlases') || ~strcmp(options.atlasset, 'Use none')
                    [cuts,allcontour] = ea_add_overlay(boundboxmm, cuts, tracor, options);
                end
            end

            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])

            % Show isovolume
            if options.d3.showisovolume
                Visoraw=spm_vol([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb,'.nii']);
                Viso=spm_vol([options.root,options.patientname,filesep,options.prefs.d2.isovolsmoothed,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb,'.nii']);
                try
                    Visostat=spm_vol([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb,'_p05.nii']);
                end

                if isstruct(elstruct)
                    for siso=1:length(ave_coords_mm)
                        coordsi{siso}=ea_mm2vox(ave_coords_mm{siso}, Viso.mat);
                    end

                    [slice,~,boundboxmm]=ea_sample_slice(Viso,dstring,options.d2.bbsize,'mm',coordsi,el);
                    try
                        [slicestat]=ea_sample_slice(Visostat,dstring,options.d2.bbsize,'mm',coordsi,el);
                    end
                else
                    coordsi{side}=ea_mm2vox([coordsmm(1),coordsmm(1),coordsmm(1)], Viso.mat);
                    [slice,~,boundboxmm]=ea_sample_slice(Viso,dstring,options.d2.bbsize,'mm',coordsi,el);
                    try
                        [slicestat]=ea_sample_slice(Visostat,dstring,options.d2.bbsize,'mm',coordsi,el);
                    end
                end

                slice(slice==0)=nan;
                rwholemap=spm_read_vols(Visoraw);
                swholemap=spm_read_vols(Viso);
                wholemap=rwholemap;
                wholemap(~isnan(rwholemap))=swholemap(~isnan(rwholemap));
                maxval=nanmax(wholemap(:));
                minval=nanmin(wholemap(:));

                % define an alpha mask
                alpha=slice;
                alpha(~isnan(alpha))=1;
                alpha(isnan(alpha))=0;

                jetlist=eval(options.prefs.d2.isovolcolormap);

                slice(~isnan(slice))=ea_contrast(slice(~isnan(slice)),cont,1);
                slice=slice-1;

                slice=round(slice.*(length(gray)-1))+1; % set min max to boundaries 1-length(gray).
                slice(slice<1)=1; slice(slice>length(gray))=length(gray);

                slicer=slice; sliceg=slice; sliceb=slice;
                slicer(~isnan(slicer))=jetlist(slicer(~isnan(slicer)),1);
                sliceg(~isnan(sliceg))=jetlist(sliceg(~isnan(sliceg)),2);
                sliceb(~isnan(sliceb))=jetlist(sliceb(~isnan(sliceb)),3);
                slicergb=cat(3,slicer,sliceg,sliceb);
                isv=imagesc(slicergb);
                set(isv,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
                set(isv,'AlphaData',alpha);

                % draw significance countour if available:
                try
                    slicestat(isnan(slicestat))=0;
                    warning('off')
                    [cmat,statcontour]=contour(slicestat,1);
                    set(statcontour,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
                    set(statcontour,'Color','w');
                    warning('on')
                end
            end

            % Plot orientation label
            % text(addsubsigned(min(boundboxmm{onedim}),2,plusminusl),mean(boundboxmm{secdim}),Ltxt,'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            % text(addsubsigned(max(boundboxmm{onedim}),2,plusminusr),mean(boundboxmm{secdim}),Rtxt,'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');

            % Plot scale
            % plot([addsubsigned(mean(boundboxmm{onedim}),2.5,'minus'),addsubsigned(mean(boundboxmm{onedim}),2.5,'plus')],[addsubsigned(min(boundboxmm{secdim}),1,'minus'),addsubsigned(min(boundboxmm{secdim}),1,'minus')],'-w','LineWidth',2.5);
            % text(mean(boundboxmm{onedim}),addsubsigned(min(boundboxmm{secdim}),2.5,'minus'),'5 mm','color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');

            % Plot slice depth
            % text(mean(boundboxmm{onedim}),addsubsigned(max(boundboxmm{secdim}),2,'minus'),[lstring,sprintf('%.2f',mean(boundboxmm{planedim})),' mm'],'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');

            % Show coordinates
            if isstruct(elstruct)
                if length(elstruct)>1
                    cmap=ea_nice_colors(length(elstruct),options);
                    for pt=1:length(elstruct)
                        ptnames{pt}=elstruct(pt).name;
                    end
                else
                    cmap=[0.9,0.9,0.9];
                end
            end

            % 1. Plot stars
            for c=1:length(elstruct)
                % prepare active/passive contacts
                if ~isfield(elstruct(c),'elmodel') % usually, elspec is defined by the GUI. In case of group analyses, for each patient, a different electrode model can be selected for rendering.
                    elspec=options.elspec;
                else % if elspec is defined for each electrode, overwrite options-struct settings here.
                    o=ea_resolve_elspec(elstruct(c));
                    elspec=o.elspec; clear o
                end

                elstruct=testifactivecontacts(elstruct,elspec,c); % small function that tests if active contacts are assigned and if not assigns them all as passive.

                if isstruct(elstruct)
                    if elstruct(c).activecontacts{side}(elcnt) && options.d3.showactivecontacts ...
                            || ~elstruct(c).activecontacts{side}(elcnt) && options.d3.showpassivecontacts
                        faceColor=[1,1,1];
                        if isfield(elstruct,'group')
                            faceColor=elstruct(c).groupcolors(elstruct(c).group,:);
                        end

                        edgeColor=faceColor./2;
                        if options.d3.hlactivecontacts
                            if elstruct(c).activecontacts{side}(elcnt)
                                edgeColor='r';
                            end
                        end

                        if options.d2.fid_overlay
                            elplt(c)=plot(elstruct(c).coords_mm{side}(elcnt,onedim),elstruct(c).coords_mm{side}(elcnt,secdim),'o','MarkerSize',10,'MarkerEdgeColor',edgeColor,'MarkerFaceColor',faceColor,'LineWidth',2);
                        end
                    end

                    if showBothSide ...
                            && (elstruct(c).activecontacts{side+1}(elcnt) && options.d3.showactivecontacts ...
                            || ~elstruct(c).activecontacts{side+1}(elcnt) && options.d3.showpassivecontacts)
                        faceColor=[1,1,1];
                        if isfield(elstruct,'group')
                            faceColor=elstruct(c).groupcolors(elstruct(c).group,:);
                        end

                        edgeColor=faceColor./2;
                        if options.d3.hlactivecontacts
                            if elstruct(c).activecontacts{side+1}(elcnt)
                                edgeColor='r';
                            end
                        end

                        if options.d2.fid_overlay
                            elplt(c)=plot(elstruct(c).coords_mm{side+1}(elcnt,onedim),elstruct(c).coords_mm{side+1}(elcnt,secdim),'o','MarkerSize',10,'MarkerEdgeColor',edgeColor,'MarkerFaceColor',faceColor,'LineWidth',2);
                        end
                    end
                end
            end

            % 2. Plot legend
            if exist('elplt','var') % if no stars have been plottet, no legend is needed.
                if numel(elplt)>1
                    if options.d2.showlegend
                        if exist('ptnames','var')
                            if numel(elplt)>50
                                cols=round(sqrt(numel(elplt(:))));
                                if cols>6; cols=6; end
                                ea_columnlegend(cols,elplt,ptnames,'Location','Middle');
                            else
                                legend(elplt,ptnames,'Location','southoutside','Orientation','Horizontal','FontSize',9,'FontWeight','bold');
                                legend('boxoff');
                            end
                        end
                    end
                end
            end
            axis xy
            axis off
            drawnow

            hold off

            axis equal
            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])

            % Save results
            if strcmp(figvisible,'on')
                set(cuts,'visible','on');
            end

            axis off

            if svfig>0
                set(cuts,'position',[100, 100, 600*aspectratio,600]);
            elseif svfig==0
                set(cuts,'position',[100, 100, 3200*aspectratio,3200]);
            end

            set(0,'CurrentFigure',cuts)
            drawnow

            set(gca,'position',[0,0,1,1],'units','normalized'); % fullscreen plot.
            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])
            drawnow

            expslice=double(frame2im(getframe(cuts))); % export plot.
            expslice=(expslice-min(expslice(:)))/(max(expslice(:))-min(expslice(:))); % set 0 to 1
            expslice=ea_crop_img(expslice);

            if svfig==1 % only export if figure needs to be saved.
                if options.d3.showisovolume % TODO: Fix export name
                    isofnadd=['_',options.prefs.d2.isovolsmoothed,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb];
                else
                    isofnadd='';
                end
                contactTag = strjoin(options.elspec.contactnames(el), '');
                contactTag = erase(contactTag, {' ', '(', ')'});
                desc = [contactTag, isofnadd];
                switch tracor
                    case 1
                        ea_screenshot([fileBaseName 'view-ax_',desc,'.png'],'myaa');
                    case 2
                        ea_screenshot([fileBaseName 'view-cor_',desc,'.png'],'myaa');
                    case 3
                        ea_screenshot([fileBaseName 'view-sag_',desc,'.png'],'myaa');
                end
            end
            axis xy
        end
    end
end

if fid>0
    fclose(fid);
end

if svfig<2
    close(cuts)
end

if svfig
    disp('Done.');
end


function meanCoords = ea_ave_elstruct(elstruct, options)
% Calculate the mean coordinates of a group of patients

coords = {elstruct.coords_mm}';

% Get max number of electrodes for the group of patients
numSides = max(cellfun(@length, coords));

% Get max number of contacts for the group of patients
numContacts = zeros(1, numSides);
for p = 1:length(coords)
    for s = 1:numSides
        if length(coords{p}) >= s && size(coords{p}{s},1) > numContacts(s)
            numContacts(s) = size(coords{p}{s},1);
        end
    end
end

% Reformat coordinates matrix to the same [maximum] size
for p = 1:length(coords)
    for s = 1:numSides
        temp = nan(numContacts(s),3);
        if length(coords{p}) < s
            coords{p}{s} = temp;
        elseif size(coords{p}{s},1) < numContacts(s)
            temp(1:size(coords{p}{s},1), :) = coords{p}{s};
            coords{p}{s} = temp;
        end
    end
end

% Reshap coordinates to P x S cell array
coords = vertcat(coords{:});

% Calculate the mean coordinates
meanCoords = cell(1, numSides);
for s = 1:numSides
    meanCoords{s} = nan(numContacts(s),3);
    temp = cat(3, coords{:,s});
    for c = 1:numContacts(s)
        for d = 1:3
           meanCoords{s}(c,d) = ea_robustmean(squeeze(temp(c,d,:)));
        end
    end
end

% Shift-up option
if options.shifthalfup
    shiftedMeanCoords = cell(size(meanCoords));
    for side = 1:length(meanCoords)
        for c = 1:size(meanCoords{side},1)-1
            shiftedMeanCoords{side}(c,:) = mean([meanCoords{side}(c,:); meanCoords{side}(c+1,:)], 1);
        end
    end
    meanCoords = shiftedMeanCoords;
end


function val=addsubsigned(val,add,command)
switch command
    case 'plus'
        if val>0
            val=val+add;
        elseif val<0
            val=val-add;
        end
    case 'minus'
        if val>0
            val=val-add;
        elseif val<0
            val=val+add;
        end
end


function elstruct=testifactivecontacts(elstruct,elspec,c)
if ~isstruct(elstruct)
    return
end

if ~isfield(elstruct(c),'activecontacts')
    elstruct(c).activecontacts{1}=zeros(elspec.numel,1);
    elstruct(c).activecontacts{2}=zeros(elspec.numel,1);
else
    if isempty(elstruct(c).activecontacts)
        elstruct(c).activecontacts{1}=zeros(elspec.numel,1);
        elstruct(c).activecontacts{2}=zeros(elspec.numel,1);
    end
end
