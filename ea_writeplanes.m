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
% defaults:
% elstruct=0;
manualtracor=0;
svfig=1;
figvisible='on';

if ~isfield(options,'shifthalfup')
    options.shifthalfup=0;
end

if nargin==1
    % load prior results
    coords_mm=ea_load_reconstruction(options);
    % if there is only one patient to show, ave_coords_mm are the same as the single entry in elstruct(1).coords_mm.
    elstruct(1).coords_mm=coords_mm;
    %fill missing sides with nans, matching it to the other present side.
    %At least one side should be present, which should be always the case.
    elstruct=ea_elstruct_match_and_nanfill(elstruct);
    clear coords_mm
    ave_coords_mm=ea_ave_elstruct(elstruct,options);
        
elseif nargin>1 % elstruct has been supplied, this is a group visualization
    if isstruct(varargin{2})
        elstruct=varargin{2};
        %fill missing sides with nans, matching it to the other present
        %side. At least one side should be present, which should be always
        %the case.
        elstruct=ea_elstruct_match_and_nanfill(elstruct);
        % average coords_mm for image slicing
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

cuts=figure('name',[options.patientname,': 2D cut views...'],'numbertitle','off','Position',[1 scrsz(4)/1.2 scrsz(3)/1.2 scrsz(4)/1.2],'Visible',figvisible,'MenuBar', 'none', 'ToolBar', 'none');
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
        coords{side}=Vtra.mat\[ave_coords_mm{side},ones(size(ave_coords_mm{side},1),1)]';
        coords{side}=coords{side}(1:3,:)';
    end
else
    elstruct=[elstruct,1]';
    coordsmm=elstruct;
    elstruct=manualV.mat\elstruct;
    planedim=ea_getdims(manualtracor,1);
    %elstruct=elstruct(planedim);
end
%XYZ_src_vx = src.mat \ XYZ_mm;

fid=fopen([options.root,options.patientname,filesep,'cuts_export_coordinates.txt'],'w');
for iside=1:length(options.sides)
    side=options.sides(iside);
    %% write out axial/coronal/sagittal images
    for tracor=find(tracorpresent)'
        for elcnt=1:(options.elspec.numel-options.shifthalfup)
            if ~isstruct(elstruct)
                coords={elstruct(1:3)'};
            end
            el=elcnt+options.elspec.numel*(side-1);
            %subplot(2,2,el);

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
                case 2 % coronal images
                    if manualtracor
                        V=manualV;
                    else
                        V=Vcor;
                    end
                case 3 % sagittal images
                    if manualtracor
                        V=manualV;
                    else
                        V=Vsag;
                    end
            end

            [planedim,onedim, secdim , dstring, lstring, Ltxt, Rtxt,plusminusc,plusminusr,plusminusl]=ea_getdims(tracor,side);

            %title(['Electrode ',num2str(el-1),', transversal view.']);

            [slice,~,boundboxmm,sampleheight]=ea_sample_slice(V,dstring,options.d2.bbsize,'mm',coords,el);
            if length(slice)==1 && isnan(slice)
                %there was no electrode here, skip this slice
                continue;
            end
            
            cont=1;
            try                cont=evalin('base','custom_cont'); end
            offs=1;
            try                offs=evalin('base','custom_offs'); end

            slice=ea_contrast(slice,cont,offs);
            level='';

            try                level=evalin('base','level_offset'); end
            if ~isempty(level)

                ms=ea_robustmean(slice(:));
                slice=slice-ms;
                maxmin=[max(slice(:)),min(slice(:))];
                mm=max((maxmin));
                slice(slice>mm)=mm;
                slice(slice<-mm)=-mm;
                %slice(slice<0)=0;
            end

            try                level=evalin('base','level_offset'); end

            printstr_el_stat=['Electrode(s) k',num2str(el-1),'/',options.elspec.contactnames{el} ', ',dstring,' view: ',lstring,'',num2str(sampleheight),' mm.'];
            disp(printstr_el_stat);
            if fid>0 % only if file exists (does sometimes not exist if called from lead anatomy or the slice-cuts feature of elvis)
                fprintf(fid,'%s\n',printstr_el_stat);
            end
            set(0,'CurrentFigure',cuts)

            %             try
            %                 hi=imagesc(slice,...
            %                     [ea_nanmean(slice(slice>0))-3*nanstd(slice(slice>0)) ea_nanmean(slice(slice>0))+3*nanstd(slice(slice>0))]);
            %             catch
            hi=imagesc(slice);
            if ~isempty(level)
                clims=1*[-mm,mm];
                caxis(clims);

            end
            %  end
            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])

            %axis equal
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
                        coordsi{siso}=Viso.mat\[ave_coords_mm{siso},ones(size(ave_coords_mm{siso},1),1)]';
                        coordsi{siso}=coordsi{siso}(1:3,:)';
                    end

                    [slice,~,boundboxmm]=ea_sample_slice(Viso,dstring,options.d2.bbsize,'mm',coordsi,el);

                    try [slicestat]=ea_sample_slice(Visostat,dstring,options.d2.bbsize,'mm',coordsi,el); end
                else
                    coordsi{side}=Viso.mat\[coordsmm(1);coordsmm(1);coordsmm(1);1];
                    coordsi{side}=coordsi{side}(1:3,:)';
                    [slice,~,boundboxmm]=ea_sample_slice(Viso,dstring,options.d2.bbsize,'mm',coordsi,el);
                    try [slicestat]=ea_sample_slice(Visostat,dstring,options.d2.bbsize,'mm',coordsi,el); end
                end
                slice(slice==0)=nan;
                rwholemap=spm_read_vols(Visoraw);
                swholemap=spm_read_vols(Viso);
                wholemap=rwholemap;
                wholemap(~isnan(rwholemap))=swholemap(~isnan(rwholemap));

                %wholemap(~logical(wholemap))=nan;
                maxval=nanmax(wholemap(:));
                minval=nanmin(wholemap(:));

                % define an alpha mask
                alpha=slice;
                alpha(~isnan(alpha))=1;%0.8;
                alpha(isnan(alpha))=0;
                % convert slice to rgb format
                %slicergb=nan([size(slice),3]);

                jetlist=eval(options.prefs.d2.isovolcolormap);
                %jetlist=summer;
                %slice=(slice-minval)/(maxval-minval); % set min max to boundaries 0-1.

                slice(~isnan(slice))=ea_contrast(slice(~isnan(slice)),cont,1);
                slice=slice-1;

                % ##
                % add some "contrast" ? remove this part for linear
                % colormapping
                %
                %                 slice=slice-0.5;
                %                 slice(slice<0)=0;
                %                 slice=slice*2;
                %                 slice(slice>1)=1;

                % ##

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

            %             % Plot L, R and sizelegend
            %             text(addsubsigned(min(boundboxmm{onedim}),2,plusminusl),mean(boundboxmm{secdim}),Ltxt,'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            %             text(addsubsigned(max(boundboxmm{onedim}),2,plusminusr),mean(boundboxmm{secdim}),Rtxt,'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            %
            %             plot([addsubsigned(mean(boundboxmm{onedim}),2.5,'minus'),addsubsigned(mean(boundboxmm{onedim}),2.5,'plus')],[addsubsigned(min(boundboxmm{secdim}),1,'minus'),addsubsigned(min(boundboxmm{secdim}),1,'minus')],'-w','LineWidth',2.5);
            %
            %             text(mean(boundboxmm{onedim}),addsubsigned(min(boundboxmm{secdim}),2,'minus'),'5 mm','color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            %
            %             % Plot slice depth legend
            %             text(mean(boundboxmm{onedim}),addsubsigned(max(boundboxmm{secdim}),2,'minus'),[lstring,sprintf('%.2f',mean(boundboxmm{planedim})),' mm'],'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            %             text(mean(boundboxmm{onedim}),addsubsigned(max(boundboxmm{secdim}),2,plusminusc),[lstring,sprintf('%.2f',mean(boundboxmm{planedim})),' mm'],'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            % Show coordinates
            if isstruct(elstruct)
                if length(elstruct)>1
                    cmap=ea_nice_colors(length(elstruct),options);
                    for pt=1:length(elstruct)
                        ptnames{pt}=elstruct(pt).name;
                    end
                    %ptnames=struct2cell(elstruct.name);
                    %ptnames=squeeze(ptnames(end,1,:))';
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

                if isstruct(elstruct) && ((elstruct(c).activecontacts{side}(elcnt) && options.d3.showactivecontacts) || (~elstruct(c).activecontacts{side}(elcnt) && options.d3.showpassivecontacts))
                    wstr=[1,1,1];
                    if isfield(elstruct,'group')
                    wstr=elstruct(c).groupcolors(elstruct(c).group,:);

                    end
                    estr=wstr./2;
                    if options.d3.hlactivecontacts
                        if elstruct(c).activecontacts{side}(elcnt)
                            %wstr='r';
                            estr='r';
                        end
                    end

                    if options.d2.fid_overlay

                        elplt(c)=plot(elstruct(c).coords_mm{side}(elcnt,onedim),elstruct(c).coords_mm{side}(elcnt,secdim),'o','MarkerSize',10,'MarkerEdgeColor',estr,'MarkerFaceColor',wstr,'LineWidth',2);
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
            drawnow % this is needed here to set alpha below.

            % 3. Dampen alpha by distance (this *has* to be performed
            % last, if not, info is erased by legend again).
            %                 try % not sure if this is supported by earlier ML versions.
            %                     for c=1:length(elstruct)
            %
            %                         dist=abs(diff([elstruct(c).coords_mm{side}(elcnt,planedim),ave_coords_mm{side}(elcnt,planedim)]));
            %                         % dampen alpha by distance
            %                         alp=2*1/exp(dist);
            %                         hMarker = elplt(c).MarkerHandle;
            %                         hMarker.EdgeColorData=uint8(255*[cmap(c,:)';alp]);
            %                     end
            %                 end



            hold off
            %end

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
            %            set(gca, 'LooseInset', [0,0,0,0]);
            %            axis equal
            drawnow
            set(gca,'position',[0,0,1,1],'units','normalized'); % fullscreen plot.
            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])
            drawnow
            expslice=double(frame2im(getframe(cuts))); % export plot.
            expslice=(expslice-min(expslice(:)))/(max(expslice(:))-min(expslice(:))); % set 0 to 1

            expslice=crop_img(expslice);

            if svfig==1 % only export if figure needs to be saved.
                if options.d3.showisovolume
                    isofnadd=['_',options.prefs.d2.isovolsmoothed,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb];
                else
                    isofnadd='';
                end
                switch tracor
                    case 1
                        %saveas(cuts,[options.root,options.patientname,filesep,options.elspec.contactnames{el},'_axial.png']);
                        ea_screenshot([options.root,options.patientname,filesep,options.elspec.contactnames{el},'_axial',isofnadd,'.png'],'myaa');
                    case 2
                        ea_screenshot([options.root,options.patientname,filesep,options.elspec.contactnames{el},'_coronal',isofnadd,'.png'],'myaa');
                    case 3
                        ea_screenshot([options.root,options.patientname,filesep,options.elspec.contactnames{el},'_sagittal',isofnadd,'.png'],'myaa');
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


function y = ea_nanmean(varargin)
if nargin==2
    x=varargin{1};
    dim=varargin{2};
elseif nargin==1
    x=varargin{1};
    dim=1;
end

N = sum(~isnan(x), dim);
y = nansum(x, dim) ./ N;

function img2 = crop_img(img,border)
%
% Crop image by removing edges with homogeneous intensity
%
%
%USAGE
%-----
% img2 = crop_img(img)
% img2 = crop_img(img,border)
%
%
%INPUT
%-----
% - IMG: MxNxC matrix, where MxN is the size of the image and C is the
%   number of color layers
% - BORDER: maximum number of pixels at the borders (default: 0)
%
%
%OUPUT
%-----
% - IMG2: cropped image
%
%
%EXAMPLE
%-------
% >> img  = imread('my_pic.png');
% >> img2 = crop_img(img,0);
% >> imwrite(img2,'my_cropped_pic.png')
%

% Guilherme Coco Beltramini (guicoco@gmail.com)
% 2013-May-29, 12:29 pm


% Input
%==========================================================================
if nargin<2
    border = 0;
end


% Initialize
%==========================================================================
[MM,NN,CC] = size(img);
edge_col   = zeros(2,CC); % image edges (columns)
edge_row   = edge_col;    % image edges (rows)


% Find the edges
%==========================================================================
for cc=1:CC % loop for the colors


    % Top left corner
    %================

    % Find the background
    %--------------------
    img_bg = img(:,:,cc) == img(1,1,cc);

    % Columns
    %--------
    cols = sum(img_bg,1);
    if cols(1)==MM % first column is background
        tmp = find(diff(cols),1,'first'); % background width on the left
        if ~isempty(tmp)
            edge_col(1,cc) = tmp + 1 - border;
        else % no background
            edge_col(1,cc) = 1;
        end
    else % no background
        edge_col(1,cc) = 1;
    end

    % Rows
    %-----
    rows = sum(img_bg,2);
    if rows(1)==NN % first row is background
        tmp = find(diff(rows),1,'first'); % background height at the top
        if ~isempty(tmp)
            edge_row(1,cc) = tmp + 1 - border;
        else % no background
            edge_row(1,cc) = 1;
        end
    else % no background
        edge_row(1,cc) = 1;
    end


    % Bottom right corner
    %====================

    % Find the background
    %--------------------
    img_bg = img(:,:,cc) == img(MM,NN,cc);

    % Columns
    %--------
    cols = sum(img_bg,1);
    if cols(end)==MM % last column is background
        tmp = find(diff(cols),1,'last'); % background width on the right
        if ~isempty(tmp)
            edge_col(2,cc) = tmp + border;
        else % no background
            edge_col(2,cc) = NN;
        end
    else % no background
        edge_col(2,cc) = NN;
    end

    % Rows
    %-----
    rows = sum(img_bg,2);
    if rows(end)==NN % last row is background
        tmp = find(diff(rows),1,'last'); % background height at the bottom
        if ~isempty(tmp)
            edge_row(2,cc) = tmp + border;
        else % no background
            edge_row(2,cc) = MM;
        end
    else % no background
        edge_row(2,cc) = MM;
    end


    % Identify homogeneous color layers
    %==================================
    if edge_col(1,cc)==1 && edge_col(2,cc)==NN && ...
            edge_row(1,cc)==1 && edge_row(2,cc)==MM && ...
            ~any(any(diff(img(:,:,cc),1)))
        edge_col(:,cc) = [NN;1];
        edge_row(:,cc) = [MM;1]; % => ignore layer
    end


end


% Indices of the edges
%==========================================================================

% Columns
%--------
tmp      = min(edge_col(1,:),[],2);
edge_col = [tmp ; max(edge_col(2,:),[],2)];
if edge_col(1)<1
    edge_col(1) = 1;
end
if edge_col(2)>NN
    edge_col(2) = NN;
end

% Rows
%-----
tmp      = min(edge_row(1,:),[],2);
edge_row = [tmp ; max(edge_row(2,:),[],2)];
if edge_row(1)<1
    edge_row(1) = 1;
end
if edge_row(2)>MM
    edge_row(2) = MM;
end


% Crop the edges
%==========================================================================
img2 = zeros(edge_row(2)-edge_row(1)+1,...
    edge_col(2)-edge_col(1)+1,CC,class(img));
for cc=1:CC % loop for the colors
    img2(:,:,cc) = img(edge_row(1):edge_row(2),edge_col(1):edge_col(2),cc);
end



function coords_mm=ea_ave_elstruct(elstruct,options)
% simply averages coordinates of a group to one coords_mm 1x2 cell
coords_mm=elstruct(1).coords_mm; % initialize mean variable
for side=1:length(coords_mm)

    for xx=1:size(coords_mm{side},1)
        for yy=1:size(coords_mm{side},2)
            vals=zeros(length(elstruct),1);
            for vv=1:length(elstruct)
                if ~isempty(elstruct(vv).coords_mm{side})
                    vals(vv)=elstruct(vv).coords_mm{side}(xx,yy);
                end
            end
            coords_mm{side}(xx,yy)=ea_robustmean(vals);

        end
    end
end
if options.shifthalfup
    for side=1:length(coords_mm)
        for c=1:length(coords_mm{side})-1
            scoords_mm{side}(c,:)=mean([coords_mm{side}(c,:);coords_mm{side}(c+1,:)],1);

        end
    end
    coords_mm=scoords_mm;
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

function str=sub2space(str) % replaces subscores with spaces
str(str=='_')=' ';








