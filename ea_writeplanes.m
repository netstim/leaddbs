function [cuts,expslice]=ea_writeplanes(varargin)

% This function exports slice views of all electrode contacts reconstructed
% priorly. Images are written as .png image files. Bot transversal and
% coronar views are being exported. Additionally, overlays from atlas-data
% can be visualized via the function ea_add_overlay which uses all atlas
% files that are found in the lead_dbs/atlases directory.
% inputs: options (struct using standard lead-dbs fields), optional:
% elstruct (for group visualization).

% usage: [cuts,expslice]=ea_writeplanes(options, elstruct or manualheight (in mm), manualtracor,manualV, visible ['on'/'off'], savefig [1/0])

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn



options=varargin{1};
% defaults:
%elstruct=0;
manualtracor=0;
svfig=1;
figvisible='on';
if nargin==1
    % load prior results
    coords_mm=ea_load_reconstruction(options);
    ave_coords_mm=coords_mm;
    clear coords_mm
    elstruct(1).coords_mm=ave_coords_mm; % if there is only one patient to show, ave_coords_mm are the same as the single entry in elstruct(1).coords_mm.
    
elseif nargin>1 % elstruct has been supplied, this is a group visualization
    if isstruct(varargin{2})
        elstruct=varargin{2};
        % average coords_mm for image slicing
        ave_coords_mm=ea_ave_elstruct(elstruct);
    else % concrete height is being supplied (without electrode star plotting).
        
        elstruct=varargin{2};
        options.elspec.numel=1; % only iterate once below.
    end
end

if nargin==7 % also has flags to hide and not to save the result (as it will be used by 3D-figure).
    manualtracor=varargin{3}; % manually specify if to export tra, cor or sag image.
    manualV=varargin{4};
    figvisible=varargin{5};
    svfig=varargin{6};
    atlases=varargin{7};
end
if svfig
    disp('Exporting 2D slice output...');
end

if strcmp(options.prefs.d2.useprepost,'pre') % use preoperative images, overwrite filenames to preoperative version
    options.prefs.gtranii=options.prefs.gprenii;
    options.prefs.tranii=options.prefs.prenii;
    options.prefs.gcornii=options.prefs.gprenii;
    options.prefs.cornii=options.prefs.prenii;
    options.prefs.gsagnii=options.prefs.gprenii;
    options.prefs.sagnii=options.prefs.prenii;
elseif strcmp(options.prefs.d2.useprepost,'template')
    
    options.modality=3;
end

% resolve 2d-options from td_options.mat
options=resolve2doptions(options);

scrsz = get(0,'ScreenSize');


cuts=figure('name',[options.patientname,': 2D cut views (figure is being saved)...'],'numbertitle','off','Position',[1 scrsz(4)/1.2 scrsz(3)/1.2 scrsz(4)/1.2],'Visible',figvisible);
axis off
set(cuts,'position',[100, 100, 800 ,800]);
set(cuts,'color','w');
tracorpresent=zeros(3,1); % check if files are present.
if ~manualtracor
    
    switch options.modality
        case 1 % MR
            try
                Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
                tracorpresent(1)=1;
            catch
                try
                    Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
                    tracorpresent(1)=1;
                end
                
            end
            try
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii));
                tracorpresent(2)=1;
                
            catch
                try
                    Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii));
                    tracorpresent(1)=1;
                end
            end
            try
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gsagnii));
                tracorpresent(3)=1;
                
            catch
                try
                    Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.sagnii));
                    tracorpresent(3)=1;
                end
            end
        case 2 % CT
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
            Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
            Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
            tracorpresent(1:3)=1;
        case 3 % use template
            Vtra=spm_vol(fullfile(options.earoot,'templates','mni_hires_bb.nii'));
            Vcor=spm_vol(fullfile(options.earoot,'templates','mni_hires_bb.nii'));
            Vsag=spm_vol(fullfile(options.earoot,'templates','mni_hires_bb.nii'));
            tracorpresent(1:3)=1;
    end
else
    tracorpresent(manualtracor)=1; % only export specified orientation.
end


if isstruct(elstruct)
    for side=1:length(ave_coords_mm)
        coords{side}=Vtra.mat\[ave_coords_mm{side},ones(size(ave_coords_mm{side},1),1)]';
        coords{side}=coords{side}(1:3,:)';
    end
    
else
    elstruct=[repmat(elstruct,3,1);1];
    coordsmm=elstruct;
    elstruct=manualV.mat\elstruct;
    planedim=getdims(manualtracor,1);
    elstruct=elstruct(planedim);
end
%XYZ_src_vx = src.mat \ XYZ_mm;


for side=1:length(options.sides)
    %% write out axial images
    for tracor=find(tracorpresent)'
        
        for elcnt=1:options.elspec.numel
            
            
            if ~isstruct(elstruct)
                coords={elstruct};
            end
            el=elcnt+options.elspec.numel*(side-1);
            %subplot(2,2,el);
            
            % Show MR-volume
            set(0,'CurrentFigure',cuts)
            colormap gray
            
            
            switch tracor
                
                case 1 % transversal images
                    if manualtracor
                        V=manualV;
                    else
                        V=Vtra;
                    end
                case 2 % coronar images
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
            
            
            [planedim,onedim, secdim , dstring, lstring, Ltxt, Rtxt,plusminusc,plusminusr,plusminusl]=getdims(tracor,side);
            
            
            
            %title(['Electrode ',num2str(el-1),', transversal view.']);
            
            [slice,~,boundboxmm,sampleheight]=ea_sample_slice(V,dstring,options.d2.bbsize,'mm',coords,el);
            disp(['Electrode(s) k',num2str(el-1),', ',dstring,' view: ',lstring,'',num2str(sampleheight),' mm.']);
            
            set(0,'CurrentFigure',cuts)
            try
                hi=imagesc(slice,...
                    [ea_nanmean(slice(slice>0))-3*nanstd(slice(slice>0)) ea_nanmean(slice(slice>0))+3*nanstd(slice(slice>0))]);
            catch
                hi=imagesc(slice);
                
            end
            set(hi,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
            axis([min(boundboxmm{onedim}),max(boundboxmm{onedim}),min(boundboxmm{secdim}),max(boundboxmm{secdim})])
            axis square
            hold on
            
            
            
            
            % Show overlays
            if options.d2.writeatlases
                try options.atlases=atlases; end
                cuts=ea_add_overlay(boundboxmm,cuts,tracor,options);
            end
            
            
            
            % Show isovolume
            
            if options.d3.showisovolume
                Visoraw=spm_vol([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb,'.nii']);
                Viso=spm_vol([options.root,options.patientname,filesep,options.prefs.d2.isovolsmoothed,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb,'.nii']);
                Visostat=spm_vol([options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb,'_p05.nii']);
                if exist('ave_coords_mm','var')
                    for siso=1:length(ave_coords_mm)
                        coordsi{siso}=Viso.mat\[ave_coords_mm{siso},ones(size(ave_coords_mm{siso},1),1)]';
                        coordsi{siso}=coordsi{siso}(1:3,:)';
                    end
                    [slice,~,boundboxmm]=ea_sample_slice(Viso,dstring,options.d2.bbsize,'mm',coordsi,el);
                    [slicestat]=ea_sample_slice(Visostat,dstring,options.d2.bbsize,'mm',coordsi,el);
                    
                else
                    coordsi{side}=Viso.mat\[coordsmm(1);coordsmm(1);coordsmm(1);1];
                    coordsi{side}=coordsi{side}(1:3,:)';
                    [slice,~,boundboxmm]=ea_sample_slice(Viso,dstring,options.d2.bbsize,'mm',coordsi,el);
                    [slicestat]=ea_sample_slice(Visostat,dstring,options.d2.bbsize,'mm',coordsi,el);
                    
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
                alpha(~isnan(alpha))=0.8;
                alpha(isnan(alpha))=0;
                % convert slice to rgb format
                %slicergb=nan([size(slice),3]);
                
                jetlist=eval(options.prefs.d2.isovolcolormap);
                %jetlist=summer;
                slice=(slice-minval)/(maxval-minval); % set min max to boundaries 0-1.
                
                % ##
                % add some "contrast" ? remove this part for linear
                % colormapping
                
                slice=slice-0.5;
                slice(slice<0)=0;
                slice=slice*2;
                slice(slice>1)=1;
                
                % ##
                
                slice=round(slice.*63)+1; % set min max to boundaries 1-64.
                slice(slice<1)=1; slice(slice>64)=64;
                
                slicer=slice; sliceg=slice; sliceb=slice;
                slicer(~isnan(slicer))=jetlist(slicer(~isnan(slicer)),1);
                sliceg(~isnan(sliceg))=jetlist(sliceg(~isnan(sliceg)),2);
                sliceb(~isnan(sliceb))=jetlist(sliceb(~isnan(sliceb)),3);
                slicergb=cat(3,slicer,sliceg,sliceb);
                isv=imagesc(slicergb);
                set(isv,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
                set(isv,'AlphaData',alpha);
                
                % draw significance countour:
                slicestat(isnan(slicestat))=0;
                warning('off')
                [cmat,statcontour]=contour(slicestat,1);
                set(statcontour,'XData',boundboxmm{onedim},'YData',boundboxmm{secdim});
                set(statcontour,'Color','w');
                warning('on')
                
            end
            
            
            
            
            % Plot L, R and sizelegend
            %text(addsubsigned(min(boundboxmm{onedim}),2,plusminusl),mean(boundboxmm{secdim}),Ltxt,'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            %text(addsubsigned(max(boundboxmm{onedim}),2,plusminusr),mean(boundboxmm{secdim}),Rtxt,'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            
            %plot([addsubsigned(mean(boundboxmm{onedim}),2.5,'minus'),addsubsigned(mean(boundboxmm{onedim}),2.5,'plus')],[addsubsigned(min(boundboxmm{secdim}),1,'minus'),addsubsigned(min(boundboxmm{secdim}),1,'minus')],'-w','LineWidth',2.5);
            
            %text(mean(boundboxmm{onedim}),addsubsigned(min(boundboxmm{secdim}),2,'minus'),'5 mm','color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            
            % Plot slice depth legend
            %text(mean(boundboxmm{onedim}),addsubsigned(max(boundboxmm{secdim}),2,'minus'),[lstring,sprintf('%.2f',mean(boundboxmm{planedim})),' mm'],'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
            %text(mean(boundboxmm{onedim}),addsubsigned(max(boundboxmm{secdim}),2,plusminusc),[lstring,sprintf('%.2f',mean(boundboxmm{planedim})),' mm'],'color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',40,'FontWeight','bold');
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
                    
                    if (elstruct(c).activecontacts{side}(elcnt) && options.d3.showactivecontacts) || (~elstruct(c).activecontacts{side}(elcnt) && options.d3.showpassivecontacts)
                        
                        wstr='w';
                        if options.d3.hlactivecontacts
                            
                            
                            if elstruct(c).activecontacts{side}(elcnt)
                                wstr='r';
                                
                            end
                        end
                        elplt(c)=plot(elstruct(c).coords_mm{side}(elcnt,onedim),elstruct(c).coords_mm{side}(elcnt,secdim),'*','MarkerSize',15,'MarkerEdgeColor',wstr,'MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',4,'LineSmoothing','on');
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
            end
            
            %set(gca,'LooseInset',get(gca,'TightInset'))
            % Save results
            if strcmp(figvisible,'on')
                set(cuts,'visible','on');
            end
            axis off
            if svfig
                set(cuts,'position',[100, 100, 600 ,600]);
            else
                set(cuts,'position',[100, 100, 3200 ,3200]);
            end
            set(0,'CurrentFigure',cuts)
            set(gca,'position',[0,0,1,1],'units','normalized'); % fullscreen plot.
            expslice=double(frame2im(getframe(cuts))); % export plot.
            if svfig % only export if figure needs to be saved.
                if options.d3.showisovolume
                    isofnadd=['_',options.prefs.d2.isovolsmoothed,options.d3.isomatrix_name,'_',options.prefs.d2.isovolsepcomb];
                else
                    isofnadd='';
                end
                switch tracor
                    case 1
                        %saveas(cuts,[options.root,options.patientname,filesep,options.elspec.contactnames{el},'_axial.png']);
                        ea_screenshot([options.root,options.patientname,filesep,options.elspec.contactnames{el},'_axial',isofnadd,'.png']);
                    case 2
                        ea_screenshot([options.root,options.patientname,filesep,options.elspec.contactnames{el},'_coronar',isofnadd,'.png']);
                    case 3
                        ea_screenshot([options.root,options.patientname,filesep,options.elspec.contactnames{el},'_sagittal',isofnadd,'.png']);
                end
            end
            
            axis xy
            
        end
    end
end




close(cuts)
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





function coords_mm=ea_ave_elstruct(elstruct)
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
            coords_mm{side}(xx,yy)=mean(vals);
            
        end
    end
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



function [planedim,onedim, secdim, dstring, lstring, Ltxt, Rtxt,plusminusc,plusminusr,plusminusl]=getdims(tracor,side)

switch tracor
    
    case 1 % transversal images
        onedim=1;
        secdim=2;
        planedim=3;
        dstring='tra';
        lstring='z = ';
        Ltxt='M';
        Rtxt='L';
        plusminusc='plus';
        switch side
            case 1
                plusminusr='minus';
                plusminusl='plus';
            case 2
                plusminusr='plus';
                plusminusl='minus';
        end
    case 2 % coronar images
        onedim=1;
        secdim=3;
        planedim=2;
        dstring='cor';
        lstring='y = ';
        Ltxt='M';
        Rtxt='L';
        plusminusc='minus';
        
        switch side
            case 1
                plusminusr='minus';
                plusminusl='plus';
            case 2
                plusminusr='plus';
                plusminusl='minus';
        end
    case 3 % saggital images
        onedim=2;
        secdim=3;
        planedim=1;
        dstring='sag';
        lstring='x = ';
        Ltxt='P';
        Rtxt='A';
        plusminusc='minus';
        switch side
            case 1
                plusminusr='plus';
                plusminusl='minus';
            case 2
                plusminusr='plus';
                plusminusl='minus';
        end
        
end


function options=resolve2doptions(options)

try
    tdhandles=load([options.earoot,'td_options.mat']);
    set(tdhandles.ea_spec2dwrite,'visible','off');
    options.d2.col_overlay=get(tdhandles.tdcolorscheck,'Value');
    options.d2.con_overlay=get(tdhandles.tdcontourcheck,'Value');
    options.d2.con_color=getappdata(tdhandles.tdcontourcolor,'color');
    if isempty(options.d2.con_color)
        options.d2.con_color=[1,1,1]; % white
    end
    options.d2.lab_overlay=get(tdhandles.tdlabelcheck,'Value');
    options.d2.bbsize=str2double(get(tdhandles.bbsize,'String'));
    delete(tdhandles.ea_spec2dwrite);
catch % defaults
    options.d2.col_overlay=1;
    options.d2.con_overlay=1;
    options.d2.con_color=[1,1,1];
    options.d2.lab_overlay=1;
    options.d2.bbsize=50;
end

