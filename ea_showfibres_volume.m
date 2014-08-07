function [PL]=ea_showfibres_volume(resultfig,options)
% This function shows fiber-connectivity from a volume defined by a nx3
% point-list (volume). If stimparams.showconnectivities is set, connected
% areas to the volume are also visualized. To do so, the function uses
% inhull.m which is covered by the BSD-license (see below).
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

PL.ht=uitoolbar(resultfig);
set(0,'CurrentFigure',resultfig)
colornames='rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk';


hold on
% get app data

stimparams=getappdata(resultfig,'stimparams');

for side=1:length(stimparams)
    VAT{side}=stimparams(side).VAT;
end
fiberthresh=stimparams.fiberthresh;


% prepare statvat exports once if needed.
if options.expstatvat.do % export statvat nifti images.
    tV=spm_vol([options.earoot,'templates',filesep,'bb.nii']);
    tnii=spm_read_vols(tV);
    tnii(:)=0;
    % generate mesh of hires MNI
    [x,y,z]=ind2sub(size(tnii),1:numel(tnii));
    templatecoords=[x;y;z;ones(1,length(x))]; clear x y z
    templatecoords=tV.mat*templatecoords;
    templatecoords=templatecoords(1:3,:)';
end


% set togglebuttons. small flag-variables are being defined in
% length(sides)xlength(labelatlases) dimensionality for the togglebuttons.

togglenames={'vaton','fibson','dcfibson','addfibson','labelson','captionson'};
for but=1:length(togglenames)
    
    eval([togglenames{but},'=getappdata(resultfig,''',togglenames{but},''');']);
    %fibson=getappdata(resultfig,'fibson');
    switch togglenames{but}
        case {'labelson','captionson','dcfibson'}
            expand=length(stimparams(1).labelatlas);
        otherwise
            expand=1;
    end
    if isempty(eval(togglenames{but}))
        eval([togglenames{but},'=repmat(1,expand,length(options.sides));']);
    end
    
    setappdata(resultfig,togglenames{but},eval(togglenames{but}));
    
end
clear expand



load([options.root,options.patientname,filesep,'ea_stats']);
try
    upriorvatlength=length(ea_stats.vat)+1;
    upriorftlength=length(ea_stats.ft)+1;

catch
    upriorvatlength=1;
    upriorftlength=1;
end


for side=1:length(VAT)
    if options.expstatvat.do;    thisvatnii=cell(length(options.expstatvat.vars),1); end
    for vat=1:length(VAT{side}.VAT)
        K(side).K{vat}=convhulln(VAT{side}.VAT{vat}+randn(size(VAT{side}.VAT{vat}))*0.000001); % create triangulation.
        
        % show vat
        
        PL.vatsurfs(side,vat)=trisurf(K(side).K{vat},VAT{side}.VAT{vat}(:,1),VAT{side}.VAT{vat}(:,2),VAT{side}.VAT{vat}(:,3),...
            abs(repmat(60,length(VAT{side}.VAT{vat}),1)...
            +randn(length(VAT{side}.VAT{vat}),1)*2)');
        
        
        
        % export vatstat if required:
        
        
        if options.expstatvat.do % export statvat nifti images.
            in=inhull(templatecoords,VAT{side}.VAT{vat},K(side).K{vat});
            for vatvar=1:length(options.expstatvat.vars)
                if isempty(thisvatnii{vatvar})
                    thisvatnii{vatvar}=tnii; % initialize on blank templatenii.
                end
                thisvatnii{vatvar}(in)=options.expstatvat.vars{vatvar}(options.expstatvat.pt); % set to clinical/regress variable score.
            end
        end
        
        % the following is some code required for
        % Web/JSON/BrainBrowser-Export.
        PL.vatfv(side,vat).vertices=[VAT{side}.VAT{vat}(:,1),VAT{side}.VAT{vat}(:,2),VAT{side}.VAT{vat}(:,3)];
        PL.vatfv(side,vat).faces=K(side).K{vat};
        PL.vatfv(side,vat).normals=get(PL.vatsurfs(side,vat),'Vertexnormals');
        PL.vatfv(side,vat).colors=repmat([1,0,0,0.7],size(PL.vatfv(side,vat).vertices,1),1);
        
        ea_spec_atlas(PL.vatsurfs(side,vat),'vat',jet,1);
        
        if options.writeoutstats
            
            
            
            if stimparams(side).U(vat)>0 % stimulation on in this VAT,
                load([options.root,options.patientname,filesep,'ea_stats']);
                try
                    priorvatlength=length(ea_stats.vat);
                catch
                    priorvatlength=0;
                end
                ea_stats.vat(priorvatlength+1).U=stimparams(side).U(vat);
                ea_stats.vat(priorvatlength+1).Im=stimparams(side).Im(vat);
                ea_stats.vat(priorvatlength+1).Contact=[side,vat];
                
                ea_stats.vat(priorvatlength+1).Side=side;
                iXYZ=getappdata(gcf,'iXYZ');
                ipixdim=getappdata(gcf,'ipixdim');
                
                for atlas=1:size(iXYZ,1)
                    
                    thisatl=iXYZ{atlas,side};
                    tpd=ipixdim{atlas,side};
                    if isempty(thisatl) % for midline or combined atlases, only the right side atlas is used.
                        thisatl=iXYZ{atlas,1};
                        tpd=ipixdim{atlas,1};
                    end
                    tpv=abs(tpd(1))*abs(tpd(2))*abs(tpd(3)); % volume of one voxel in mm^3.
                    
                    
                    ea_stats.vat(priorvatlength+1).AtlasIntersection(atlas)=sum(inhull(thisatl,VAT{side}.VAT{vat},K(side).K{vat}))*tpv;
                    ea_stats.vat(priorvatlength+1).nAtlasIntersection(atlas)=ea_stats.vat(priorvatlength+1).AtlasIntersection(atlas)/stimparams(1,side).volume(vat);
                end
                
                
                
                
                
                save([options.root,options.patientname,filesep,'ea_stats'],'ea_stats');
                
                
            end
            
        end
        
        
    end
    
    
    % export vatvar stats if needed:
    if options.expstatvat.do
        for vatvar=1:length(options.expstatvat.vars)
            mkdir([options.expstatvat.dir,'statvat_results']);
            switch side
                case 1
                    si='rh';
                case 2
                    si='lh';
            end
            tV.fname=[options.expstatvat.dir,'statvat_results',filesep,'s',num2str(options.expstatvat.pt),'_',si,'.nii'];
            spm_write_vol(tV,thisvatnii{vatvar});
        end
    end
    
    vatbutton(side)=uitoggletool(PL.ht,'CData',ea_get_icn('vat',options),'TooltipString','Volume of activated tissue','OnCallback',{@objvisible,PL.vatsurfs(side,:),resultfig,'vaton',[],side,1},'OffCallback',{@objvisible,PL.vatsurfs(side,:),resultfig,'vaton',[],side,0},'State',getstate(vaton(side)));

    
end



if stimparams(1).showfibers
    
    
    % load fiberset
    
    
    switch stimparams(1).usefiberset
        case 'Patient-specific DTI-Data'
            fs=load(fullfile(options.root,options.patientname,[options.prefs.FTR_normalized]));
        otherwise
            fs=load(fullfile(options.earoot,'fibers',[lower(stimparams(1).usefiberset),'.mat']));
    end
    fn = fieldnames(fs);
    
    eval(sprintf('normalized_fibers_mm = fs.%s;',fn{1}));
    if size(normalized_fibers_mm,1)>size(normalized_fibers_mm,2)
        normalized_fibers_mm=normalized_fibers_mm';
    end
    
    
    
    
    
    
    
    
    dispercent(0,'Selecting fibers');
    [idx,~]=cellfun(@size,normalized_fibers_mm);
    normalized_fibers_mm=cell2mat(normalized_fibers_mm');
    idxv=zeros(size(normalized_fibers_mm,1),1);
    lid=1; cnt=1;
    for id=idx
        idxv(lid:lid+id-1)=cnt;
        lid=lid+id;
        cnt=cnt+1;
        
    end
    cnt=1;
    
    % Select fibers
    
    
    maxvat=length(VAT);
    for side=1:maxvat
        for vat=1:length(VAT{side}.VAT)
            
            if stimparams(side).U(vat) % check if U ~= 0
                
                in=inhull(normalized_fibers_mm,VAT{side}.VAT{vat},K(side).K{vat})';
                
                selectedfibs{vat,side}=unique(idxv(in));
            end
        end
        dispercent(side/maxvat);
        
    end
    
    dispercent(1,'end');
    
    normalized_fibers_mm=mat2cell(normalized_fibers_mm,idx,3)';
    for side=options.sides
        sideselectedfibs{side}=unique(cell2mat(selectedfibs(:,side)));
        
        connectingfibs{side}=normalized_fibers_mm(sideselectedfibs{side});
    end
    selectedfibs=sideselectedfibs; clear sideselectedfibs
    
    dispercent(1,'end');
    
    
    
    % check which areas are connected to VAT by fibers:
    doubleconnectingfibs=cell(0);
    % check how many fibertrack stat-entries are there already (appending
    % to the analysis)
    load([options.root,options.patientname,filesep,'ea_stats']);
    
    try
        vll=length(ea_stats.ft);
    catch % if none, set to zeros.
        vll=0;
    end
    
    
    if stimparams(1).showconnectivities
        for la=1:length(stimparams(1).labelatlas)
            
            for side=options.sides
                
                todelete{la,side}=[];
                
                cnt=1; % reset cnt.
                
                %% extract areas connected by fibres.
                
                atlas=load_nii(fullfile(options.earoot,'templates','labeling',[stimparams(1).labelatlas{la},'.nii']));
                if options.writeoutpm && ~exist('pm','var')
                    pm=atlas;
                    pm.img(:)=0;
                end
                
                V=spm_vol(fullfile(options.earoot,'templates','labeling',[stimparams(1).labelatlas{la},'.nii']));
                aID = fopen(fullfile(options.earoot,'templates','labeling',[stimparams(1).labelatlas{la},'.txt']));
                atlas_lgnd=textscan(aID,'%d %s');
                allcareas=[];
                
                fibmax=length(connectingfibs{side});
                dispercent(0,'Gathering region information');
                for fib=1:fibmax
                    dispercent(fib/fibmax);
                    
                    thisfibendpoints=[connectingfibs{side}{fib}(1,1:3);connectingfibs{side}{fib}(end,1:3)];
                    thisfibendpoints=V.mat\[thisfibendpoints,ones(2,1)]'; % mm 2 vox
                    thisfibendpoints=double(thisfibendpoints(1:3,:));
                    
                    conareas=spm_sample_vol(V,thisfibendpoints(1,:),thisfibendpoints(2,:),thisfibendpoints(3,:),0);
                    if any(conareas)
                        doubleconnectingfibs{side}{la,cnt}=connectingfibs{side}{fib};
                        todelete{la,side}=[todelete{la,side},fib];
                        cnt=cnt+1;
                    end
                    if options.writeoutpm
                        try
                            pm.img(round(thisfibendpoints(1,1)),round(thisfibendpoints(2,1)),round(thisfibendpoints(3,1)))=...
                                pm.img(round(thisfibendpoints(1,1)),round(thisfibendpoints(2,1)),round(thisfibendpoints(3,1)))+1;
                            pm.img(round(thisfibendpoints(1,2)),round(thisfibendpoints(2,2)),round(thisfibendpoints(3,2)))=...
                                pm.img(round(thisfibendpoints(1,2)),round(thisfibendpoints(2,2)),round(thisfibendpoints(3,2)))+1;
                        end
                    end
                    allcareas=[allcareas,conareas];
                end
                allcareas=round(allcareas);
                dispercent(100,'end');
                
                atlength=length(atlas_lgnd{1});
                howmanyfibs{side}=zeros(atlength,1);
                tareas{side}=[];
                tcnt=1;
                for reg=1:atlength
                    howmanyfibs{side}(reg)=sum(allcareas==reg); % how many fibers connect VAT and anat. region.
                    if howmanyfibs{side}(reg)>=(fiberthresh(1))
                        tareas{side}(tcnt)=reg;
                        tcnt=tcnt+1;
                    end
                end
                tareas{side}=unique(tareas{side});
                
                
                % Write out connectivity stats
                if options.writeoutstats
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    ea_stats.ft(vll+side).fibercounts{la}=howmanyfibs{side};
                    
                    ea_stats.ft(vll+side).nfibercounts{la}=howmanyfibs{side}/sum(stimparams(1,side).volume);
                    ea_stats.ft(vll+side).labels{la}=atlas_lgnd{2};
                    % log which vats are used for this side of ft results.
                    cnt=1;
                    for vatsused=upriorvatlength:length(ea_stats.vat)
                        if ea_stats.vat(vatsused).Side==side
                            ea_stats.ft(vll+side).vatsused{la}(cnt)=vatsused;
                            cnt=cnt+1;
                        end
                    end
                    
                    
                end
            
            
            save([options.root,options.patientname,filesep,'ea_stats'],'ea_stats');
            
            
            
            clear allcareas conareas
            %% now show areas
            atlas.img=round(atlas.img);
            %tareas=1:4;
            if ~isempty(tareas{side})
                V=spm_vol(fullfile(options.earoot,'templates','labeling',[stimparams(1).labelatlas{la},'.nii']));
                for anatarea=1:length(tareas{side})
                    
                    [xx,yy,zz]=ind2sub(size(atlas.img),find(atlas.img==tareas{side}(anatarea)));
                    XYZ=[xx,yy,zz];
                    
                    XYZ=map_coords_proxy(XYZ,V);
                    %XYZ=XYZ';
                    if options.prefs.lhullmethod==0
                        k=convhulln(XYZ);
                    elseif options.prefs.lhullmethod==1
                        k=ea_concavehull(XYZ,1.5);
                        
                    end
                    [~,centroid]=kmeans(XYZ,1);
                    centroid=centroid(1,:);
                    
                    
                    %%
                    
                    if options.prefs.lhullmethod==2 % use isosurface
                        
                        bb=[0,0,0;size(atlas.img)];
                        
                        bb=map_coords_proxy(bb,V);
                        gv=cell(3,1);
                        for dim=1:3
                            gv{dim}=linspace(bb(1,dim),bb(2,dim),size(atlas.img,dim));
                        end
                        [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
                        
                        thisatlas=round(atlas.img);
                        thisatlas(thisatlas~=tareas{side}(anatarea))=0;
                        thisatlas(thisatlas==tareas{side}(anatarea))=1;
                        if options.prefs.lhullsmooth
                            thisatlas = smooth3(thisatlas,'gaussian',options.prefs.lhullsmooth);
                        end
                        fv=isosurface(X,Y,Z,permute(thisatlas,[2,1,3]),0.3);
                        if ischar(options.prefs.lhullsimplify)
                            
                            % get to 300 faces
                            simplify=300/length(fv.faces);
                            fv=reducepatch(fv,simplify);
                            
                        else
                            if options.prefs.lhullsimplify<1 && options.prefs.lhullsimplify>0
                                fv=reducepatch(fv,options.prefs.lhullsimplify);
                            end
                        end
                        % set cdata
                        
                        
                        
                        
                        if ~isfield(stimparams,'group')
                            cdat=abs(repmat(anatarea*(64/length(tareas{side})),length(fv.vertices),1)... % C-Data for surface
                                +randn(length(fv.vertices),1)*2)';
                        else % if more than one group is analyzed, coloring info will be off the group color.
                            RGB=zeros(1,1,3);
                            
                            RGB(:,:,1)=stimparams(1).groupcolors(stimparams(1).group,1);
                            RGB(:,:,2)=stimparams(1).groupcolors(stimparams(1).group,2);
                            RGB(:,:,3)=stimparams(1).groupcolors(stimparams(1).group,3);
                            
                            Rind=double(rgb2ind(RGB,jet));
                            cdat=abs(repmat(Rind,length(fv.vertices),1)... % C-Data for surface
                                +randn(length(fv.vertices),1)*2)';
                        end
                        
                        PL.regionsurfs(la,side,anatarea)=patch(fv,'CData',cdat,'FaceColor',[0.8 0.8 1.0],'facealpha',0.7,'EdgeColor','none','facelighting','phong');
                        
                    else
                        
                        
                        PL.regionsurfs(la,side,anatarea)=trisurf(k,XYZ(:,1),XYZ(:,2),XYZ(:,3),...
                            abs(repmat(anatarea*(64/length(tareas{side})),length(XYZ),1)...
                            +randn(length(XYZ),1)*0.1*length(tareas{side}))');
                    end
                    
                    
                    
                    
                    
                    %%
                    
                    
                    
                    
                    %% shading etc.
                    colorc=colornames(anatarea);
                    colorc=rgb(colorc);
                    ea_spec_atlas(PL.regionsurfs(la,side,anatarea),'labeling',jet,1);
                    
                    
                    %% put a label to it
                    thislabel=sub2space(atlas_lgnd{2}{atlas_lgnd{1}==tareas{side}(anatarea)});
                    PL.conlabels(la,side,anatarea)=text(centroid(1),centroid(2),centroid(3),thislabel,'VerticalAlignment','Baseline');
                end
                
            end
            end
                        clear tareas

        end
        
        
        
        % Write out probability map of fiber terminals
        if options.writeoutpm
            save_nii(pm,[options.root,options.patientname,filesep,'ea_pm','.nii']);
            
        end
        
    end
    
    
    
    % plot fibers that do connect to electrode VAT:
    for side=options.sides
        try % since this is only defined if using show_connectivities, too.
            alltodelete{side}=[];
            for la=1:length(stimparams(1).labelatlas)
                alltodelete{side}=[alltodelete{side},todelete{la,side}];
            end
            connectingfibs{side}(alltodelete{side})=[]; % clear doubleconnected fibers (will be plotted lateron) from vatconnected fibers.
        end
        
        if ~isempty(connectingfibs{side})
            fibmax=length(connectingfibs{side});
            dispercent(0,'Plotting fibers that connect to VAT of electrode');
            for fib=1:fibmax
                dispercent(fib/fibmax);
                %for segment=1:length(connectingfibs{fib})-1;
                connectingfibs{side}{la,fib}=connectingfibs{side}{la,fib}';
                
                
                if ~isfield(stimparams,'group')
                    connectingfibs{side}{la,fib}(4,:)=detcolor(connectingfibs{side}{la,fib}); % add coloring information to the 4th column.
                else % if more than one group is analyzed, coloring info will be off the group color.
                    RGB=zeros(1,1,3);
                    RGB(:,:,1)=stimparams(1).groupcolors(stimparams(1).group,1);
                    RGB(:,:,2)=stimparams(1).groupcolors(stimparams(1).group,2);
                    RGB(:,:,3)=stimparams(1).groupcolors(stimparams(1).group,3);
                    
                    connectingfibs{side}{la,fib}(4,:)=rgb2ind(RGB,jet);
                end
                
                
                
                
                for dim=1:4
                    thisfib(dim,:)=double(interp1q([1:size(connectingfibs{side}{fib},2)]',connectingfibs{side}{fib}(dim,:)',[1:0.1:size(connectingfibs{side}{fib},2)]')');
                end
                
                PL.fib_plots.fibs(side,fib)=surface([thisfib(1,:);thisfib(1,:)],...
                    [thisfib(2,:);thisfib(2,:)],...
                    [thisfib(3,:);thisfib(3,:)],...
                    [thisfib(4,:);thisfib(4,:)],'facecol','no','edgecol','interp','linew',1.5);
                
                % store for webexport
                
                jetlist=jet;
                
                try
                    PL.bbfibfv(fib).vertices=thisfib(1:3,:)';
                    PL.bbfibfv(fib).faces=[1:size(thisfib,2)-1;2:size(thisfib,2)]';
                    PL.bbfibfv(fib).normals=zeros(size(PL.bbfibfv(fib).vertices,1),3);
                    PL.bbfibfv(fib).colors=[squeeze(ind2rgb(round(thisfib(4,:)),jetlist)),repmat(0.7,size(thisfib,2),1)];
                end
                
                
                
                
                
                
                clear thisfib
                
            end
            dispercent(100,'end');
            
            set(PL.fib_plots.fibs(side,logical(PL.fib_plots.fibs(side,:))),'EdgeAlpha',0.05);
            
            
            
            try
                fiberbutton=uitoggletool(PL.ht,'CData',ea_get_icn('fibers_vat',options),'TooltipString','Fibers (Electrode only)','OnCallback',{@objvisible,PL.fib_plots.fibs(side,:),resultfig,'fibson',[],side,1},'OffCallback',{@objvisible,PL.fib_plots.fibs(side,:),resultfig,'fibson',[],side,0},'State',getstate(fibson(side)));
            end
        end
    end
    
    % plot fibers that connect to both a region of the labeling atlas and the electrode VAT:
    try cnt=length(PL.bbfibfv)+1; end
    for la=1:size(doubleconnectingfibs{side},1)
        for side=options.sides
            try
            fibmax=length(doubleconnectingfibs{side}(la,:));
            
            dispercent(0,'Plotting fibers that connect to both the VAT and a region within the labeling atlas');
            for fib=1:fibmax
                
                    dispercent(fib/fibmax);
                    doubleconnectingfibs{side}{la,fib}=doubleconnectingfibs{side}{la,fib}';
                    
                    if ~isfield(stimparams,'group')
                        
                        doubleconnectingfibs{side}{la,fib}(4,:)=detcolor(doubleconnectingfibs{side}{la,fib}); % add coloring information to the 4th column.
                        
                    else % if more than one group is analyzed, coloring info will be off the group color.
                        RGB=zeros(1,1,3);
                        RGB(:,:,1)=stimparams(1).groupcolors(stimparams(1).group,1);
                        RGB(:,:,2)=stimparams(1).groupcolors(stimparams(1).group,2);
                        RGB(:,:,3)=stimparams(1).groupcolors(stimparams(1).group,3);
                        
                        doubleconnectingfibs{side}{la,fib}(4,:)=rgb2ind(RGB,jet);
                    end
                    
                    for dim=1:4
                        thisfib(dim,:)=double(interp1q([1:size(doubleconnectingfibs{side}{la,fib},2)]',doubleconnectingfibs{side}{la,fib}(dim,:)',[1:0.1:size(doubleconnectingfibs{side}{la,fib},2)]')');
                    end
                    
                    PL.fib_plots.dcfibs(la,side,fib)=surface([thisfib(1,:);thisfib(1,:)],...
                        [thisfib(2,:);thisfib(2,:)],...
                        [thisfib(3,:);thisfib(3,:)],...
                        [thisfib(4,:);thisfib(4,:)],'facecol','no','edgecol','interp','linew',1.5);
                    
                    % store for webexport
                    jetlist=jet;
                    
                    try
                        PL.bbfibfv(cnt).vertices=thisfib(1:3,:)';
                        PL.bbfibfv(cnt).faces=[1:size(thisfib,2)-1;2:size(thisfib,2)]';
                        PL.bbfibfv(cnt).normals=zeros(size(PL.bbfibfv(cnt).vertices,1),3);
                        PL.bbfibfv(cnt).colors=[squeeze(ind2rgb(round(thisfib(4,:)),jetlist)),repmat(0.7,size(thisfib,2),1)];
                        
                        
                        cnt=cnt+1;
                    end
                    
                    
                    clear thisfib
                
                    
                end
            end
            dispercent(100,'end');
            
            try
            set(PL.fib_plots.dcfibs(la,side,logical(PL.fib_plots.dcfibs(la,side,:))),'EdgeAlpha',0.05);
            end
            
            try
                dcfiberbutton(la,side)=uitoggletool(PL.ht,'CData',ea_get_icn('fibers_both',options),'TooltipString','Fibers (Electrode and Labeling Atlas)','OnCallback',{@objvisible,PL.fib_plots.dcfibs(la,side,:),resultfig,'dcfibson',la,side,1},'OffCallback',{@objvisible,PL.fib_plots.dcfibs(la,side,:),resultfig,'dcfibson',la,side,0},'State',getstate(dcfibson(la,side)));
            
            end
            
            try
                regionbutton(la,side)=uitoggletool(PL.ht,'CData',ea_get_icn('connectivities',options),'TooltipString','Connected Regions','OnCallback',{@objvisible,PL.regionsurfs(la,side,:),resultfig,'labelson',la,side,1},'OffCallback',{@objvisible,PL.regionsurfs(la,side,:),resultfig,'labelson',la,side,0},'State',getstate(labelson(la,side)));
                captionbutton(la,side)=uitoggletool(PL.ht,'CData',ea_get_icn('labels',options),'TooltipString','Captions of Connected Regions','OnCallback',{@objvisible,PL.conlabels(la,side,:),resultfig,'captionson',la,side,1},'OffCallback',{@objvisible,PL.conlabels(la,side,:),resultfig,'captionson',la,side,0},'State',getstate(captionson(la,side)));
                
            end
        end
    end
    
  
    
end
% correct togglestates


if options.writeoutstats
    load([options.root,options.patientname,filesep,'ea_stats']);
    
    try
        ea_stats.vatanalyses(end+1).vatsused=upriorvatlength:length(ea_stats.vat);
    catch
        ea_stats.vatanalyses(1).vatsused=upriorvatlength:length(ea_stats.vat);
    end
    if stimparams(1).showfibers
        try
            ea_stats.vatanalyses(end).fibersused=upriorftlength:length(ea_stats.ft);
        catch
            ea_stats.vatanalyses(1).vatsused=upriorvatlength:length(ea_stats.vat);

        end
    end
    save([options.root,options.patientname,filesep,'ea_stats'],'ea_stats');
    
    
end

for side=options.sides
    if ~vaton(side)
        try
            objvisible([],[],PL.vatsurs(side,:),resultfig,'vaton',[],side,0)
        end
    end
    
    if ~fibson(side)
        try
            objvisible([],[],PL.fib_plots.fibs(side,:),resultfig,'fibson',[],side,0)
        end
    end
    
    
    
    for la=1:length(stimparams(1).labelatlas)
        
        if ~dcfibson(la,side)
            try
                objvisible([],[],PL.fib_plots.dcfibs(la,side,:),resultfig,'dcfibson',la,side,0)
            end
        end
        
        if ~labelson(la,side)
            try
                objvisible([],[],PL.regionsurfs(la,side,:),resultfig,'labelson',la,side,0)
            end
        end
        
        if ~captionson(la,side)
            try
                objvisible([],[],PL.conlabels(la,side,:),resultfig,'captionson',la,side,0)
            end
            
        end
    end
    
end


setappdata(resultfig,'PL',PL);


function indcol=detcolor(mat) % determine color based on traversing direction.

xyz=abs(diff(mat,1,2));
rgb=xyz/max(xyz(:));

rgb=[rgb,rgb(:,end)];
rgbim=zeros(1,size(rgb,2),3);
rgbim(1,:,:)=rgb';
indcol=double(rgb2ind(rgbim,jet));

function objvisible(hobj,ev,atls,resultfig,what,la,side,onoff)
% set visibility
set(atls, 'Visible', getstate(onoff));

% log visibility
tvalue=getappdata(resultfig,what);


if isempty(la)
    tvalue(side)=onoff;
else
    tvalue(la,side)=onoff;
end

setappdata(resultfig,what,tvalue);
%disp([atls,'visible clicked']);




function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);


function str=sub2space(str) % replaces subscores with spaces
str(str=='_')=' ';


function in = inhull(testpts,xyz,tess,tol)

% Copyright (c) 2009, John D'Errico
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% inhull: tests if a set of points are inside a convex hull
% usage: in = inhull(testpts,xyz)
% usage: in = inhull(testpts,xyz,tess)
% usage: in = inhull(testpts,xyz,tess,tol)
%
% arguments: (input)
%  testpts - nxp array to test, n data points, in p dimensions
%       If you have many points to test, it is most efficient to
%       call this function once with the entire set.
%
%  xyz - mxp array of vertices of the convex hull, as used by
%       convhulln.
%
%  tess - tessellation (or triangulation) generated by convhulln
%       If tess is left empty or not supplied, then it will be
%       generated.
%
%  tol - (OPTIONAL) tolerance on the tests for inclusion in the
%       convex hull. You can think of tol as the distance a point
%       may possibly lie outside the hull, and still be perceived
%       as on the surface of the hull. Because of numerical slop
%       nothing can ever be done exactly here. I might guess a
%       semi-intelligent value of tol to be
%
%         tol = 1.e-13*mean(abs(xyz(:)))
%
%       In higher dimensions, the numerical issues of floating
%       point arithmetic will probably suggest a larger value
%       of tol.
%
%       DEFAULT: tol = 0
%
% arguments: (output)
%  in  - nx1 logical vector
%        in(i) == 1 --> the i'th point was inside the convex hull.
%
% Example usage: The first point should be inside, the second out
%
%  xy = randn(20,2);
%  tess = convhulln(xy);
%  testpoints = [ 0 0; 10 10];
%  in = inhull(testpoints,xy,tess)
%
% in =
%      1
%      0
%
% A non-zero count of the number of degenerate simplexes in the hull
% will generate a warning (in 4 or more dimensions.) This warning
% may be disabled off with the command:
%
%   warning('off','inhull:degeneracy')
%
% See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 10/26/06

% get array sizes
% m points, p dimensions
p = size(xyz,2);
[n,c] = size(testpts);
if p ~= c
    error 'testpts and xyz must have the same number of columns'
end
if p < 2
    error 'Points must lie in at least a 2-d space.'
end

% was the convex hull supplied?
if (nargin<3) || isempty(tess)
    tess = convhulln(xyz);
end
[nt,c] = size(tess);
if c ~= p
    error 'tess array is incompatible with a dimension p space'
end

% was tol supplied?
if (nargin<4) || isempty(tol)
    tol = 0;
end

% build normal vectors
switch p
    case 2
        % really simple for 2-d
        nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * [0 1;-1 0];
        
        % Any degenerate edges?
        del = sqrt(sum(nrmls.^2,2));
        degenflag = (del<(max(del)*10*eps));
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate edges identified in the convex hull'])
            
            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
    case 3
        % use vectorized cross product for 3-d
        ab = xyz(tess(:,1),:) - xyz(tess(:,2),:);
        ac = xyz(tess(:,1),:) - xyz(tess(:,3),:);
        nrmls = cross(ab,ac,2);
        degenflag = false(nt,1);
    otherwise
        % slightly more work in higher dimensions,
        nrmls = zeros(nt,p);
        degenflag = false(nt,1);
        for i = 1:nt
            % just in case of a degeneracy
            % Note that bsxfun COULD be used in this line, but I have chosen to
            % not do so to maintain compatibility. This code is still used by
            % users of older releases.
            %  nullsp = null(bsxfun(@minus,xyz(tess(i,2:end),:),xyz(tess(i,1),:)))';
            nullsp = null(xyz(tess(i,2:end),:) - repmat(xyz(tess(i,1),:),p-1,1))';
            if size(nullsp,1)>1
                degenflag(i) = true;
                nrmls(i,:) = NaN;
            else
                nrmls(i,:) = nullsp;
            end
        end
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate simplexes identified in the convex hull'])
            
            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
end

% scale normal vectors to unit length
nrmllen = sqrt(sum(nrmls.^2,2));
% again, bsxfun COULD be employed here...
%  nrmls = bsxfun(@times,nrmls,1./nrmllen);
nrmls = nrmls.*repmat(1./nrmllen,1,p);

% center point in the hull
center = mean(xyz,1);

% any point in the plane of each simplex in the convex hull
a = xyz(tess(~degenflag,1),:);

% ensure the normals are pointing inwards
% this line too could employ bsxfun...
%  dp = sum(bsxfun(@minus,center,a).*nrmls,2);
dp = sum((repmat(center,nt,1) - a).*nrmls,2);
k = dp<0;
nrmls(k,:) = -nrmls(k,:);

% We want to test if:  dot((x - a),N) >= 0
% If so for all faces of the hull, then x is inside
% the hull. Change this to dot(x,N) >= dot(a,N)
aN = sum(nrmls.*a,2);

% test, be careful in case there are many points
in = false(n,1);

% if n is too large, we need to worry about the
% dot product grabbing huge chunks of memory.
memblock = 1e6;
blocks = max(1,floor(n/(memblock/nt)));
aNr = repmat(aN,1,length(1:blocks:n));
for i = 1:blocks
    j = i:blocks:n;
    if size(aNr,2) ~= length(j),
        aNr = repmat(aN,1,length(j));
    end
    in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
end


function  dispercent(varargin)
%
percent=round(varargin{1}*100);

if nargin==2
    if strcmp(varargin{2},'end')
        fprintf('\n')
        fprintf('\n')
        
        fprintf('\n')
        
    else
        fprintf(1,[varargin{2},':     ']);
        
        
    end
else
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
end




function str=getstate(val)
switch val
    case 1
        str='on';
    case 0
        str='off';
end


function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';
