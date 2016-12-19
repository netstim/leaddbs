% interface program for batch program
% Operations with streamline tracts (Mori Tracts).
% Author: Susanne Schnell
% PC 13.02.2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = dti_ftrop_ui(cmd, P)

switch cmd,
    case 'selectbyROI'
        autofname = '_selFTR.mat';

        %load fibers and ROI
        FTRStruct = ftrstruct_read(P.filename1{1});        
        ROI = maskstruct_read(P.roidef.roiname{1});
        
        %% transform fibers into space defined ROI
        t = inv(ROI.mrsProp.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   
        
        if P.seltype == 0,
            SelType = 'inArea';
        else
            SelType = 'endPointInArea'  ; %'inArea';
        end;
        
        % select fibers
        if isfield(P.roidef.mask,'masknumber')
            ROInames = maskstruct_query(ROI,'maskNames');
            if isfield(P.fibersubset,'fibernumber')
                FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', ROInames{P.roidef.mask.masknumber}), FiberNames{P.fibersubset.fibernumber}, 1);
            else
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), P.fibersubset.fibername, 1);
            end
        else
            if isfield(P.fibersubset,'fibernumber')
                FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), FiberNames{P.fibersubset.fibernumber}, 1);
            else
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), P.fibersubset.fibername, 1);
            end
        end
        if ~isempty(errStr)
            error(errStr)
        end
        if isempty(fibersNew)
            error('No fibers selected!');
            return;
        end
        roiNames= ftrstruct_query(FTRStruct, 'roiNames', []);
        if isfield(P.roidef.mask,'masknumber')
            roiNames{end + 1, 1}= ROInames{P.roidef.mask.masknumber};
        else
            roiNames{end + 1, 1}= P.roidef.mask.maskname;
        end
        fibersNew.roiName= roiNames;
        fibersNew.name= P.fibername;
        [FTRStruct, errStr]= ftrstruct_modify(FTRStruct, 'insertFiber', fibersNew);  

        %% transform fibers into fiber space 
        invt = inv(t);
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*invt(1:3,1:3)' + repmat(invt(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);           
        
        if ~isempty(errStr)
            error(errStr)
        end
        
        
    case 'selectbyOverlapROI'
        autofname = '_selFTR.mat';

        %load fibers and ROI
        FTRStruct = ftrstruct_read(P.filename1{1});        
        ROI = maskstruct_read(P.roidef.roiname{1});
        
        %% transform fibers into space defined ROI
        t = inv(ROI.mrsProp.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   
        SelType = 'inArea';
                
        fuzzy.sigma = 1;
        fuzzy.thres = P.overlapratio;
        fuzzy.relative = true;
        
        % select fibers
        if isfield(P.roidef.mask,'masknumber')
            ROInames = maskstruct_query(ROI,'maskNames');
            if isfield(P.fibersubset,'fibernumber')
                FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', ROInames{P.roidef.mask.masknumber}), FiberNames{P.fibersubset.fibernumber}, 1,[],fuzzy);
            else
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), P.fibersubset.fibername, 1,[],fuzzy);
            end
        else
            if isfield(P.fibersubset,'fibernumber')
                FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), FiberNames{P.fibersubset.fibernumber}, 1,[],fuzzy);
            else
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,SelType, maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), P.fibersubset.fibername, 1,[],fuzzy);
            end
        end
        if ~isempty(errStr)
            error(errStr)
        end
        if isempty(fibersNew)
            error('No fibers selected!');
            return;
        end
        roiNames= ftrstruct_query(FTRStruct, 'roiNames', []);
        if isfield(P.roidef.mask,'masknumber')
            roiNames{end + 1, 1}= ROInames{P.roidef.mask.masknumber};
        else
            roiNames{end + 1, 1}= P.roidef.mask.maskname;
        end
        fibersNew.roiName= roiNames;
        fibersNew.name= P.fibername;
        [FTRStruct, errStr]= ftrstruct_modify(FTRStruct, 'insertFiber', fibersNew);  

        %% transform fibers into fiber space 
        invt = inv(t);
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*invt(1:3,1:3)' + repmat(invt(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);           
        
        if ~isempty(errStr)
            error(errStr)
        end
        
        
        
        
        
        
        
    case 'eliminatebyROI'
        autofname = '_elimFTR.mat';
        %eliminate fibers
        FTRStruct = ftrstruct_read(P.filename1{1});
        ROI = maskstruct_read(P.roidef.roiname{1});
        
        %% transform fibers into space defined ROI
        t = inv(ROI.mrsProp.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   
        
        
        if isfield(P.roidef.mask,'masknumber')
            ROInames = maskstruct_query(ROI,'maskNames');
            if isfield(P.fibersubset,'fibernumber')
                FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,'inArea',maskstruct_query(ROI, 'getMRMask', ROInames{P.roidef.mask.masknumber}), FiberNames{P.fibersubset.fibernumber}, 0);
            else
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,'inArea',maskstruct_query(ROI, 'getMRMask', ROInames{P.roidef.mask.masknumber}), P.fibersubset.fibername, 0);
            end
        else
            if isfield(P.fibersubset,'fibernumber')
                FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,'inArea',maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), FiberNames{P.fibersubset.fibernumber}, 0);
            else
                [fibersNew, errStr]= ftrstruct_query(FTRStruct,'inArea',maskstruct_query(ROI, 'getMRMask', P.roidef.mask.maskname), P.fibersubset.fibername, 0);
            end
        end
        if isempty(fibersNew)
            return;
        end
        if ~isempty(errStr)
            error(errStr)
        end
        roiNames= ftrstruct_query(FTRStruct, 'roiNames', []);
        if isfield(P.roidef.mask,'masknumber')
            roiNames{end + 1, 1}= ROInames{P.roidef.mask.masknumber};
        else
            roiNames{end + 1, 1}= P.roidef.mask.maskname;
        end
        fibersNew.roiName= roiNames;
        fibersNew.name= P.fibername;
        [FTRStruct, errStr]= ftrstruct_modify(FTRStruct, 'insertFiber', fibersNew);

        %% transform fibers into fiber space 
        invt = inv(t);
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*invt(1:3,1:3)' + repmat(invt(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);           
        
        
        if ~isempty(errStr)
            error(errStr)
        end
    case 'visitMapext'
        autofname = '_visitMap.mat';
        
        %load data
        FTRStruct = ftrstruct_read(P.filename1{1});
        
        
        ref = P.filename3{1};
        [pathname finame ext] = fileparts(ref);
        if strcmp(ext,'.mat');
            mr = load(ref);
            if isfield(mr,'mrStruct')
                mr = mr.mrStruct;
            elseif isfield(mr,'b0_image_struc')
                mr = mr.b0_image_struc;
            end;
        else
          mr = nifti_to_mrstruct('volume',{ref});   
        end;
        
        
        
        
        vox = mr.vox;
        t = inv(mr.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   
        if isfield(P.fibersubset,'fibernumber')
            fnum = P.fibersubset.fibernumber;
        else
            FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
            fnum = find(cellfun(@(x) strcmp(x,P.fibersubset.fibername),FiberNames));
        end;
        if fnum == 1,
            fnum = [];
        else
            fnum = fnum -1;
        end;

        reparam_step = 0.9*min(vox(1:3))/P.oversampling; 
        for k = 1:length(FTRStruct.curveSegCell),
            FTRStruct.curveSegCell{k} = reparametrize_arclen(single(FTRStruct.curveSegCell{k})',double(reparam_step))';
        end;        
        
        [fdrgb,fd,ep] = ftr2FDmaps(FTRStruct,size(mr.dataAy),fnum,P.oversampling);
        
        if strcmp(P.newfilename,'.mat') || strcmp(P.newfilename,'.nii')            
           [pathname finame ] = fileparts(P.filename1{1});
           outname = fullfile(pathname,finame);
           ext = P.newfilename;
        else
           [pathname outname ext] = fileparts(P.newfilename);
           outname = fullfile(pathname,outname);
        end;
        
        fdrgbname = [outname '_fdrgb'];
        fdname = [outname '_fd'];
        epname = [outname '_ep'];
        
        mr.vox = mr.vox/P.oversampling;
        mI = eye(4); mI(1:3,4) = -1;
        sC = eye(3)/P.oversampling; sC(4,4) = 1;
        mr.edges = mr.edges*inv(mI)*sC*mI;

        if strcmp(ext,'.mat')
           mr.memoryType = 'series3D';
           mr.dataAy = fdrgb; mrstruct_write(mr,[fdrgbname '.mat']);    out.fdrgb = {[fdrgbname '.mat']};
           mr.memoryType = 'volume';
           mr.dataAy = fd; mrstruct_write(mr,[fdname '.mat']);       out.fd = {[fdname '.mat']};
           mr.memoryType = 'volume';
           mr.dataAy = ep;  mrstruct_write(mr,[epname '.mat']);       out.ep = {[epname '.mat']};
        elseif strcmp(ext,'.nii')
           mr.memoryType = 'series3D';
           mr.dataAy = fdrgb; [res, errStr] = mrstruct_to_nifti(mr,[fdrgbname '.nii'],'float32'); out.fdrgb = {[fdrgbname '.nii']};
           if ~isempty(errStr)
               error(errStr)
           end
           mr.memoryType = 'volume';
           mr.dataAy = fd; [res, errStr] = mrstruct_to_nifti(mr,[fdname '.nii'],'float32');  out.fd = {[fdname '.nii']};
           if ~isempty(errStr)
               error(errStr)
           end
           mr.memoryType = 'volume';
           mr.dataAy = ep; [res, errStr] = mrstruct_to_nifti(mr,[epname '.nii'],'float32');  out.ep = {[epname '.nii']};
           if ~isempty(errStr)
               error(errStr)
           end
        end;
        
         
    case 'visitMap'
        autofname = '_visitMap.mat';
        
        %load data
        FTRStruct = ftrstruct_read(P.filename1{1});
        DTD = dtdstruct_read(P.filename2{1});        
        t = inv(DTD.b0_image_struc.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   
                
        
        if isfield(P.fibersubset,'fibernumber')
            FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
            %visit Mask
            [map, errStr]= ftrstruct_query(FTRStruct, 'getVisitMap', FiberNames{P.fibersubset.fibernumber},size(DTD.b0_image_struc.dataAy));
        else
            [map, errStr]= ftrstruct_query(FTRStruct, 'getVisitMap', P.fibersubset.fibername,size(DTD.b0_image_struc.dataAy));
        end
        if ~isempty(errStr)
            error(errStr)
        end
    case 'visitMask'
        autofname = '_visitMask.mat';

        %load data
        FTRStruct = ftrstruct_read(P.filename1{1});
        DTD = dtdstruct_read(P.filename2{1});

        t = inv(DTD.b0_image_struc.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   
        
        
        if isfield(P.fibersubset,'fibernumber')
            FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
            %visit Mask
            [mask, errStr]= ftrstruct_query(FTRStruct, 'getVisitMask', FiberNames{P.fibersubset.fibernumber},size(DTD.b0_image_struc.dataAy),P.thresh);
        else
            %visit Mask
            [mask, errStr]= ftrstruct_query(FTRStruct, 'getVisitMask',P.fibersubset.fibername,size(DTD.b0_image_struc.dataAy),P.thresh);
        end
        if ~isempty(errStr)
            error(errStr)
        end
    case 'endpointMask'
        autofname = '_endMask.mat';
        
        % load data
        FTRStruct = ftrstruct_read(P.filename1{1});
        DTD = dtdstruct_read(P.filename2{1});
        
        t = inv(DTD.b0_image_struc.edges)*FTRStruct.hMatrix;        
        FTRStruct.curveSegCell = cellfun(@(x) 1+ (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) , FTRStruct.curveSegCell,'uniformoutput',false);   

        
        if isfield(P.fibersubset,'fibernumber')
            FiberNames = ftrstruct_query(FTRStruct, 'fiberNames');
            % endpointMask
            [map, errStr]= ftrstruct_query(FTRStruct, 'getEndPointMap', FiberNames{P.fibersubset.fibernumber},size(DTD.b0_image_struc.dataAy));
        else
            [map, errStr]= ftrstruct_query(FTRStruct, 'getEndPointMap', P.fibersubset.fibername,size(DTD.b0_image_struc.dataAy));
        end
        if ~isempty(errStr)
            error(errStr)
        end
    case 'deformFTR'

        autofname = '_defFTR.mat';
        try
            FTRStruct = ftr_warp( P.filedefi{1},P.filedef{1}, P.filename1{1});
        catch        
            error(lasterr);
        end     
        
    case 'tractStats'
        
        FTRStruct = ftrstruct_read(P.filename1{1});
        mrProp = mrPropFromFTR(FTRStruct);
        con = getContrasts(P,mrProp);
        
        contrasts = {};
        if not(isempty(con)),
            contrasts = {con.mr};
            contrastnames = {con.name};
        end;
                        
        bundles = stat_along_tract(FTRStruct,contrasts);
        
        
        [p n e] = fileparts(P.filename1{1});
        txtfile = fullfile(p,[P.newfilename2 '.txt']);
        fid = fopen(txtfile,'w');
                
     
        fprintf(fid,'name fibcnt mean_length sdev_length max_length min_length prctile05_length prctile95_length\n');
        for k = 1:length(bundles),
            fprintf(fid,'%s\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n',bundles(k).name ,bundles(k).fibcnt ,bundles(k).mean_length ,bundles(k).sdev_length ,bundles(k).max_length ,bundles(k).min_length ,bundles(k).prctile05_length ,bundles(k).prctile95_length);
        end
        fprintf(fid,'\n');
        if isfield(bundles(1),'contrast')
            fprintf(fid,'bundle contrast mean sdev min max prctile05 prctile95\n');
            for k = 1:length(bundles),
                for j = 1:length(bundles(k).contrast)
                    c = bundles(k).contrast(j);
                     bundles(k).contrast(j).name = con(j).name;
                     fprintf(fid,'%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n',bundles(k).name,con(j).name,c.mean,c.sdev,c.min,c.max,c.prctile05,c.prctile95);
                end;
            end;
            fprintf(fid,'\n');        
        end;
        fclose(fid);
        
        matfile = fullfile(p,[P.newfilename2 '.mat']);
        if exist('contrastnames'),
            save(matfile,'bundles','contrastnames');
        else
            save(matfile,'bundles');
        end;
        out.files{1} = matfile;
        
    case 'tractStatsROI'  
                
        % get FTR
        ftr = ftrstruct_read(P.filename1{1});
        mrProp = mrPropFromFTR(ftr);
        
        
        % get Contrasts for tractometry
        con = getContrasts(P,mrProp);        
            
        % get ROIs
        [p n e] = fileparts(P.filename5{1});
        mastru = false;
        if strcmp(e,'.mat'),
            dat = load(P.filename5{1});
            if maskstruct_istype(dat),
                roi = dat;
                mastru = true;
            elseif isfield(dat,'mrStruct'),
                roi = dat;
            end;        
        elseif strcmp(e,'.nii') | strcmp(e,'.hdr')                
             resliceFlags.mask= 0;
             resliceFlags.mean= 0;
             resliceFlags.interp= 0; 
             resliceFlags.which=1;
             roi = nifti_to_mrstruct('volume', {P.filename5{1}},mrProp,pwd,resliceFlags);             
        end;

        % selection type
        if P.seltype == 0,
            fibQuery = @(ftr,mask,i,fuzzy) ftrstruct_query(ftr,'inArea',mask,i,1,[],fuzzy);     
        else
            fibQuery = @(ftr,mask,i,fuzzy) ftrstruct_query(ftr,'endPointInAreaFuzzy',mask,i,fuzzy);     
        end;
        
        fuzzy = [];
        if isfield(P.fuzzysel,'fuzzy'),
            fuzzy.sigma = P.fuzzysel.fuzzy.sigmafuz;
            fuzzy.thres = P.fuzzysel.fuzzy.threshfuz;
        end;
        
        
        if not(mastru),
            if isempty(P.roinumber),
                roiidx = unique(roi.dataAy(roi.dataAy(:)>0));
            else
                roiidx = P.roinumber;
            end;
        else
            if isempty(P.roinumber),                
                roiidx = 1:length(roi.maskCell);
            else
                roiidx = P.roinumber;
            end;
        end;
        Nrois = length(roiidx);
        
        
        contrasts = {};
        if not(isempty(con)),
            contrasts = {con.mr};
            contrastnames = {con.name};
        end;
        
        for k = 1:Nrois
            
            mask = mrProp;
            if not(mastru),                
                 mask.dataAy = zeros(size(mask.dataAy));
                 mask.dataAy(roi.dataAy==roiidx(k)) = 1;
                 maskname = ['ROI' num2str(roiidx(k))];
            else
                 mask = roi.mrsProp;
                 mask.dataAy = roi.maskCell{roiidx(k)};
                 maskname = roi.maskNamesCell{roiidx(k)};
            end;
            ftr.fiber = {};
            fibsA =  fibQuery(ftr,mask,[],fuzzy);  %ftrstruct_query(ftr,SelType,mask,[],1,[],fuzzy);     
            bundles(k,k) = stat_along_tract([]);
            if not(isempty(fibsA)),
                fibsA.name = maskname;
                ftr.fiber = {fibsA};
                bundles(k,k) = stat_along_tract(ftr,contrasts);

                for j = k+1:Nrois
                    fprintf('.');
                    if not(mastru),                
                         mask.dataAy = zeros(size(mask.dataAy));
                         mask.dataAy(roi.dataAy==roiidx(j)) = 1;
                         maskname2 = ['ROI' num2str(roiidx(j))];
                    else
                         mask = roi.mrsProp;
                         mask.dataAy = roi.maskCell{roiidx(j)};
                         maskname2 = roi.maskNamesCell{roiidx(j)};
                    end;

                    fibs = fibQuery(ftr,mask,1,fuzzy); %ftrstruct_query(ftr,SelType,mask,1,1,[],fuzzy);
                    bundles(k,j) = stat_along_tract([]);
                    if not(isempty(fibs)),
                        fibs.name = [maskname '/' maskname2];
                        ftrn = ftr;
                        ftrn.fiber = {fibs};
                        bundles(k,j) = stat_along_tract(ftrn,contrasts);
                    end;

                    bundles(j,k) = bundles(k,j);

                end   
            end;
            fprintf('\n');
        end
    
        res.bundles = bundles;
        bundles = arrayfun(@(x) setfield(x,'fibcnt',emptyToZero(x.fibcnt)),bundles);
        res.Cmat = reshape([bundles(:).fibcnt],[Nrois Nrois]);
        res.totnumfibs = length(ftr.curveSegCell);

        matfile = fullfile(p,[P.newfilename2 '.mat']);        
        if exist('contrastnames'),
            res.contrastnames = contrastnames;
        else
            save(matfile,'bundles');
        end;                   
        save(matfile,'-struct','res');
        
        out.files{1} = matfile;
        
end


if not(strcmp(cmd,'visitMapext')) & not(strcmp(cmd,'tractStats'))  & not(strcmp(cmd,'tractStatsROI')),
    % newfilename
    if strcmp(P.newfilename,'.mat') || isempty(P.newfilename)
        [path,name] = fileparts(P.filename1{1});
        filename = fullfile(path, [name(1:end-4), autofname]);
    else
        [path,name,ext] = fileparts(P.filename1{1});
        [newpath,name,ext] = fileparts(P.newfilename);
        if isempty(newpath)
            filename = fullfile(path,[name ext]);
        else
            filename = fullfile(newpath,[name ext]);
        end
    end

    %save result
    if strcmp(cmd,'selectbyROI') || strcmp(cmd,'selectbyOverlapROI') ||strcmp(cmd ,'eliminatebyROI') || strcmp(cmd ,'deformFTR')
        [res, errStr] = ftrstruct_write(FTRStruct,filename);
        if ~isempty(errStr)
            error(errStr)
        end
    elseif strcmp(cmd,'visitMap')
        map = mrstruct_init('volume',map,DTD.b0_image_struc);
        mrstruct_write(map,filename);
    elseif strcmp(cmd,'visitMask')
        [res, errStr] = maskstruct_write(mask,filename);
        if ~isempty(errStr)
            error(errStr)
        end
    else
        [mask, errStr] = maskstruct_init(DTD.b0_image_struc);
        if ~isempty(errStr)
            error(errStr)
        end
        [mask, errStr] = maskstruct_modify(mask,'createMask','endPointMask');
        if ~isempty(errStr)
            error(errStr)
        end
        [mask, errStr] = maskstruct_modify(mask,'setMask',map,'endPointMask');
        if ~isempty(errStr)
            error(errStr)
        end
        [res, errStr] = maskstruct_write(mask,filename);
        if ~isempty(errStr)
            error(errStr)
        end
    end
    out.files{1} = filename;
end;



%%%% function to get gather contrasts for tractometry

function con = getContrasts(P,mrprop)
          
        if isfield(P.contrastsel,'nocontrast'),
            con = [];
            return;
        end;

        cnt = 1;
        for k = 1:length(P.contrastsel.filename4);
            [p n e] = fileparts(P.contrastsel.filename4{k});
            if strcmp(e,'.mat'),
                dat = load(P.contrastsel.filename4{k});
                if dtdstruct_istype(dat),
                    con(cnt).name = 'FA';
                    con(cnt).mr = dtdstruct_query(dat,'getFA');
                    cnt = cnt + 1;
                    con(cnt).name = 'Trace';
                    con(cnt).mr = dtdstruct_query(dat,'getTrace');
                    con(cnt).mr.dataAy = con(cnt).mr.dataAy*10^3;
                    cnt = cnt + 1;
                    con(cnt).name = 'EigAx';
                    con(cnt).mr = dtdstruct_query(dat,'getEigVal1');
                    con(cnt).mr.dataAy = con(cnt).mr.dataAy*10^3;
                    cnt = cnt + 1;
                    con(cnt).name = 'EigRad';
                    con(cnt).mr = dtdstruct_query(dat,'getEigVal2');
                    ev3 =  dtdstruct_query(dat,'getEigVal3');
                    con(cnt).mr.dataAy = (con(cnt).mr.dataAy + ev3.dataAy)*0.5;
                    con(cnt).mr.dataAy = con(cnt).mr.dataAy*10^3;
                    cnt = cnt + 1;
                    
                elseif isfield(dat,'mrStruct'),
                    con(cnt).name = n;
                    con(cnt).mr = dat.mrStruct;
                    cnt = cnt + 1;                
                end;
                
            elseif strcmp(e,'.nii') | strcmp(e,'.hdr')                
                 mrStruct = nifti_to_mrstruct('volume', {P.contrastsel.filename4{k}},mrprop);
                 con(cnt).name = n;
                 con(cnt).mr = mrStruct;
                 cnt = cnt + 1;
            end;
            
            
        end;

function x = emptyToZero(x)        
    if isempty(x)
        x = 0;
    end;
        


 function mrProp = mrPropFromFTR(ftr)
        mrProp = mrstruct_init;
        mrProp.edges = ftr.hMatrix;
        mrProp.vox = ftr.vox;
        sz =  ceil(max( cat(1,ftr.curveSegCell{:})));
        mrProp.dataAy = zeros(sz);
        mrProp.memoryType = 'volume';
        mrProp.dim3 = 'size_z';







