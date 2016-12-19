% interface program for batch program
% Operations with FODs
% Author: Marco Reisert
% PC 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = dti_ftrop_ui(cmd, P)

switch cmd,
    
    case 'PItrails'
         
        picon = load(P.filenamePIcon{1});
        
        [pa fi] = fileparts(P.filenamePIcon{1});
        fi = fi(1:end-3);
        
        roip = P.roipairs;
        
        template = picon.maps.mrstructProp;
        template.dim4 = 'unused';
        template.user = [];
        maps = picon.maps;
        Im = maps.InterpolationMat;
        crop = maps.cropRegion;
        sz = [crop(2)-crop(1)+1 crop(4)-crop(3)+1 crop(6)-crop(5)+1];
        sz(4) = size(maps.mrstructProp.user.bDir,2)*2;
        for k = 1:size(roip,2),            
            t1 = reshape(Im*double(maps.data(:,roip(1,k))),sz);
            t2 = reshape(Im*double(maps.data(:,roip(2,k))),sz);
            map = sum(t1.*circshift(t2,[1 1 1 sz(4)/2]),4);
            template.memoryType = 'volume';
            template.dataAy = zeros(maps.size);
            template.dataAy(crop(1):crop(2),crop(3):crop(4),crop(5):crop(6)) = map;
            fname{k} = fullfile(pa,[fi sprintf('roi%i-%i.nii',roip(1,k),roip(2,k))]);
            mrstruct_to_nifti(template, fname{k}, 'float32');            
        end;
            
        out.mapnames = fname;
    
    case 'PItrack'
    
        fodname = P.filename2{1};
        roiname =  P.filenameROIs{1};
        wmname = P.roiname{1};
        thres = P.thresh;
        
        
        roisubset = P.roisubset;
        
        fod = mrstruct_read(fodname);
        
        resliceFlags.mask= 0;
        resliceFlags.mean= 0;
        resliceFlags.interp= 0; 
        resliceFlags.which=1;
        
        roi = loadsegmentation(roiname,1,fod,resliceFlags);
        WM = loadsegmentation(wmname,1,fod,resliceFlags);
        
        thres_toMax = P.thres_toMax;
        
        P.saveftr = not(isfield(P.streamftr,'noftr'));
        P.saveftrpairs = not(isfield(P.streamftrpaired,'noftr'));
        
        
        [fods mask crop] = prepFODdata(fod.dataAy,fod.user.bDir,WM.dataAy>thres,P.sh,thres_toMax,P.normtype);
        fod = rmfield(fod,'dataAy'); % save memory
        rois = roi.dataAy(crop(1):crop(2),crop(3):crop(4),crop(5):crop(6));       
        
        result = PItracking(fods,mask,[fod.user.bDir -fod.user.bDir],rois,P,roisubset);
        result.maps.mrstructProp = fod;
        result.maps.cropRegion = crop;
        result.maps.size = size(roi.dataAy);
        result.roifilename = roiname;
        result.roiidx = roisubset;
        
        if P.savePImaps == 0,
            result = rmfield(result,'maps');
        end;
        
         [p n e] = fileparts(P.filename2{1});
         if strcmpi('fod',n(end-2:end)),
             n = n(1:end-3);
         end;                

        
        if P.saveftr
            %%
             ftr = ftrstruct_init;
             ftr.hMatrix = fod.edges;
             ftr.vox = fod.vox;
             ftr.curveSegCell = cat(2,result.fibs{:})';
             ftr.curveD = cellfun(@(x) [x(:,4)'; log(eps+abs(x(:,4)'))], ftr.curveSegCell,'uniformoutput',false);
             ftr.curveSegCell = cellfun(@(x) x(:,[2 1 3]) + repmat([crop(1) crop(3) crop(5)],size(x,1),1), ftr.curveSegCell,'uniformoutput',false);
             runf = 1;
             for k = 1:length(result.fibs),
                 ftr.fiber{k}.name = ['ROI_' num2str(k)];
                 ftr.fiber{k}.curveID = runf:(runf+length(result.fibs{k})-1);
                 ftr.fiber{k}.roiName = ['ROI_' num2str(k)];
                 ftr.fiber{k}.user = [];
                 runf = runf + length(result.fibs{k});
             end;
             ftr.fiber = ftr.fiber';
             ftr.connectCell = arrayfun(@(x)x ,(1:length(ftr.curveSegCell))','uniformoutput',false);
             
             outnameftr = buildoutname(P.streamftr.outftr.newfilenamemat,p,n,'FTR','.mat');                          
             save(outnameftr,'-struct','ftr');
             
        end;
        if P.saveftrpairs
             %%
             ftr = ftrstruct_init;
             ftr.hMatrix = fod.edges;
             ftr.vox = fod.vox;
             ftr.curveSegCell = cat(2,result.fibsPaired{:})';
             ftr.curveD = cellfun(@(x) [x(:,4)'; log(eps+abs(x(:,4)'))], ftr.curveSegCell,'uniformoutput',false);
             ftr.curveSegCell = cellfun(@(x) x(:,[2 1 3]) + repmat([crop(1) crop(3) crop(5)]-2,size(x,1),1), ftr.curveSegCell,'uniformoutput',false);
             runf = 1;
             cnt = 1;
             for k = 1:length(result.fibsPaired),
                 for j = 1:length(result.fibsPaired),
                     if k ~= j,
                         fibp = result.fibsPaired{j,k};
                         if not(isempty(fibp)),
                             ftr.fiber{cnt}.name = ['ROI_' num2str(k) '_' num2str(j)];
                             ftr.fiber{cnt}.curveID = runf:(runf+length(fibp)-1);
                             ftr.fiber{cnt}.roiName = ['ROI_' num2str(k) '_' num2str(j)];
                             ftr.fiber{cnt}.user = [];
                             runf = runf + length(fibp);
                             cnt = cnt + 1;
                         end;
                     end;
                 end;
             end;
             ftr.fiber = ftr.fiber';
             ftr.connectCell = arrayfun(@(x)x ,(1:length(ftr.curveSegCell))','uniformoutput',false);
             
             outnameftrpair = buildoutname(P.streamftr.outftr.newfilenamemat,p,n,'paired_FTR','.mat');                          
             save(outnameftrpair,'-struct','ftr');
                          
             
        end;
        
        
        savename = buildoutname(P.newfilenamemat,p,n,'PICON','.mat');                   
        save(savename,'-struct','result');
        
        out.PIinfo = {savename};
        
        if P.saveftr,
            out.PIinfo = [out.PIinfo {outnameftr}];
        end;
        
        if P.saveftrpairs,
            out.PIinfo = [out.PIinfo {outnameftrpair}];
        end;
        
        
case 'PXtrack'
    
        fodname = P.filename2{1};
        roiname =  P.filenameROIs{1};
        
        roisubset = P.roisubset;
        
        fod = mrstruct_read(fodname);

        
        if strcmp(P.newfilenamemat,'.mat')
            [p n e] = fileparts(P.filename2{1});
            n = n(1:end-4);
            savename = fullfile(p,[n '_PXCON.mat']);
        else
            savename = P.newfilenamemat;
        end;                
        
        resliceFlags.mask= 0;
        resliceFlags.mean= 0;
        resliceFlags.interp= 0; 
        resliceFlags.which=1;
        
        roi = loadsegmentation(roiname,1,fod,resliceFlags);
        usermask = sum(fod.dataAy,4)>0;
        P.rev = false;
        [cm_pt ROIsize Nfac] = probtrax_roi(repmat(fod.dataAy,[1 1 1 2]),[fod.user.bDir -fod.user.bDir]',usermask,roi.dataAy,P,roisubset);            
        
        result.roiidx = roisubset;
        result.Cmat = cm_pt;
        result.ROIsize = ROIsize;
        result.Nfac = Nfac;
        result.params = P;
        
        save(savename,'-struct','result');
        

        
        out.PXinfo = {savename};
        
    case 'doL1CSD'
      
        hrname = P.filename1{1};
        wmname = P.roidef.roiname{1};
        outname = P.newfilename;
       

        params.csdL1.maxit = P.maxit;
        params.csdL1.lambda = P.lambdaL1;
        params.csdL1.powsm = P.powsm;
        params.csdL1.numdir = P.numdir;
        params.csdL1.D = [P.Dax P.Drad];


        [mrfod] = computeCSDL1(hrname,wmname,params);

        [p n e] = fileparts(hrname);
        if strcmpi('hardi',n(end-4:end)),
            n = n(1:end-5);
        end;
                
        outname = buildoutname(outname,p,n,'FOD','.mat');
            
        mrstruct_write(mrfod,outname);
        
        out.FOD = {outname};
     case 'doCSDTournier'
         
        hrname = P.filename1{1};
        outname = P.newfilename;
        shcut = P.shcutoff;
        numdirs = P.numdirsout;
        Dax = P.Dax;
        wmname = P.roidef.roiname{1};
        thres = P.roidef.thresh;
        mr = mrstruct_read(hrname);
        
        wm = loadsegmentation(wmname,1,mr);
        mrfod = computeCSDSH(mr,wm.dataAy>thres,Dax,shcut,numdirs);

        [p n e] = fileparts(hrname);
        if strcmpi('hardi',n(end-4:end)),
            n = n(1:end-5);
        end;
                
        outname = buildoutname(outname,p,n,'FOD','.mat');
            
        mrstruct_write(mrfod,outname);
        
        out.FOD = {outname};

     case 'doCSDFC'
         
        hrname = P.filename1{1};
        outname = P.newfilename;
        numdirs = P.numoutdir;
        wmname = P.roidef.roiname{1};
        thres = P.roidef.thresh;
        numit = P.numit;
        numround = P.numround;
        Daxial = P.Daxial;
        lambda = P.lambda;
        alpha = P.alpha;
        purefc_sym = P.purefc_sym;
        purefc_asym = P.purefc_asym;
        iso_sym = P.iso_sym;
        iso_asym = P.iso_asym;
        lambda_lb = P.lambda_lb;
        gamma_sym = P.gamma_sym;
        gamma_asym = P.gamma_asym;
        symout = P.symout;
        
        
        mr = mrstruct_read(hrname);        
        wm = loadsegmentation(wmname,1,mr);
        
        if not(isfield(mr.user,'bfactor')),
            mr.user.bfactor = 1000;
        end;
                
        if size(mr.user.bDir,1) ~= 3,
            mr.user.bDir = mr.user.bDir';
        end
        
        mrfod = computeCSDFC(mr.dataAy ,mr.user.bDir,'mask',wm.dataAy>thres,'bDir',numdirs, ...
                             'numround',numround,'numit',numit,'lambda',lambda,'iso_sym',iso_sym,'iso_asym',iso_asym,'gamma_asym',gamma_asym,'gamma_sym',gamma_sym,'lambda_lb',lambda_lb, ...
                              'alpha',alpha,'kernwid',Daxial*mr.user.bfactor(1)/1000,'purefc_sym',purefc_sym,'purefc_asym',purefc_asym,'output_sym',symout,'mrprop',mr);
               
                          
                          
                          
        [p n e] = fileparts(hrname);
        if strcmpi('hardi',n(end-4:end)),
            n = n(1:end-5);
        end;
                
        outname = buildoutname(outname,p,n,'FOD','.mat');
            
        mrstruct_write(mrfod,outname);
        
        out.FOD = {outname};
        
   case 'doCSA'
      
        hrname = P.filename1{1};
        outname = P.newfilename;
        shcut = P.shcutoff;
        penal = P.laplacebeltrami;
        numdirs = P.numdirsout;

        mr = mrstruct_read(hrname);
        [mrfod] =  computeCSAdODF(mr,numdirs,shcut,penal);

        [p n e] = fileparts(hrname);
        if strcmpi('hardi',n(end-4:end)),
            n = n(1:end-5);
        end;
                
        outname = buildoutname(outname,p,n,'FOD','.mat');
            
        mrstruct_write(mrfod,outname);
        
        out.FOD = {outname};
            
    case 'doL1CSD_TFD'
      
        hrname = P.filename1{1};
        wmname = P.roidef.roiname{1};
        outname = P.newfilename;
       
       if isfield(P.fiout,'newfilename'), % mrstruct
           savemr_or_nifti = 1;
           outnameTFD = P.fiout.newfilename;
       elseif isfield(P.fiout,'newfinameniftireslice'), % nifti + reslice
           savemr_or_nifti = 0;
           ref = P.fiout.newfinameniftireslice.filename3{1};
           outnameTFD =  P.fiout.newfinameniftireslice.newfilename;
       elseif isfield(P.fiout,'newfinamenifti'),  % nifti
           savemr_or_nifti = 0;
           ref = [];
           outnameTFD =  P.fiout.newfinamenifti.newfilename;           
       end;
               
        
        
        params.osamp = P.osamp;

        params.csdL1.maxit = P.maxit;
        params.csdL1.lambda = P.lambdaL1;
        params.csdL1.powsm = P.powsm;
        params.csdL1.numdir = P.numdir;
        params.csdL1.D = [P.Dax P.Drad];

        params.TFD.wmthres = P.roidef.thresh;
        params.TFD.csfweight = 0;
        params.TFD.gmweight = 0;
        params.TFD.globalreg = 1;
        params.TFD.lambda = [1 P.lambda1 P.lambda2];
        params.TFD.planorder = P.planorder;        

        [mrres mrfod] = computeTensorFOD(hrname,wmname,[],[],params);

        [p n e] = fileparts(hrname);
        if strcmpi('hardi',n(end-4:end)),
            n = n(1:end-5);
        end;
                
        outname = buildoutname(outname,p,n,'FOD','.mat');                    
        mrstruct_write(mrfod,outname);
        
        out.FOD = {outname};
        
        if savemr_or_nifti == 1,
            outnameTFD = buildoutname(outnameTFD,p,n,'TFD','.mat');
            mrstruct_write(mrres,outnameTFD);                
        else,
           outnameTFD = buildoutname(outnameTFD,p,n,'TFD','.nii');
           if not(isempty(ref)),
               tmpname = '/tmp/niftitomrstructconversiontmp.nii';           
               tmpnamer = '/tmp/rniftitomrstructconversiontmp.nii';           
               [res, errStr]= mrstruct_to_nifti(mrres, tmpname , 'float32');
               flags.mean = false;
               flags.which = 1;
               spm_reslice({ ref tmpname },flags);
               system(['cp ' tmpnamer ' ' outnameTFD]);
           else,
               [res, errStr]= mrstruct_to_nifti(mrres, outnameTFD , 'float32');
           end;
        end;
        
        out.TFD = {outnameTFD};
        
    case 'doTFD'
      
        hrname = P.filename1{1};
        wmname = P.roidef.roiname{1};
       
       if isfield(P.fiout,'newfilename'), % mrstruct
           savemr_or_nifti = 1;
           outnameTFD = P.fiout.newfilename;
       elseif isfield(P.fiout,'newfinameniftireslice'), % nifti + reslice
           savemr_or_nifti = 0;
           ref = P.fiout.newfinameniftireslice.filename3{1};
           outnameTFD =  P.fiout.newfinameniftireslice.newfilename;
       elseif isfield(P.fiout,'newfinamenifti'),  % nifti
           savemr_or_nifti = 0;
           ref = [];
           outnameTFD =  P.fiout.newfinamenifti.newfilename;           
       end;
               
        
        
        params.osamp = P.osamp;

        params.TFD.wmthres = P.roidef.thresh;
        params.TFD.csfweight = 0;
        params.TFD.gmweight = 0;
        params.TFD.globalreg = 1;
        params.TFD.lambda = [1 P.lambda1 P.lambda2];
        params.TFD.planorder = P.planorder;        
        params.TFD.modifyFOD = P.modifyFOD;
        
        [mrres] = computeTensorFD(hrname,wmname,[],[],params);

        [p n e] = fileparts(hrname);
        if strcmpi('fod',n(end-2:end)),
            n = n(1:end-3);
        end;
                        
        if savemr_or_nifti == 1,
            outnameTFD = buildoutname(outnameTFD,p,n,'TFD','.mat');
            mrstruct_write(mrres,outnameTFD);                
        else,
           outnameTFD = buildoutname(outnameTFD,p,n,'TFD','.nii');
           if not(isempty(ref)),
               tmpname = '/tmp/niftitomrstructconversiontmp.nii';           
               tmpnamer = '/tmp/rniftitomrstructconversiontmp.nii';           
               [res, errStr]= mrstruct_to_nifti(mrres, tmpname , 'float32');
               flags.mean = false;
               flags.which = 1;
               spm_reslice({ ref tmpname },flags);
               system(['cp ' tmpnamer ' ' outnameTFD]);
           else,
               [res, errStr]= mrstruct_to_nifti(mrres, outnameTFD , 'float32');
           end;
        end;
        
       
        
        out.TFD = {outnameTFD};
                
    case 'deformFOD'    
        
       fodname =  P.filename2{1};
       iy = P.filedefi{1};
       y = P.filedef{1};
       interp = P.interp;
       modulate = P.modulate;
       outname = P.newfilename;        
       noutdir = P.noutdir;
       
       deformedFOD = mrFOD_warp(iy,y,fodname,'modulate',modulate,'method',interp,'noutdir',noutdir);
        
       [p n e] = fileparts(fodname);
       if strcmpi('fod',n(end-2:end)),
           n = n(1:end-3);
       end;
       
       outname = buildoutname(outname,p,n,'defFOD','.mat');
              
       mrstruct_write(deformedFOD,outname);           
       
       out.FOD = {outname};

    case 'deformFDalong'    
        
       fodname =  P.filename2{1};
       fdname =  P.filename4{1};
       iy = P.filedefi{1};
       y = P.filedef{1};
       interp = P.interp;
       modulate = P.modulate;
       outname = P.newfilename;        

       fod = mrstruct_read(fodname);
       mfod = mean(fod.dataAy,4);

       fd = mrstruct_read(fdname);
       
       ratio = fd.dataAy ./ (eps+mfod);
       for k = 1:size(fod.dataAy,4),
           fod.dataAy(:,:,:,k) = fod.dataAy(:,:,:,k) .* ratio;
       end
              
       fod = mrFOD_warp(iy,y,fod,'modulate',modulate==1,'method',interp);
       
       fd.dataAy = mean(fod.dataAy,4);
       
       [p n e] = fileparts(fodname);
       
       outname = buildoutname(outname,p,n,'def','.mat');
       
       mrstruct_write(fd,outname);           

       
    case 'smoothFOD'    
       fodname =  P.filename2{1};
       axsig = P.axialsigma;
       rasig = P.radialsigma;
       outname = P.newfilename;        
      
       if rasig > axsig,
           warning('smaller axial width than radial not yet implemented');
           return;
       end;
       
       mr = mrstruct_read(fodname);
       mrsm = mr;
       bDir = mr.user.bDir;
       sz = size(mr.dataAy);
       [X Y Z] = ndgrid(-ceil(sz(1)/2):floor(sz(1)/2)-1,-ceil(sz(2)/2):floor(sz(2)/2)-1,-ceil(sz(3)/2):floor(sz(3)/2)-1);
       X = fftshift(X); Y = fftshift(Y); Z = fftshift(Z);
       R2 = X.^2 + Y.^2 + Z.^2;
       isog = (exp(-R2/(2*axsig^2)));       
       mrsm.dataAy(isnan(mrsm.dataAy(:))) = 0;
       for k = 1:size(mr.dataAy,4)
           fgauss = fftn(exp(- (R2-(X*bDir(1,k) +  Y*bDir(2,k) + Z*bDir(3,k)).^2) * (1/(2*rasig.^2) - 1/(2*axsig.^2))) .* isog);
           mrsm.dataAy(:,:,:,k) = real(ifftn(fftn(squeeze(mrsm.dataAy(:,:,:,k))).*fgauss));
           fprintf('.');
       end;
                                
       [p n e] = fileparts(fodname);

       
       ex = '';
       if strcmpi('fod',n(end-2:end)),
           ex = '_FOD';
           n = n(1:end-3);
       end;
       if strcmpi('hardi',n(end-4:end)),
           ex = '_HARDI';
           n = n(1:end-5);
       end;  
       
       
       
       outname = buildoutname(outname,p,n,['smoothed' ex],'.mat');
       mrsm.dataAy = single(mrsm.dataAy);
       
       mrstruct_write(mrsm,outname);    
              
       out.FOD = {outname};
       
   
    case 'smoothFODRic'    
       fodname =  P.filename2{1};
       axsig = P.axialsigma;
       rasig = P.radialsigma;
       ksz = P.ksz;
       noiselevel = P.noiselevel;
       outname = P.newfilename;        
      
       if rasig > axsig,
           warning('smaller axial width than radial not yet implemented');
           return;
       end;
       
       mr = mrstruct_read(fodname);
       
       mrsm = smoothAnisoRice(mr,axsig,rasig,ksz,noiselevel);
                                               
       [p n e] = fileparts(fodname);

       ex = '';
       if strcmpi('fod',n(end-2:end)),
           ex = '_FOD';
           n = n(1:end-3);
       end;
       if strcmpi('hardi',n(end-4:end)),
           ex = '_HARDI';
           n = n(1:end-5);
       end;
       
       outname = buildoutname(outname,p,n,['smoothed' ex],'.mat');
       mrsm.dataAy = single(mrsm.dataAy);
       
       mrstruct_write(mrsm,outname);    
              
       out.FOD = {outname};     
       
       
case 'resampleHARDI'    
     hardiname =  P.filename1{1};
     hardiref =  P.filenameHARDI1{1};
    
     mr = mrstruct_read(hardiname);
     ref = mrstruct_read(hardiref);
     bDir = ref.user.bDir;
     clear ref;
     
     mr = resampleHARDI(mr,bDir);
     
     [p n e] = fileparts(hardiname);
     outname = fullfile(p,['r' n e]);
     
     mrstruct_write(mr,outname);
     
     out.HARDI = {outname};
    
    
case 'getsiFOD'    
       fodname =  P.filename2{1};
       if isfield(P.fiout,'newfilename'), % mrstruct
           savemr_or_nifti = 1;
           outname = P.fiout.newfilename;
       elseif isfield(P.fiout,'newfinameniftireslice'), % nifti + reslice
           savemr_or_nifti = 0;
           ref = P.fiout.newfinameniftireslice.filename3{1};
           outname =  P.fiout.newfinameniftireslice.newfilename;
       elseif isfield(P.fiout,'newfinamenifti'),  % nifti
           savemr_or_nifti = 0;
           ref = [];
           outname =  P.fiout.newfinamenifti.newfilename;           
       end;
       
              
       mr = mrstruct_read(fodname);
       mrmfod = mr;
       if P.sitype == 0,  % GFA
            mrmfod.dataAy = sqrt(squeeze(var(mr.dataAy,[],4)) ./ squeeze(sum(mr.dataAy.^2,4)) * size(mr.dataAy,4));
            mrfod.user.type = 'GFA';
       elseif P.sitype == 1,  % mean FOD
            mrmfod.dataAy = squeeze(mean(mr.dataAy,4));
            mrfod.user.type = 'meanFOD';
       elseif P.sitype == 2,  % mean FOD b0norm
            mrmfod.dataAy = squeeze(mean(mr.dataAy,4))./(eps+mr.user.b0avg);
            mrfod.user.type = 'meanFODb0norm';
       elseif P.sitype == 3, %max FOD
            mrmfod.dataAy = max(mr.dataAy,[],4);
            mrfod.user.type = 'maxFOD';
       elseif P.sitype == 4, %max FOD b0norm
            mrmfod.dataAy = max(mr.dataAy,[],4)./(eps+mr.user.b0avg);
            mrfod.user.type = 'maxFODb0norm';
       end;
            
       mrmfod.dim4 = 'unused';
              
       [p n e] = fileparts(fodname);

       if strcmpi('fod',n(end-2:end)),
           n = n(1:end-3);
       end;       
       
       if savemr_or_nifti == 1,
           outname = buildoutname(outname,p,n,mrfod.user.type,'.mat');
           mrstruct_write(mrmfod,outname);                
       else,
           outname = buildoutname(outname,p,n,mrfod.user.type,'.nii');
           if not(isempty(ref)),
               tmpname = '/tmp/niftitomrstructconversiontmp.nii';           
               tmpnamer = '/tmp/rniftitomrstructconversiontmp.nii';           
               [res, errStr]= mrstruct_to_nifti(mrmfod, tmpname , 'float32');
               flags.mean = false;
               flags.which = 1;
               spm_reslice({ ref tmpname },flags);
               system(['cp ' tmpnamer ' ' outname]);
           else,
               [res, errStr]= mrstruct_to_nifti(mrmfod, outname , 'float32');
           end;
       end;
           
       out.SIndex = {outname};    
           
end



function outname = buildoutname(outname,p,n,str,ext)
    if length(outname) >= 4,
        if strcmpi(outname(1:4),ext),
            if length(outname) > 5 & outname(5) == ',',
                suffix = outname(6:end);
                outname = fullfile(p,[n suffix str ext]);
                return;
            else,
                outname = fullfile(p,[n str ext]);
                return;
            end;
        end
    end;
    outname = fullfile(p,outname);


function m = loadsegmentation(name,osamp,mrprop,resliceflags)
mrprop.dataAy = zeros(size(mrprop.dataAy,1),size(mrprop.dataAy,2),size(mrprop.dataAy,3)); mrprop.dim4 = 'unused';

[p n e] = fileparts(name);

if isempty(mrprop.edges),
    mrprop.edges = eye(4);
end;

mrprop.edges(1:3,1:3) = mrprop.edges(1:3,1:3)/osamp;
mrprop.vox = mrprop.vox/osamp;
mrprop.dataAy = zeros(size(mrprop.dataAy)*osamp);

if strcmpi(e(1:4),'.mat'),
    if osamp ~= 1,
        warning('oversampling > 1 not yet supported for masks of type mat');
    end;
    if length(e) > 4,
        name = fullfile(p,[n e(1:4)]);
    end;
    ms = load(name);
    if isfield(ms,'mrStruct');
        m = ms.mrStruct;
    elseif isfield(ms,'maskCell');    
        num = 1;
        if length(e) > 4,
            num = str2num(e(6:end));
        end;
        m = ms.mrsProp;
        m.dataAy = ms.maskCell{num};
    end;
    
elseif strcmpi(e(1:4),'.nii'),
    if exist('resliceflags'),
        m = nifti_to_mrstruct('volume',{name},mrprop,pwd,resliceflags);    
    else
        m = nifti_to_mrstruct('volume',{name},mrprop);    
    end;
end;


