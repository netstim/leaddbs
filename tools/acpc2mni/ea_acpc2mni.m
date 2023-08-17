function fid=ea_acpc2mni(varargin)
% leadfig can be either a handle to the lead figure with selected patient
% directories or a cell of patient directories.
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


fidpoints_mm=[ 0.25,   1.298,   -5.003    % AC
              -0.188, -24.756,  -2.376    % PC
               0.25,   1.298,    55];     % Midsag

cfg=varargin{1};
leadfig=varargin{2};

if iscell(leadfig)
    uidir=leadfig;
else
    uidir=getappdata(leadfig,'uipatdir');
end

% prompt for ACPC-coordinates:
acpc=[cfg.xmm,cfg.ymm,cfg.zmm];
if cfg.mapmethod
    if nargin<5
        [FileName,PathName] = uiputfile('ACPC2MNI_Mapping.nii','Save Mapping...');
    else
        [PathName,FileName,Ext] = fileparts(varargin{5});
        FileName=[FileName,Ext];
    end
end

if isempty(uidir)
    ea_error('Please choose and normalize patients first.');
end

%disp('*** Converting ACPC-coordinates to MNI based on normalizations in selected patients.');
%ea_dispercent(0,'Iterating through patients');
for pt=1:length(uidir)
 %   ea_dispercent(pt/length(uidir));
    if ~strcmp(uidir{pt}(end), filesep)
        directory = [uidir{pt},filesep];
    else
        directory = uidir{pt};
    end
    
    options = ea_getptopts(directory);
    coreg_files = cellfun(@(x) options.subj.coreg.anat.preop.(x), fieldnames(options.subj.coreg.anat.preop), 'uni', 0);

    if nargin>5 % determine whether to use manually or automatically defined AC/PC
        automan=varargin{6};
    else
        if exist(options.subj.acpc.acpcManual,'file') % manual AC/PC definition present
            automan='manual';
        else
            automan='auto';
        end
    end
    
    if nargin>2
        whichnormmethod=varargin{3};
        template=varargin{4};
    else
        [whichnormmethod,template]=ea_whichnormmethod(directory);
    end

    fidpoints_vox=ea_getfidpoints(fidpoints_mm,template);

    switch automan
        case 'auto' % auto AC/PC detection
            if ~exist(options.subj.acpc.acpcAutodetect,'file')
                % warp into patient space
                [fpinsub_mm] = ea_map_coords(fidpoints_vox', ...
                                            template,...
                                            [options.subj.subjDir,filesep,'forwardTransform'],...
                                            coreg_files{1},...
                                            whichnormmethod);

                fpinsub_mm=fpinsub_mm';

                try
                    fid(pt).AC=fpinsub_mm(1,:);
                catch
                    keyboard
                end

                fid(pt).PC=fpinsub_mm(2,:);
                fid(pt).MSP=fpinsub_mm(3,:);
                acpc_fiducials=fid(pt); % save for later use
                if ~isfolder(options.subj.acpcDir)
                    mkdir(options.subj.acpcDir);
                end
                save(options.subj.acpc.acpcAutodetect,'-struct','acpc_fiducials'); clear acpc_fiducials
            else
                tmp=load(options.subj.acpc.acpcAutodetect);
                fid(pt).AC=tmp.AC;
                fid(pt).PC=tmp.PC;
                fid(pt).MSP=tmp.MSP;
            end

        case {'manual'} % manual AC/PC definition, assume F.fcsv file inside pt folder
            copyfile(options.subj.acpc.acpcManual,strrep(options.subj.acpc.acpcManual,'.fcsv','.dat'))
            Ct=readtable(strrep(options.subj.acpc.acpcManual,'.fcsv','.dat'));

            % AC
            cnt=1;
            fid(pt).AC=zeros(1,3);
            for i=17:19
                thisval=table2array(Ct(i,1));
                fid(pt).AC(cnt)=str2double(thisval{1});
                cnt=cnt+1;
            end

            % PC
            cnt=1;
            fid(pt).PC=zeros(1,3);
            for i=31:33
                thisval=table2array(Ct(i,1));
                fid(pt).PC(cnt)=str2double(thisval{1});
                cnt=cnt+1;
            end

            % MSP
            cnt=1;
            fid(pt).MSP=zeros(1,3);
            for i=45:47
                thisval=table2array(Ct(i,1));
                fid(pt).MSP(cnt)=str2double(thisval{1});
                cnt=cnt+1;
            end

        case 'mnidirect'
            fid(pt).AC=[0.25,   1.298,   -5.003];
            fid(pt).PC=[-0.188, -24.756,  -2.376];
            fid(pt).MSP=[0.25,   1.298,    55];

    end

    % x-dimension
    A=fid(pt).MSP-fid(pt).AC;
    B=fid(pt).PC-fid(pt).AC;
    xvec=cross(A,B); %normal to given plane
    xvec=xvec/norm(xvec);
    % y-dimension (just move from ac to pc and scale by y dimension):
    yvec=(fid(pt).PC-fid(pt).AC);
    yvec=yvec/norm(yvec);
    yvec=-yvec; % this vector points to anterior of the AC!! (but acpc y coordinates are given with - sign).

    zvec=cross(xvec,yvec);
    zvec=zvec/norm(zvec);

    switch cfg.acmcpc
        case 1 % relative to AC:
            warpcoord_mm=fid(pt).AC+acpc(1)*xvec+acpc(2)*yvec+acpc(3)*zvec;
        case 2 % relative to midcommissural point:
            warpcoord_mm=mean([fid(pt).AC;fid(pt).PC],1)+acpc(1)*xvec+acpc(2)*yvec+acpc(3)*zvec;
        case 3 % relative to PC:
            warpcoord_mm=fid(pt).PC+acpc(1)*xvec+acpc(2)*yvec+acpc(3)*zvec;
    end

    anat=ea_load_nii(coreg_files{1});
    warpcoord_mm=[warpcoord_mm';1];
    warpcoord_vox=anat.mat\warpcoord_mm;
    warpcoord_vox=warpcoord_vox(1:3);
    fid(pt).WarpedPointNative=warpcoord_mm(1:3)';

    % re-warp into MNI
    if ~isfield(cfg,'native') || ~cfg.native % when working in native space, no need to warp acpc back to mni at all.
        switch automan
            case 'mnidirect'
                fid(pt).WarpedPointMNI=warpcoord_mm(1:3)';
            otherwise
                [warpinmni_mm] = ea_map_coords(warpcoord_vox,...
                                 coreg_files{1},...
                                 [options.subj.subjDir,filesep,'inverseTransform'],...
                                 template,...
                                 whichnormmethod);
                try
                    warppts(pt,:)=warpinmni_mm';
                catch
                    keyboard
                end
                fid(pt).WarpedPointMNI=warppts(pt,:);
        end
    end

    if cfg.mapmethod==2
        anat.img(:)=0;
        anat.img(round(warpcoord_vox(1)),round(warpcoord_vox(2)),round(warpcoord_vox(3)))=1;
        anat.fname=[options.subj.acpcDir,'ACPCquerypoint.nii'];
        spm_write_vol(anat,anat.img);

        % warp into nativespace
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.subj.norm.transform.forwardBaseName, 'spm.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[options.subj.acpcDir,'ACPCquerypoint.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {options.subj.acpcDir};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        spm_jobman('run',{matlabbatch});
        clear matlabbatch
        wfis{pt}=[options.subj.acpcDir,'wACPCquerypoint.nii'];
    end
end

% create clear cut version:
if cfg.mapmethod==1
    bb=ea_load_nii([ea_space,'bb.nii']);
    bb.img(:)=0;
    warppts_vox=[warppts';ones(1,size(warppts,1))];
    warppts_vox=round(bb.mat\warppts_vox);

    for pnt=1:size(warppts_vox,2)
        try
            bb.img(warppts_vox(1,pnt),warppts_vox(2,pnt),warppts_vox(3,pnt))=1;
        end
    end

    bb.fname=[PathName,FileName];
    spm_write_vol(bb,bb.img);
elseif cfg.mapmethod==2
    % create innativespacemapped files:
    matlabbatch{1}.spm.util.imcalc.input = wfis;
    matlabbatch{1}.spm.util.imcalc.output = [FileName];
    matlabbatch{1}.spm.util.imcalc.outdir = {PathName};
    matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch
end

if cfg.mapmethod
    % smooth clear version:
    matlabbatch{1}.spm.spatial.smooth.data = {[PathName,FileName,',1']};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [1 1 1];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    [pth,fn,ext]=fileparts(bb.fname);
    ea_crop_nii(fullfile(pth,[fn,ext]));
    ea_crop_nii(fullfile(pth,['s',fn,ext]));
end

assignin('base','fid',fid);


function fidpoints_vox=ea_getfidpoints(fidpoints_mm,tempfile)

V=spm_vol(tempfile);
fidpoints_vox=V(1).mat\[fidpoints_mm,ones(size(fidpoints_mm,1),1)]';
fidpoints_vox=fidpoints_vox(1:3,:)';
