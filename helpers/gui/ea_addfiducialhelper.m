function ea_addfiducialhelper(~,~,handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
cnt=1;
        spacedef=ea_getspacedef;

while 1
    try close(TPpoint); end
    uuid{cnt}=ea_generate_uuid;
    
    % choose in template:
    spm('defaults', 'fmri');
    Fgraph = spm_figure('GetWin', 'Graphics');
    movegui(Fgraph,'northwest');
    Finter = spm('FnUIsetup','Select Fiducials', 0);
    movegui(Finter,'south');
    
    figure(Fgraph); clf;
    
    spm_orthviews('Reset');
    spm_orthviews('Image', [ea_space,spacedef.templates{1},'.nii']);
    colormap('gray');
    cameratoolbar('resetcamera')
    cameratoolbar('close')
    rotate3d off;
    
    if spm_input(['# ',num2str(cnt),': Select point on template and click'] , 1,'OK|Stop', [1,0], 1)
        tppos = spm_orthviews('Pos')';
        FR=getframe(Fgraph);
        TPpoint=figure('Name','Fiducial on template','NumberTitle','off','Position',Fgraph.Position,'Visible','off','Toolbar','none','Menubar','none');
        movegui(TPpoint,'northeast');
        imshow(FR.cdata);
        TPpoint.Visible='on';
    else
        break
    end
    
    
    % choose in patient(s):
    
    for pt=1:length(uipatdir)
        
        options=ea_getptopts([uipatdir{pt},filesep]);
        [~,presentfiles]=ea_assignpretra(options);
        ptspace{pt}=[uipatdir{pt},filesep,presentfiles{1}];
        
        figure(Fgraph); clf;
        spm_orthviews('Reset');
        spm_orthviews('Image', ptspace{pt});
        colormap('gray');
        cameratoolbar('resetcamera')
        cameratoolbar('close')
        rotate3d off;
        
        if spm_input(['# ',num2str(cnt),': Select corresponding point on patient anchor modality and click'] , 1,'OK|Stop', [1,0], 1)
            patpos(pt,:) = spm_orthviews('Pos')';
        else
            exit=1; % to be able to break out of while loop as well.
            break
        end
    end
    if exist('exit','var')
        break
    end
    try close(Fgraph); end
    try close(Finter); end
    try close(TPpoint); end
    
    fid(cnt).template=tppos;
    fid(cnt).patient=patpos;
    cnt=cnt+1;
end
try close(Fgraph); end
try close(Finter); end
try close(TPpoint); end
if ~exist('fid','var')
    disp('No fiducials defined');
    return
else
    disp('Fiducials defined. Logging & smoothing them to be used in next ANTs-based transform.');
end
for f=1:length(fid)
    for pt=1:length(uipatdir)        
        % define in template:
        ea_mkdir([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace]);
        ea_spherical_roi([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,uuid{f},'.nii'],fid(f).template,ea_species_adjustsize(10),0,[ea_space,spacedef.templates{1},'.nii']);
        tfis{pt}{f}=[uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,uuid{f},'.nii'];
        %smoothgzip([ea_space,'fiducials'],[uuid{f},'.nii']);
        
        % define in pt:
        ea_mkdir([uipatdir{pt},filesep,'fiducials',filesep,'native']);
        ea_spherical_roi([uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,uuid{f},'.nii'],fid(f).patient(pt,:),ea_species_adjustsize(10),0,ptspace{pt});
        pfis{pt}{f}=[uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,uuid{f},'.nii'];
        %smoothgzip([uipatdir{pt},filesep,'fiducials'],[uuid{f},'.nii']);
    end
end

% flatten ROI:
for pt=1:length(uipatdir)
    if length(tfis{pt})>1
        fguid=ea_generate_uuid;
        matlabbatch{1}.spm.util.imcalc.input = tfis{pt}';
        matlabbatch{1}.spm.util.imcalc.output = [fguid,'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[uipatdir{pt},filesep,'fiducials',filesep,ea_getspace]};
        matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 512;
        spm_jobman('run',{matlabbatch});
        smoothgzip([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace],[fguid,'.nii']);
    else
        [pathn,filenn]=fileparts(tfis{pt}{1});
        smoothgzip(pathn,[filenn,'.nii']);
    end
        ea_delete(tfis{pt});

    if length(pfis{pt})>1
        matlabbatch{1}.spm.util.imcalc.input = pfis{pt}';
        matlabbatch{1}.spm.util.imcalc.output = [fguid,'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[uipatdir{pt},filesep,'fiducials',filesep,'native']};
        matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 512;
        spm_jobman('run',{matlabbatch});
        smoothgzip([uipatdir{pt},filesep,'fiducials',filesep,'native'],[fguid,'.nii']);
    else
        [pathn,filenn]=fileparts(pfis{pt}{1});
        smoothgzip(pathn,[filenn,'.nii']);
    end
    ea_delete(pfis{pt});
end


function smoothgzip(pathn,filen)

kernel=ea_species_adjustsize(8);

matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pathn,filen)};
matlabbatch{1}.spm.spatial.smooth.fwhm = [kernel kernel kernel];
matlabbatch{1}.spm.spatial.smooth.dtype = 512;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
movefile(fullfile(pathn,['s',filen]),fullfile(pathn,filen));
gzip(fullfile(pathn,filen));
delete(fullfile(pathn,filen));



