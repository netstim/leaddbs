function ea_gensurfice_temps
disp('Generating surfaces of template space (left/right hemisphere)');
setenv('FSLOUTPUTTYPE', 'NIFTI');
cmd = [ea_getExec([ea_getearoot,'ext_libs',filesep,'fsl',filesep,'bet2'], escapePath = 1) ' ' ea_space 't1.nii ' ea_space 't1b.nii'];
ea_dispt('Skullstripping T1');

ea_runcmd(cmd);
ea_dispt('Splitting T1 into left & right hemispheres');

ea_split_nii_lr([ea_space,'t1b.nii']);
delete([ea_space,'t1b.nii']);

vizz=0;
sidest={'r','l'};
for side=1:2
    % LH
    disp(['Processing ',sidest{side},' hemisphere']);
    ea_dispt('Generating surface');
    left=ea_load_nii([ea_space,'t1b_',sidest{side},'.nii']);
    left.img=left.img>51;
    CC=bwconncomp(left.img);
    sizes=cellfun(@length,CC.PixelIdxList,'UniformOutput',false);
    sizes=cell2mat(sizes');
    [~,mix]=max(sizes);
    left.img(:)=0;
    left.img(CC.PixelIdxList{mix})=1;
    left.img=imfill(left.img,'holes');

    clear fv
    tfv=isosurface(left.img,0.5);
    fv(1).vertices=left.mat*[tfv.vertices(:,2),tfv.vertices(:,1),tfv.vertices(:,3),ones(length(tfv.vertices),1)]';
    fv(1).vertices=fv(1).vertices(1:3,:)';
    fv(1).faces=tfv.faces;
    tfv=isocaps(left.img,0.5);
    fv(2).vertices=left.mat*[tfv.vertices(:,2),tfv.vertices(:,1),tfv.vertices(:,3),ones(length(tfv.vertices),1)]';
    fv(2).vertices=fv(2).vertices(1:3,:)';
    fv(2).faces=tfv.faces;

    fv=ea_concatfv(fv);
    ea_dispt('Smoothing surface');

    fvs=ea_smoothpatch(fv,1,10,1);
    fvs.facevertexcdata=repmat(2,length(fvs.vertices),1);
    if vizz
        figure
        h=patch('vertices',fvs.vertices,'faces',fvs.faces,'FaceColor','interp','EdgeColor','none','FaceVertexCdata',fvs.facevertexcdata);
        a=camlight('headlight');
        axis equal
        % fvs.facevertexnormals=h.FaceNormals;
    end
    ea_dispt('Exporting surface');


    fvs=ea_mapcolvert2face(fvs);
    ea_stlwrite([ea_space,'surf.',sidest{side},'h.stl'],fvs,'FACECOLOR',fvs.facevertexcdata);
end
ea_dispt('Cleaning up');

delete([ea_space,'t1b_l.nii']);
delete([ea_space,'t1b_r.nii']);
ea_dispt('');
