function ea_importfssegmentations(options)

if options.prefs.fs.subcorticalseg.thalamus
    % thalamus
    thaldict=getthalmap;

    temp=ea_getleadtempdir;

    system([options.prefs.fs.dir,filesep,'bin',filesep,...
        'mri_convert',...
        ' ',ea_path_helper(fullfile(options.subj.freesurferDir,['sub-',options.subj.subjId],'mri','ThalamicNuclei.mgz')),...
        ' ',ea_path_helper(fullfile(temp,'thalamic.nii'))]);

    thalnii=ea_load_nii(fullfile(temp,'thalamic.nii'));

    outdir=fullfile(options.subj.atlasDir,'FreeSurfer_Segmentations');
    ea_mkdir(fullfile(outdir,'lh'));
    ea_mkdir(fullfile(outdir,'rh'));
    ea_mkdir(fullfile(outdir,'midline'));

    for thalnuc=1:length(thaldict.labels)
        thisnuc=thalnii;

        % right
        thisnuc.img=thalnii.img==thaldict.idx.right(thalnuc);
        if any(thisnuc.img(:))
            thisnuc.fname=fullfile(outdir,'rh',[thaldict.labels{thalnuc},'.nii']);
            ea_write_nii(thisnuc);
            gzip(thisnuc.fname);
            ea_delete(thisnuc.fname);
        end
        % left
        thisnuc.img=thalnii.img==thaldict.idx.left(thalnuc);
        if any(thisnuc.img(:))

            thisnuc.fname=fullfile(outdir,'lh',[thaldict.labels{thalnuc},'.nii']);
            ea_write_nii(thisnuc);
            gzip(thisnuc.fname);
            ea_delete(thisnuc.fname);
        end
    end
end
if options.prefs.fs.subcorticalseg.hippo_amygdala
    % hippocampal subfields & amygdala
    hippdict=gethipmap;

    temp=ea_getleadtempdir;

    system([options.prefs.fs.dir,filesep,'bin',filesep,...
        'mri_convert',...
        ' ',ea_path_helper(fullfile(options.subj.freesurferDir,['sub-',options.subj.subjId],'mri','rh.hippoAmygLabels.mgz')),...
        ' ',ea_path_helper(fullfile(temp,'rh_hippo.nii'))]);

    rhippnii=ea_load_nii(fullfile(temp,'rh_hippo.nii'));

    system([options.prefs.fs.dir,filesep,'bin',filesep,...
        'mri_convert',...
        ' ',ea_path_helper(fullfile(options.subj.freesurferDir,['sub-',options.subj.subjId],'mri','lh.hippoAmygLabels.mgz')),...
        ' ',ea_path_helper(fullfile(temp,'lh_hippo.nii'))]);

    lhippnii=ea_load_nii(fullfile(temp,'lh_hippo.nii'));


    outdir=fullfile(options.subj.atlasDir,'FreeSurfer_Segmentations');
    ea_mkdir(fullfile(outdir,'lh'));
    ea_mkdir(fullfile(outdir,'rh'));

    for hippnuc=1:length(hippdict.labels)
        thisnuc=rhippnii;

        % right
        thisnuc.img=rhippnii.img==hippdict.idx(hippnuc);
        if any(thisnuc.img(:))
            thisnuc.fname=fullfile(outdir,'rh',[hippdict.labels{hippnuc},'.nii']);
            ea_write_nii(thisnuc);
            gzip(thisnuc.fname);
            ea_delete(thisnuc.fname);
        end
        thisnuc=lhippnii;

        % left
        thisnuc.img=lhippnii.img==hippdict.idx(hippnuc);
        if any(thisnuc.img(:))

            thisnuc.fname=fullfile(outdir,'lh',[hippdict.labels{hippnuc},'.nii']);
            ea_write_nii(thisnuc);
            gzip(thisnuc.fname);
            ea_delete(thisnuc.fname);
        end
    end

end
if options.prefs.fs.subcorticalseg.brainstem
    % brainstem
    bsdict=getbsmap;

    temp=ea_getleadtempdir;

    system([options.prefs.fs.dir,filesep,'bin',filesep,...
        'mri_convert',...
        ' ',ea_path_helper(fullfile(options.subj.freesurferDir,['sub-',options.subj.subjId],'mri','brainstemSsLabels.mgz')),...
        ' ',ea_path_helper(fullfile(temp,'brainstem.nii'))]);

    bsnii=ea_load_nii(fullfile(temp,'brainstem.nii'));

    outdir=fullfile(options.subj.atlasDir,'FreeSurfer_Segmentations');
    ea_mkdir(fullfile(outdir,'lh'));
    ea_mkdir(fullfile(outdir,'rh'));

    for bs=1:length(bsdict.labels)
        thisnuc=bsnii;

        % midline
        thisnuc.img=bsnii.img==bsdict.idx(bs);
        if any(thisnuc.img(:))

            thisnuc.fname=fullfile(outdir,'midline',[bsdict.labels{bs},'.nii']);
            ea_write_nii(thisnuc);
            gzip(thisnuc.fname);
            ea_delete(thisnuc.fname);
        end
    end
end

function bsdict=getbsmap
bsdict.labels={'Midbrain','Pons','Medulla','SCP'};
bsdict.idx=[173,174,175,178];

function hippdict=gethipmap
hippdict.labels={'lateral nucleus','basolateral nucleus','basal nucleus','centromedial nucleus','central nucleus','medial nucleus','cortical nucleus','accessory basal nucleus',...
    'corticoamygdaloid transition zone','anterior amygdaloid area','fusion amygdala HP FAH','hippocampal amygdala transition','endopiriform nucleus','lateral nucleus olfactory tract',...
    'paralaminar nucleus','intercalated nucleus','prepiriform cortex','periamygdaloid cortex','envelope amygdala','extranuclear amygdala'};
hippdict.idx=[7001,7002,7003,7004,7005,7006,7007,7008,7009,7010,7011,7012,7013,7014,7015,7016,7017,7018,7019,7020];

function thaldict=getthalmap
thaldict.labels={'AV','CeM','CL','CM','LD','LGN','LP','L-Sg','MDl','MDm','MGN','MV(Re)','Pc','Pf','Pt','PuA','PuI','PuL','PuM','R','VA','VAmc','VLa','VLp','VM','VPL','PaV'};

thaldict.idx.left=[8103,8104,8105,8106,8108,8109,8110,8111,8112,8113,8115,8116,8117,8118,8119,8120,8121,8122,8123,8125,8126,8127,8128,8129,8130,8133,8134];
thaldict.idx.right=[8203,8204,8205,8206,8208,8209,8210,8211,8212,8213,8215,8216,8217,8218,8219,8220,8221,8222,8223,8225,8226,8227,8228,8229,8230,8233,8234];
