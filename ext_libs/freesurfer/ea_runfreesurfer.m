function ea_runfreesurfer(options)

hastb=ea_hastoolbox('freesurfer');

if ~hastb
    ea_error('FreeSurfer needs to be installed and connected to Lead-DBS');
end
hastb=ea_hastoolbox('fsl');
if ~hastb
    ea_error('FSL needs to be installed and connected to Lead-DBS');
end


presentfiles=fieldnames(options.subj.coreg.anat.preop);
ea_mkdir(options.subj.freesurferDir)

if options.prefs.fs.reconall.do
    % run recon-all:
    system([options.prefs.fs.dir,filesep,'bin',filesep,...
        'recon-all',...
        ' -subjid ',['sub-',options.subj.subjId],...
        ' -i ',[options.subj.coreg.anat.preop.(presentfiles{1})],...
        ' -all ',...
        ' -sd ',[options.subj.freesurferDir]]);
end

if options.prefs.fs.subcorticalseg.do
    % run subcortical segmentations: http://surfer.nmr.mgh.harvard.edu/fswiki/BrainstemSubstructures

    % segment_subregions
    setenv('SUBJECTS_DIR',[options.subj.freesurferDir]);

    if options.prefs.fs.subcorticalseg.thalamus
        system([options.prefs.fs.dir,filesep,'bin',filesep,...
            'segment_subregions',...
            ' --cross ',['sub-',options.subj.subjId],...
            ' thalamus']);
        [~,fsver]=fileparts(options.prefs.fs.dir(1:end-1));
        parsestr=['Thalamic nuclei were automatically segmented using FreeSurfer version ',fsver,' following the approach introduced by Iglesias et al. 2018.'];
        refs={'Iglesias JE, Insausti R, Lerma-Usabiaga G, Bocchetta M, Van Leemput K, Greve DN, van der Kouwe A, Fischl B, Caballero-Gaudes C, Paz-Alonso PM. A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology. NeuroImage. 2018;183:314-326. doi:10.1016/j.neuroimage.2018.08.012'};
        ea_methods(options,parsestr,refs)
    end
    if options.prefs.fs.subcorticalseg.hippo_amygdala
        system([options.prefs.fs.dir,filesep,'bin',filesep,...
            'segment_subregions',...
            ' --cross ',['sub-',options.subj.subjId],...
            ' hippo-amygdala']);
        [~,fsver]=fileparts(options.prefs.fs.dir(1:end-1));
        parsestr=['Hippocampal subfields (Iglesias et al., 2015) and nuclei of the amygdala (Saygin et al., 2017) were automatically segmented using freesurfer version ',fsver,'.'];
        refs={'Iglesias JE, Augustinack JC, Nguyen K, Player CM, Player A, Wright M, Roy N, Frosch MP, McKee AC, Wald LL, Fischl B, Van Leemput K. A computational atlas of the hippocampal formation using ex vivo , ultra-high resolution MRI: Application to adaptive segmentation of in vivo MRI. NeuroImage. 2015;115:117-137. doi:10.1016/j.neuroimage.2015.04.042',...
            'Saygin ZM, Kliemann D, Iglesias JE, van der Kouwe AJW, Boyd E, Reuter M, Stevens A, Van Leemput K, McKee A, Frosch MP, Fischl B, Augustinack JC. High-resolution magnetic resonance imaging reveals nuclei of the human amygdala: manual segmentation to automatic atlas. NeuroImage. 2017;155:370-382. doi:10.1016/j.neuroimage.2017.04.046'};
        ea_methods(options,parsestr,refs)
    end
    if options.prefs.fs.subcorticalseg.brainstem
        system([options.prefs.fs.dir,filesep,'bin',filesep,...
            'segment_subregions',...
            ' --cross ',['sub-',options.subj.subjId],...
            ' brainstem']);
        [~,fsver]=fileparts(options.prefs.fs.dir(1:end-1));
        parsestr=['Brainstem structures were automatically segmented using FreeSurfer version ',fsver,' following the approach introduced by Iglesias et al. 2015.'];
        refs={'Iglesias JE, Van Leemput K, Bhatt P, Casillas C, Dutt S, Schuff N, Truran-Sacrey D, Boxer A, Fischl B. Bayesian segmentation of brainstem structures in MRI. NeuroImage. 2015;113:184-195. doi:10.1016/j.neuroimage.2015.02.065'};
        ea_methods(options,parsestr,refs)
    end

    ea_importfssegmentations(options);
end

if options.prefs.fs.samseg.do
    % run samseg:

    allvols=[];
    for vol=1:length(presentfiles)
        allvols=[allvols,options.subj.coreg.anat.preop.(presentfiles{vol}),' '];
    end
    allvols=allvols(1:end-1); % rm last trailing space
    ea_mkdir([options.subj.subjDir,filesep,'fs',filesep,'samseg']);
    system([options.prefs.fs.dir,filesep,'bin',filesep,...
        'run_samseg',...
        ' --input ',allvols,...
        ' --output ',[options.subj.subjDir,filesep,'fs',filesep,'samseg'],...
        ' --pallidum-separate ']);
end
