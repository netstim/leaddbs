function tmat = ea_getantsrawct2preniimat(options)
% Gets the ANTS transformation from options.prefs.rawctnii_unnormalized to options.prefs.prenii_unnormalized
% and extracts the transformation matrix by calling_ea_antrsmat2mat. The
% matrix is read regardless of
%
% Returns: transformation matrix in LPI- (= RAS+)
%
% 2018, Andreas Husch, University of Luxembourg, Intenventional
% Neuroscience Group

directory=[options.root,options.patientname,filesep];
try
    load([directory 'ea_coregctmethod_applied.mat']);
    switch coregct_method_applied{end}
        case {'ea_coregctmri_ants','ea_coregctmri_ants_refine'}
            antsmts=dir([directory,ea_stripex(options.prefs.prenii_unnormalized),'2',ea_stripex(options.prefs.rawctnii_unnormalized),'*','.mat']);
            t=load([directory,antsmts(end).name]);
            % The affine field name in tfields{1} differs depending on the ants call, its often
            % "AffineTransform_float_3_3", but alternativley "AffineTransform_double_3_3"
            % or "CompositeTransform_double_3_3" could occour. Note that composite
            % transforms would need further handling (multipling of the resulting
            % matrices) as they store multiple ants transforms in one file *wihthout* combining them directly.
            tfields = fieldnames(t);
            % affine         % fixed
            tmat=ea_antsmat2mat(t.(tfields{1}),t.(tfields{2}));
        case 'ea_coregctmri_brainsfit'
            reg2org.fixed = h5read([directory 'postop_ct2anat_t1_brainsfit_Inverse.h5'],'/TransformGroup/0/TranformFixedParameters');
            reg2org.AffineTransform_float_3_3 = h5read([folder 'postop_ct2anat_t1_brainsfit_Inverse.h5'],'/TransformGroup/0/TranformParameters');
            tmat = ea_antsmat2mat(reg2org.AffineTransform_float_3_3,reg2org.fixed);
        case 'ea_coregctmri_fsl'
            %             tmat_reg2org = dlmread([folder 'anat_t12postop_ct_flirt1.mat']));
            disp(['Warning: Temporary fix to use DiODe algorithm with FLIRT. rpostop_ct is used so results may be slightly less accurate.'])
            ct = ct_reg;
    end
catch
    warning(['ea_getantsrawct2preniimat: Failure reading ' filepath]);
end