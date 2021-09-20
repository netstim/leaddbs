function [tmat,postopct] = ea_getrawct2preniimat(options,inverse)
% Gets the transformation from options.prefs.rawctnii_unnormalized to
% options.prefs.prenii_unnormalized and extracts the transformation matrix
%
% Returns: transformation matrix in LPI- (= RAS+)
%
% 2018, Andreas Husch, University of Luxembourg, Intenventional
% Neuroscience Group

if ~exist('inverse','var')
    inverse=0;
end

if ~inverse
    prefix = 'forward';
else
    prefix = 'inverse';
end

switch options.prefs.reco.mancoruse
    case 'postop'
        postopct = options.subj.preproc.anat.postop.CT;
        coreg_log = loadjson(options.subj.coreg.log.method);
        
        if contains(coreg_log.method.CT, 'ANTs')
            transform_file_name = [options.subj.coreg.transform.CT.([prefix 'BaseName']) 'ants.mat'];
            t = load(transform_file_name);
            % The affine field name in tfields{1} differs depending on the ants call, its often
            % "AffineTransform_float_3_3", but alternativley "AffineTransform_double_3_3"
            % or "CompositeTransform_double_3_3" could occour. Note that composite
            % transforms would need further handling (multipling of the resulting
            % matrices) as they store multiple ants transforms in one file *wihthout* combining them directly.
            tfields = fieldnames(t);
            % affine         % fixed
            tmat=ea_antsmat2mat(t.(tfields{1}),t.(tfields{2}));
            
        elseif contains(coreg_log.method.CT, 'BRAINSFit')
            transform_file_name = [options.subj.coreg.transform.CT.([prefix 'BaseName']) 'brainsfit.mat'];
            dataset_guess = {'/TransformGroup/0/TranformFixedParameters', '/TransformGroup/0/TranformParameters';...
                             '/TransformGroup/0/TransformFixedParameters', '/TransformGroup/0/TransformParameters';...
                             '/fixed', '/AffineTransform_double_3_3'};
            for i = 1:size(dataset_guess, 1)                
                try
                    reg2org.fixed = h5read(transform_file_name, dataset_guess{i,1});
                    reg2org.AffineTransform_float_3_3 = h5read(transform_file_name, dataset_guess{i,2});
                    break
                end
            end
            tmat = ea_antsmat2mat(double(reg2org.AffineTransform_float_3_3), reg2org.fixed);

        elseif contains(coreg_log.method.CT, 'FLIRT')
            transform_file_name = [options.subj.coreg.transform.CT.forwardBaseName 'flirt.mat'];
            tmat_flirt = dlmread(transform_file_name);
            %TODO check if add + 1 is nessesary, check the inverse case
            tmat = flirtmat2worldmatPaCER(tmat_flirt, options.subj.coreg.anat.preop.(options.subj.AnchorModality), postopct, false );
            % use inv() function to return inverse when queried (for example DiODe)
            if inverse
                tmat = inv(tmat);
            end
        end

    case 'rpostop'
        tmat = eye(4);
        postopct = options.subj.coreg.anat.postop.CT;
end
     