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
        if isfield(coreg_log.method,'CT')
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
                % Tranformation from BRAINSFit has been converted to MAT file
                transform_file_name = [options.subj.coreg.transform.CT.([prefix 'BaseName']) 'brainsfit.mat'];
                t = load(transform_file_name);
                tmat = ea_antsmat2mat(t.AffineTransform_double_3_3, t.fixed);

            elseif contains(coreg_log.method.CT, 'FLIRT')
                transform_file_name = [options.subj.coreg.transform.CT.forwardBaseName 'flirt.mat'];
                tmat_flirt = readmatrix(transform_file_name, 'FileType', 'text');
                %TODO check if add + 1 is nessesary, check the inverse case
                tmat = flirtmat2worldmatPaCER(tmat_flirt, postopct, options.subj.coreg.anat.preop.(options.subj.AnchorModality), 0);
                tmat = inv(tmat); % for some reasen the way flirtmat2worldmatPaCER is coded, it returns the inverse matrix by default
                % use inv() function to return inverse when queried (for example DiODe)
                if inverse
                    tmat = inv(tmat);
                end
            end
        else
            transform_file_name = [options.subj.coreg.transform.CT.([prefix 'BaseName']) 'ants.mat'];
            if isfile(transform_file_name)
                t = load(transform_file_name);
                % The affine field name in tfields{1} differs depending on the ants call, its often
                % "AffineTransform_float_3_3", but alternativley "AffineTransform_double_3_3"
                % or "CompositeTransform_double_3_3" could occour. Note that composite
                % transforms would need further handling (multipling of the resulting
                % matrices) as they store multiple ants transforms in one file *wihthout* combining them directly.
                tfields = fieldnames(t);
                % affine         % fixed
                tmat=ea_antsmat2mat(t.(tfields{1}),t.(tfields{2}));
            else
                tmat = eye(4);
                postopct = options.subj.coreg.anat.postop.CT;
            end
        end
    case 'rpostop'
        tmat = eye(4);
        postopct = options.subj.coreg.anat.postop.CT;
end
     