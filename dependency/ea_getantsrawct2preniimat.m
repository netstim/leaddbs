function tmat = ea_getantsrawct2preniimat(options)
% Gets the ANTS transformation from options.prefs.rawctnii_unnormalized to options.prefs.prenii_unnormalized
% and extracts the transformation matrix by calling_ea_antrsmat2mat. The
% matrix is read regardless of
%
% Returns: transformation matrix in LPI- (= RAS+)
%
% 2018, Andreas Husch, University of Luxembourg, Intenventional
% Neuroscience Group
try
    filepath = [options.root,options.patientname,filesep,stripext(options.prefs.rawctnii_unnormalized),'2',stripext(options.prefs.prenii_unnormalized),'_ants1.mat'];
    t=load(filepath);
    % The affine field name in tfields{1} differs depending on the ants call, its often
    % "AffineTransform_float_3_3", but alternativley "AffineTransform_double_3_3"
    % or "CompositeTransform_double_3_3" could occour. Note that composite
    % transforms would need further handling (multipling of the resulting
    % matrices) as they store multiple ants transforms in one file *wihthout* combining them directly.
    tfields = fieldnames(t);
    % affine         % fixed
    tmat=ea_antsmat2mat(t.(tfields{1}),t.(tfields{2}));
catch
    warning(['ea_getantsrawct2preniimat: Failure reading ' filepath]);
end

%% Nested helper functions
    function fn=stripext(fn)
        [~,fn]=fileparts(fn);
        
    end
end