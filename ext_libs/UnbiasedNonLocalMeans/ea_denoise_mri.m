function ea_denoise_mri(input,output,sigma)

basedir = [fileparts(mfilename('fullpath')), filesep];
ea_dispt(['Denoising: ',input]);
UNLM = ea_getExec([basedir, 'UnbiasedNonLocalMeans'], escapePath = 1);


if ~exist('sigma','var')
    sigma=5;
end

cmd=[UNLM,' --sigma ',num2str(sigma),' ',ea_path_helper(input),' ',ea_path_helper(output)];
ea_runcmd(cmd);

[pth,inputfname]=fileparts(input);

%% add methods dump:
cits={
    '   A. Tristan-Vega, V. Garcia Perez, S. Aja-Fenandez, and C.-F. Westin, Efficient and Robust Nonlocal Means Denoising of MR Data Based on Salient Features Matching, Computer Methods and Programs in Biomedicine. 2011.'
};

ea_methods(pth,[inputfname,' was denoised using an unbiased nonlocal means filtering approach (Tristan-Vega et al. 2011; https://www.nitrc.org/projects/unlmeans/)'],...
    cits);
ea_dispt([]);


