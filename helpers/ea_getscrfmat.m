function mat=ea_getscrfmat(directory)

load([directory,'scrf',filesep,'scrf_instore.mat']);
mmat=ea_antsmat2mat(AffineTransform_float_3_3,fixed); % analytical solution
mat=ea_antsmat2mat_empirical(directory); % do an extra empirical check
if sum(abs(mmat(:)-mat(:)))<1e-06 % precision
    mat=mmat;
else
    % keep empirical mat since surely correct
    if options.prefs.env.dev==1 % debug this case.
        keyboard
    end
end