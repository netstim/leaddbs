function mat=ea_getscrfmat(directory)

load([directory,'scrf',filesep,'scrf_instore.mat']);
mmat=ea_antsmat2mat(AffineTransform_float_3_3,fixed); % analytical solution
prefs=ea_prefs;
if prefs.env.dev==1 % perform additional check    
    mat=ea_antsmat2mat_empirical(directory); % do an extra empirical check
    if sum(abs(mmat(:)-mat(:)))<1e-04 % precision
        mat=mmat;
    else % debug this case.
             keyboard
    end
else
    mat=mmat;
end