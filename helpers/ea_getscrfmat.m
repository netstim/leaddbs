function mat=ea_getscrfmat(options)

load(options.subj.brainshift.transform.instore);

mmat = ea_antsmat2mat(AffineTransform_float_3_3,fixed); % analytical solution

prefs = ea_prefs;
if prefs.env.dev==1 % perform additional check    
    mat = ea_antsmat2mat_empirical(options); % do an extra empirical check
    if sum(abs(mmat(:)-mat(:)))<1e-04 % precision
        mat = mmat;
    else % debug this case.
    	keyboard
    end
else
    mat=mmat;
end