function mat=ea_getscrfmat(directory, varargin)

load([directory,'scrf',filesep,'scrf_instore.mat']);
mmat=ea_antsmat2mat(AffineTransform_float_3_3,fixed); % analytical solution

if ~isempty(varargin) && isstruct(varargin{1})...
        && varargin{1}.prefs.env.dev==1 % perform additional check    
    mat=ea_antsmat2mat_empirical(directory); % do an extra empirical check
    if sum(abs(mmat(:)-mat(:)))<1e-06 % precision
        mat=mmat;
    else % debug this case.
            keyboard
    end
else
    mat=mmat;
end