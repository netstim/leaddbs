function ea_conformspaceto(spacefn,toreslicefn,interp,mask,newfn,headermod,force)

% set headermod to 1 by default
if ~exist('headermod', 'var')
    headermod = 1;
end
if ~exist('force','var')
    force=0;
end

sphdr = ea_open_vol(spacefn);
tohdr = ea_open_vol(toreslicefn);

if ~isequal(sphdr.mat,tohdr.mat) || force % volumes have different dimensions & hdr matrices.

    flags.mean=0;
    flags.which=1;

    if exist('interp','var') && ~isempty(interp)
        flags.interp = interp;
    end

    if exist('mask','var') && ~isempty(mask)
        flags.mask = mask;
    else
        flags.mask = 0;
    end

    if exist('newfn','var') && ~isempty(newfn)
        flags.prefix = 'r';
    else
        flags.prefix = '';
    end

    spm_reslice({spacefn,toreslicefn}, flags);

    if ~isempty(flags.prefix)
        [pth, fn, ext] = fileparts(toreslicefn);
        movefile(fullfile(pth, [flags.prefix,fn,ext]), newfn);
        toreslicefn = newfn;
    end

    nii = ea_load_nii(toreslicefn);
    delete(toreslicefn);
    ea_write_nii(nii);
end

if headermod
    % make sure headers of images are exactly identical (also corrects for qform/sform issues).
    sp = ea_load_untouch_nii(spacefn);
    tr = ea_load_untouch_nii(toreslicefn);
    sp.img = eval([class(tr.img),'(tr.img);']); % make sure to save data in same class as used before
    sp.hdr.dime.bitpix = tr.hdr.dime.bitpix;
    sp.hdr.dime.scl_slope = tr.hdr.dime.scl_slope;
    sp.hdr.dime.datatype = tr.hdr.dime.datatype; % keep datatype of original image.
    ea_save_untouch_nii(sp, toreslicefn);
end
