function [fibers,idx,voxmm,mat,vals]=ea_loadfibertracts(cfile, ref)
if ~exist('ref','var')
    ref = 'mni';
end

if endsWith(cfile, {'.trk', '.trk.gz'})
    ea_trk2ftr(cfile, ref, 1);
    cfile = replace(erase(cfile, '.gz'), '.trk', '.mat');
end

fibinfo = load(cfile);
if ~isfield(fibinfo,'ea_fibformat')
    ea_convertfibs2newformat(fibinfo,cfile);
    fibinfo = load(cfile);
end

fibers = fibinfo.fibers;
idx = fibinfo.idx;

if isfield(fibinfo,'vals')
    vals=fibinfo.vals;
else
    vals=ones(size(idx));
end

if nargout>2
    if isfield(fibinfo, 'voxmm')
        voxmm = fibinfo.voxmm;
    elseif any(fibers<0,'all')
        voxmm = 'mm';
    else % assume voxel
        voxmm = 'vox';
    end

    if isfield(fibinfo, 'mat')
        mat = fibinfo.mat;
    else
        mat = [];
    end
end


function ea_convertfibs2newformat(fibinfo,cfile)

disp('Converting fibers...');

fn=fieldnames(fibinfo);
if isfield(fibinfo,'normalized_fibers_mm')
    fibers=fibinfo.normalized_fibers_mm;
    voxmm='mm';
elseif isfield(fibinfo,'curveSegCell') % original Freiburg format
    fibers=fibinfo.curveSegCell;
    voxmm='vox';
    freiburgconvert=1;
elseif isfield(fibinfo,'fibs') % original Freiburg format
    fibers=fibinfo.fibs;
    voxmm='mm';
else
    fibers=eval(['fibinfo.',fn{1},';']);
    voxmm='mm';
end

c=size(fibers);
if c(1)<c(2)
    fibers=fibers';
end

[idx,~]=cellfun(@size,fibers);
fibers=cell2mat(fibers);
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'
    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
fibers=[fibers,idxv];

if exist('freiburgconvert','var')
    [pth, ~]=fileparts(cfile);
    prefs=ea_prefs('');
    if isempty(pth)
        b0fi=[prefs.b0];
    else
        b0fi=[pth,filesep,prefs.b0];
    end
    try
        ver=str2double(fibinfo.version(2:end));
    catch
        ver=1.1;
    end
    if ver<1.1
        disp('Flip fibers...');

        dim=getfield(spm_vol(b0fi),'dim');

        % Freiburg2World transform: xy-swap and y-flip
        tfibs=fibers;
        tfibs(:,1)=fibers(:,2);
        tfibs(:,2)=dim(2)+1-fibers(:,1);
        fibers=tfibs;
        clear tfibs
    end
    mat=getfield(spm_vol(b0fi),'mat');
    ea_savefibertracts(cfile,fibers,idx,voxmm,mat);
else
    ea_savefibertracts(cfile,fibers,idx,voxmm);
end
