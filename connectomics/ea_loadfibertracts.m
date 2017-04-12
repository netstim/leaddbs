function [fibers,idx,voxmm,mat]=ea_loadfibertracts(cfile)

if ~strcmp(cfile(end-3:end),'.mat')
    ea_trk2ftr(cfile);
    [pth,fn,~]=fileparts(cfile);
    cfile=fullfile(pth,[fn,'.mat']);
end

fibinfo=load(cfile);
if ~isfield(fibinfo,'ea_fibformat')
    ea_convertfibs2newformat(fibinfo,cfile);
    fibinfo=load(cfile);
end
fibers=fibinfo.fibers;
idx=fibinfo.idx;
if nargout>2
    try
        voxmm=fibinfo.voxmm;
    catch % assume voxel
        voxmm='vox';
    end
    try
        mat=fibinfo.mat;
    catch
        mat=[];
    end
end


function ea_trk2ftr(cfile)

[pth,fn,~]=fileparts(cfile);
[hdr,trks]=ea_trk_read([cfile]);
ftr.vox=hdr.voxel_size;

%V=spm_vol([spm('dir'),filesep,'canonical',filesep,'avg152T2.nii']); % assume MNI152 official to be voxel space of tracts
%V=spm_vol([ea_space,'t2.nii']);
for i=1:length(trks)
    ftr.curveSegCell{i}=[trks(i).matrix(:,1)/ftr.vox(1),trks(i).matrix(:,2)/ftr.vox(2),trks(i).matrix(:,3)/ftr.vox(3)];

    %% use conversion from DSI_studio:
    ftr.curveSegCell{i}(:,1)=78.0-ftr.curveSegCell{i}(:,1);
    ftr.curveSegCell{i}(:,2)=76.0-ftr.curveSegCell{i}(:,2);
    ftr.curveSegCell{i}(:,3)=-50.0+ftr.curveSegCell{i}(:,3);

    %% use conversion from a file:
%     ftr.curveSegCell{i}=[ftr.curveSegCell{i},ones(size(ftr.curveSegCell{i},1),1)]';
%     ftr.curveSegCell{i}=V.mat*ftr.curveSegCell{i};
%     ftr.curveSegCell{i}=ftr.curveSegCell{i}(1:3,:)';
end

save(fullfile(pth,[fn,'.mat']),'-struct','ftr','-v7.3');
delete(cfile);

function ea_convertfibs2newformat(fibinfo,cfile)

display('Converting fibers...');

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
    ver=str2double(fibinfo.version(2:end));
    if ver<1.1
        display('Flip fibers...');

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
