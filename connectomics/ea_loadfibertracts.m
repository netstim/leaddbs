function [fibers,idx,voxmm,mat]=ea_loadfibertracts(cfile)

if strcmp(cfile(end-3:end),'.mat')
    fibinfo=load(cfile);
    if ~isfield(fibinfo,'ea_fibformat')
        fibinfo=ea_convertfibs2newformat(fibinfo,cfile);
    end

    fibers=fibinfo.fibers;
    idx=fibinfo.idx;
    if nargout>2
        try
            voxmm=fibinfo.voxmm;
        catch % assume voxel
            voxmm='vox';
        end
        mat=[];
        try
           mat=fibinfo.mat;
        end
    end
else
    ea_trk2ftr(cfile);
    [pth,fn,~]=fileparts(cfile);
    cfile=fullfile(pth,[fn,'.mat']);
    [fibers,idx,voxmm]=ea_loadfibertracts(cfile);
end



function ea_trk2ftr(cfile)

[pth,fn,~]=fileparts(cfile);
[hdr,trks]=ea_trk_read([cfile]);
ftr.vox=hdr.voxel_size;

%V=spm_vol([spm('dir'),filesep,'canonical',filesep,'avg152T2.nii']); % assume MNI152 official to be voxel space of tracts
V=spm_vol([ea_getearoot,'templates',filesep,'mni_hires.nii']);
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

function ftr=ea_convertfibs2newformat(fibinfo,cfile)

fn=fieldnames(fibinfo);
if isfield(fn,'normalized_fibers_mm')
    fibers=fibinfo.normalized_fibers_mm;
    voxmm='mm';
else
    fibers=eval(['fibinfo.',fn{1},';']);
    voxmm='mm'
end

c=size(fibers);
if c(1)<c(2)
    fibers=fibers';
end

ea_dispercent(0,'Converting fibers');
[idx,~]=cellfun(@size,fibers);

fibers=cell2mat(fibers);
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'
    ea_dispercent(cnt/length(idx));
    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
ea_dispercent(1,'end');

fibers=[fibers,idxv];

ftr.fibers=fibers;
ftr.idx=idx;
fibinfo.voxmm=voxmm;

ea_savefibertracts(cfile,fibers,idx,'mm')
