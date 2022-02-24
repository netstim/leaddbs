function nii=ea_exportgenemap(geneidx,smoothkernel,intensnorm,outputfolder,mirror,resolution,db)
% Function to export nifti image for a particular gene expression from the
% Allen institute database.
% (c) Andreas Horn 2021


prefs=ea_prefs;

if exist('db','var') && ~isempty(db)
    load(db);
else
    load([prefs.genetics.dbdir,'genedb.mat']);
end
if ~exist('smoothkernel','var')
    smoothkernel=3;
end
if ~exist('outputfolder','var')
    outputfolder=pwd;
end
if ~exist('mirror','var')
   mirror=1;
end
if ~exist('intensnorm','var')
    intensnorm='z';
end
if ~exist('resolution','var')
   resolution='222'; 
else
    switch resolution
        case {'2','222','2 mm','2mm'}
            resolution='222';
        case {'1','111','1 mm','1mm'}
            resolution='111';
        case {'0.5','05','5','555','0.5mm','0.5 mm'}
            resolution='555';
        otherwise
            resolution='222';
    end
end

if ~iscell(geneidx)
    geneidx={geneidx};
end
genedb=double(genedb);
load([prefs.genetics.dbdir,'geneinfo.mat'],'geneinfo');
for g=1:length(geneidx)
    id{g}=[];
    tableheaders=fieldnames(geneinfo);
    for header=[2,4,5]
        id{g}=ismember(geneinfo.(tableheaders{header}),geneidx{g});
        if any(id{g})
            break
        end
    end
    if ~any(id{g})
        ea_error(['Could not find gene: ',geneidx{g},'.']);
    end
    id{g}=find(id{g});

    
    nii=ea_load_nii(fullfile(ea_getearoot,'templates','spacedefinitions',[resolution,'.nii.gz']));
    linidx=find(nii.img(:));
    [xx,yy,zz]=ind2sub(size(nii.img),linidx);
    XYZ=[xx,yy,zz];
    
    if mirror
        useidxmm=[idxmm,idxmm];
        useidxmm(1,1:(size(useidxmm,2)/2))=-useidxmm(1,1:(size(useidxmm,2)/2)); % flip half of duplicated indices
    else
        useidxmm=idxmm;
    end

    if size(useidxmm,1)==3
        useidxmm=[useidxmm;ones(1,size(useidxmm,2))];
    end

    querypoints=nii.mat\useidxmm;
    querypoints=querypoints(1:3,:)';
 
    F=scatteredInterpolant(querypoints,...
        repmat(ea_nanmean((genedb(:,id{g})),2),double(logical(mirror))+1,1),...
        'natural');
    switch intensnorm
        case ''
            nii.img(linidx)=F(XYZ);
        case 'k'
            nii.img(linidx)=ea_normal(F(XYZ));
        case 'z'
            nii.img(linidx)=ea_nanzscore(F(XYZ));
    end
    nii.dt=[16 0];
    nii.fname=fullfile(outputfolder,[geneidx{g},'.nii']);
    ea_write_nii(nii);
    
    if smoothkernel
        spm_smooth(fullfile(outputfolder,[geneidx{g},'.nii']),fullfile(outputfolder,['s',geneidx{g},'.nii']),[smoothkernel smoothkernel smoothkernel]);
        movefile(fullfile(outputfolder,['s',geneidx{g},'.nii']),fullfile(outputfolder,[geneidx{g},'.nii']));
    end
end

