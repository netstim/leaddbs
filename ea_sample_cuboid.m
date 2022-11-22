function [cimat,reldist,mat]=ea_sample_cuboid(varargin)
% This function samples image points alongside a trajectory specified.
% usage: [cimat,~,mat]=ea_sample_cuboid(trajvox (trajectory in voxel space),options (regular lead options struct),volume_to_be_sampled_from,interpolation_mode (see spm_sample_vol),width_of_sampling,distance,spacing_of_sampling);

trajectory=varargin{1};
trajvector=-mean(diff(trajectory));
options=varargin{2};
vizz=0;

if vizz
    figure;
    hold on
end
    
switch options.subj.postopModality
    case 'MRI'
        if isfield(options.subj.coreg.anat.postop, 'cor_MRI') && isfile(options.subj.coreg.anat.postop.cor_MRI)
            niifn = options.subj.coreg.anat.postop.cor_MRI;
        else
            niifn = options.subj.coreg.anat.postop.ax_MRI;
        end
    case 'CT'
        niifn = options.subj.coreg.anat.postop.CT;
end

if nargin>=3 % custom sampling image specified.
    niifn=varargin{3};
end

if nargin>=4
    interp=varargin{4};
else
    interp=3;
end

if nargin>=5
    width=varargin{5};
else
    width=15;
end

if nargin>=6
    distance=varargin{6};
end

if nargin>=7
    spacing=varargin{7};
else
    spacing=1;
end

V=spm_vol(niifn);

top_vx=[trajectory(1:2,:),ones(2,1)]';
top_mm=V.mat*top_vx;
top_vx=top_vx(1:3,:)';
top_mm=top_mm(1:3,:)';

dvox=ea_pdist(top_vx);
dmm=ea_pdist(top_mm);
mm2vx=dmm/dvox; % -> 1 mm equals mm2vx voxels.
reldist=options.elspec.eldist/mm2vx; % real measured distance between electrodes in voxels.

ntrajvector=trajvector/norm(trajvector);
trajvector=trajvector*((reldist/10)/(norm(trajvector))); % normed trajvector.
if trajvector(3)<0
    trajvector=trajvector*-1; % now going from dorsal to ventral.
    ntrajvector=ntrajvector*-1;
end

% now set trajvector to a size that 10 of it will be the size of eldist.
if vizz
    plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3),'b')
end

startpt=trajectory(end,:);            % 3d point of starting line (5 vox below coord 1).

orth=null(trajvector);
orthx=orth(:,1)'*((reldist/10)/norm(orth(:,1))); % vector going perpendicular to trajvector by 0.5 mm in x dir.
orthy=orth(:,2)'*((reldist/10)/norm(orth(:,2))); % vector going perpendicular to trajvector by 0.5 mm in y dir.

xdim=width; % default is 15
ydim=width; % default is 15

if exist('distance','var')
    zdim=distance;
else
    zdim=150; % will be sum up to 5 times reldist (three between contacts and two at borders).
end

if orthx(1)>0
    orthx=-orthx;
end

if orthy(1)>0
    orthy=-orthy;
end

xv=-xdim:spacing:xdim;
yv=-ydim:spacing:ydim;

zv=1:spacing:zdim;
imat=nan(length(yv),length(xv),length(zv));

cnt=1;
coord2write=zeros(numel(imat),3);
coord2extract=zeros(numel(imat),3);

if vizz
    pcoord2extract=zeros(1,3);
end

zcnt=1;

for zz=zv
    xcnt=1;
    for xx=xv
        ycnt=1;
        for yy=yv
            pt=startpt+zz*trajvector;
            coord2extract(cnt,:)=[pt(1)+orthx(1)*xx+orthy(1)*yy; ...
                pt(2)+orthx(2)*xx+orthy(2)*yy; ...
                pt(3)+orthx(3)*xx+orthy(3)*yy]';
            coord2write(cnt,:)=[xcnt,ycnt,zcnt];
            ycnt=ycnt+1;
            cnt=cnt+1;
        end
        xcnt=xcnt+1;
    end
    zcnt=zcnt+1;
end

if vizz
    plot3(coord2extract(:,1),coord2extract(:,2),coord2extract(:,3),'r.'); 
    plot3(coord2write(:,1),coord2write(:,2),coord2write(:,3),'b.');    
end

imat(sub2ind(size(imat),coord2write(:,1),coord2write(:,2),coord2write(:,3)))=spm_sample_vol(V,coord2extract(:,1),coord2extract(:,2),coord2extract(:,3),interp);
cimat=squeeze(imat(:,:,:));

if nargout>2 % also export V.mat
    % put coord2extract to mm coordinates
    coord2extract=V.mat*[coord2extract';ones(1,size(coord2extract,1))];
    coord2write=[coord2write,ones(size(coord2write,1),1)];
    mat=mldivide(coord2write,coord2extract');
end
