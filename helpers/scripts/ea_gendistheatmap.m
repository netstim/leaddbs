function ea_gendistheatmap(M,regno,pts,fovimg)
if ~exist('fovimg','var')
    fovimg=[ea_space,'atlas.nii'];
end
if ~exist('regno','var') || isempty(regno)
    regno=1;
end
if ~exist('pts','var') || isempty(pts)
    pts=1:length(M.patient.list);
end

if size(M.clinical.vars{regno},2)==2
    bilateral=1;
    sides=1:2;
elseif size(M.clinical.vars{regno},2)==1
    bilateral=0;
    sides=1;
else
    ea_error('This is only supported for variables with 1 entry per patient or pr hemisphere');
end

fovimg=ea_load_nii(fovimg);
distimg=fovimg;
for side=sides
    switch side
        case 1
            if bilateral
                sidestr='_rh';
            else
                sidestr='';
            end
        case 2
            sidestr='_lh';
    end
    
    fovimg.fname=[M.root,M.clinical.labels{regno},'_dist_heatmap',sidestr,'.nii'];
    distimg.fname=[M.root,M.clinical.labels{regno},'_dist_heatmap',sidestr,'_distance.nii'];

    [xx,yy,zz]=ind2sub(size(fovimg.img),1:numel(fovimg.img));
    XYZ=[xx;yy;zz;ones(1,length(xx))];
    XYZ=fovimg.mat*XYZ;
    XYZ=XYZ(1:3,:)';
    
    for pt=1:length(pts)
        if bilateral
            switch side
                case 1
                    acs(pt,:)=mean(M.elstruct(pts(pt)).coords_mm{1}(logical(M.S(pts(pt)).activecontacts{1}(1:4)),:),1);
                case 2
                    acs(pt,:)=mean(M.elstruct(pts(pt)).coords_mm{2}(logical(M.S(pts(pt)).activecontacts{2}(1:4)),:),1);
            end
        else % average points
            try
                acs(pt,:)=mean([mean(M.elstruct(pts(pt)).coords_mm{1}(logical(M.S(pts(pt)).activecontacts{1}(1:4)),:),1);...
                    ea_flip_lr_nonlinear(mean(M.elstruct(pts(pt)).coords_mm{2}(logical(M.S(pts(pt)).activecontacts{2}(1:4)),:),1))]);
            catch
                keyboard
            end
        end
    end
    
    % figure, plot3(acs(:,1),acs(:,2),acs(:,3),'r*');
    % axis equal
    ea_dispercent(0,'Iterating voxels');
    dimen=length(XYZ);
    N=length(M.patient.list(pts));
    I=M.clinical.vars{regno}(pts,side);
    chunk=200000;
    for vx=1:chunk:dimen
        
        if (vx+(chunk-1))>dimen
            chunk=(dimen-vx)+1;
        end
        %D=squareform(pdist([XYZ(vx,:);acs]));
        %D=-(pdist([XYZ(vx,:);acs]));
        D=pdist2(acs,XYZ(vx:vx+(chunk-1),:));
        %D=-D(1:N)';
        %D=1./exp(D);
        fovimg.img(vx:(vx+chunk-1))=corr(D,I,'rows','pairwise','type','Pearson');
        %        distimg.img(vx)=nansum(D);
        %distimg.img(vx:(vx+chunk-1))=nanmean(D);
                D=pdist2(mean(acs,1),XYZ(vx:vx+(chunk-1),:));
        distimg.img(vx:(vx+chunk-1))=D;
        %         b=glmfit(D,I);
        %         if isnan(b)
        %             keyboard
        %         end
        %
        %         fovimg.img(vx)=b(2)/overone(abs(b(1)));
        ea_dispercent(vx/dimen);
    end
    ea_dispercent(1,'end');
    
    
    fovimg.dt=[64,0];
    distimg.dt=[64,0];
    ea_write_nii(fovimg);
    ea_write_nii(distimg);
end

function val=overone(val)
if val<1
    val=1;
end
