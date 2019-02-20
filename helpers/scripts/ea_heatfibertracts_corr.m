function [fibsweighted,fibsin,fibsval,iaix]=ea_heatfibertracts_corr(cfile,roilist,patselection,vals,thresh,minpercent)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.

disp('ROI fiber analysis');

fibers=load(cfile);
fn=fieldnames(fibers);
try
    fibers=fibers.fibers;
catch
    fibers=fibers.(fn{1});
end

if ~iscell(roilist)
    roilist={roilist};
    vals={vals};
end

if ~exist('minpercent','var') % minimum of percentage fibers that need to be connected to a VTA.
    minpercent=0.2;
end

allroilist=cat(2,roilist{:});
%tree=KDTreeSearcher(fibers(:,1:3));
% load in all ROI
ea_dispercent(0,'Aggregating ROI');
cnt=1;
for roi=1:length(allroilist)
    if size(allroilist,2)==2 % left and right entered separately, combine.
        nii{roi,1}=ea_load_nii(allroilist{roi,1});

        [xx,yy,zz]=ind2sub(size(nii{roi,1}.img),(find(nii{roi,1}.img(:)>thresh)));
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XY=nii{roi,1}.mat*XYZvx;
        XYZmm{roi}=XY(1:3,:)';
           valsmm{roi}=nii{roi,1}.img((nii{roi,1}.img(:))>thresh);

        nii{roi,2}=ea_load_nii(allroilist{roi,2});
        [xx,yy,zz]=ind2sub(size(nii{roi,2}.img),find(nii{roi,2}.img(:)>thresh));
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XY=nii{roi,2}.mat*XYZvx;
        XYZmm{roi}=[XYZmm{roi};XY(1:3,:)'];
        
        valsmm{roi}=[valsmm{roi};nii{roi,2}.img((nii{roi,2}.img(:))>thresh)];
    else
        keyboard % this is not maintained and currently not functioning.
        nii{roi}=ea_load_nii(allroilist{roi});

        if exist('thresh','var')
            nii{roi}.img=double(nii{roi}.img>thresh);
        end
        keyboard
        [xx,yy,zz]=ind2sub(size(nii{roi}.img),find(nii{roi}.img));
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XYZmm{roi}=nii{roi}.mat*XYZvx;
        XYZmm{roi}=XYZmm{roi}(1:3,:)';
    end

    if ~exist('AllXYZ','var')
        AllXYZ=XYZmm{roi};
       % Allvals=vals{roi};
        
    else
        AllXYZ=[AllXYZ;XYZmm{roi}];
       % Allvals=[Allvals;vals{roi}];
    end

    if ismember(roi,patselection)
        XYZmmSel{cnt}=XYZmm{roi};
        niiSel{cnt,:}=nii{roi,:};

        if ~exist('AllXYZSel','var')
            AllXYZSel=XYZmmSel{cnt};
        else
            AllXYZSel=[AllXYZSel;XYZmm{cnt}];
        end
        cnt=cnt+1;
    end

    ea_dispercent(roi/length(allroilist));
end
ea_dispercent(1,'end');
ea_dispt('Selecting connected fibers');
% isolate fibers that are connected to any ROI:
AllXYZ=unique(round(AllXYZ),'rows');

[~,D]=knnsearch(AllXYZ(:,1:3),fibers(:,1:3),'Distance','chebychev');

in=D<(mean(nii{end}.voxsize)*2); % larger threshold here since seed vals have been rounded above. Also, this is just to reduce size of the connectome for speed improvements, so we can be more liberal here.

fibsin=fibers(ismember(fibers(:,4),unique(fibers(in,4))),:);
fibsval=[zeros(size(fibsin,1),length(patselection))]; % 5th column will add up values, 6th will take note how many entries were summed.
if size(fibsin,2)>4
   fibsin(:,5:end)=[];
end

AllXYZ=AllXYZSel; % now that connected fibers were selected, replace with selected VTAs only to avoid confusion.
XYZmm=XYZmmSel;
nii=niiSel;
% now color fibsin based on predictive value of improvement
ea_dispt('');

ea_dispercent(0,'Iterating ROI');
fibunique=unique(fibsin(:,4));
for roi=1:length(patselection)
    [IX,D]=knnsearch(XYZmm{roi}(:,1:3),fibsin(:,1:3),'Distance','chebychev');
    in=D<mean(nii{end}.voxsize);
    for fib=fibunique'
        fibsval(fibsin(:,4)==fib,roi)=sum(valsmm{roi}(IX(and(in,(fibsin(:,4)==fib)))));
    end
    ea_dispercent(roi/length(patselection));
end
ea_dispercent(1,'end');
ea_dispt('Correlating fibers with values');
cnt=1;
for group=1:length(roilist) % groups currently not implemented, should always be one within Lead-DBS
    thisgroupidx=cnt:(cnt+length(patselection))-1;
    cnt=cnt+length(patselection);

    % reduce to one entry per fiber:
    [~,fibidx,iaix]=unique(fibsin(:,4));
    fibsval=fibsval(fibidx,:);
    repvals=repmat(vals{group}',size(fibsval,1),1);

    [R]=corr(vals{group},fibsval','rows','pairwise','type','Spearman');
  
    fibsweighted=fibsin;

    fibsweighted=[fibsweighted,R(iaix)'];
end
