function ea_generate_probmaps(varargin)
% function that will generate atlas maps from an image based on single
% points.
%
% inputs: pts in mm space or image defining a probability map
% labels for each point as a cellstring
% src images to use as a cellstring
% c 2016 Andreas Horn

% parameters:
othresh=0.9;
numpasses=1;
mrun=50;
mildnessfactor=2000; % increase to get less noisy but less specific images.
useweights=0;

if nargin>=3
    pts=varargin{1};
    label=varargin{2};
    srcs=varargin{3};
    directory=varargin{4};

else
%     srcs={'bbt2sinc.nii','bbt1sinc.nii','bbpdsinc.nii','bbrlxsinc.nii'};
%
%     pts=ea_readcsv('GPeGPi.fcsv');
%     pts=pts(1:3,:)';
%     label='GP_onepass';
end

for src=1:length(srcs)
    % ea_reslice_nii(srcs{src},['hd',srcs{src}],[0.22 0.22 0.22]);
    try
        S{src}=ea_load_nii(srcs{src});
    catch
        keyboard
    end

end


%% part one: generate multiple points for estimation of mahalanobis-distribution based on this single point:

if ~ischar(pts)
    pts=[pts,ones(size(pts,1),1)];
    pts=pts';
    pts=S{1}.mat\pts;
    pts=pts(1:3,:)';
    pts=round(pts);
    for dim=1:3
        pts(pts(:,dim)>size(S{1}.img,dim),:)=[];
    end
    for pass=1:numpasses

        in=nan(size(S{1}.img));

        for lab=1:size(pts,1)
            in(pts(lab,1),pts(lab,2),pts(lab,3))=1;
        end
        %ea_dispercent(0,['Region growing ',labs{lab}]);
        %ea_dispercent(0,[label,', pass number ',num2str(pass)]);
        for maxrun=1:mrun

            %% determine which voxel to assign next:
            assigned=find(~isnan(in));


            %ea_dispercent(size(assigned,1)/numel(in));

            s=size(in);
            [c1{1:3}]=ndgrid(1:3);
            c2(1:3)={2};
            offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});


            toassign=zeros(27,size(assigned,1));
            for a=1:size(assigned,1)
                thisassign=assigned(a)+offsets;
                toassign(:,a)=thisassign(:);
            end
            toassign=unique(toassign(:));

            % delete entries already assigned
            d=ismember(toassign,assigned);
            toassign(d,:)=[];

            if max(toassign(:))>numel(in)
                toassign(toassign>numel(in))=[];
            end

            if min(toassign(:))<0
                toassign(toassign<0)=[];
            end

            %% calculate weighted average value of assigned voxels so far:

            if pass==1
                W=in(assigned);
                W=W/sum(W);
                meanvalue=zeros(1,length(srcs));
                for src=1:length(srcs)
                    A=S{src}.img(assigned);
                    meanvalue(src) = mean(W.'*A,2);
                end
            end


            %% calculate values of toassigned voxel
            assignvalue=zeros(size(toassign,1),length(srcs));
            for src=1:length(srcs)
                assignvalue(:,src)=S{src}.img(toassign);
            end

            %% calculate and assign similarity:


            s=size(assignvalue,1);
            avals=zeros(1,s);
            % normalize data values:

            vals=[meanvalue;assignvalue];
            vals=zscore(vals);
            meanvalue=vals(1,:);
            assignvalue=vals(2:end,:);
            % end normalize
            for a=1:s
                sim=1/exp(0.05*ea_pdist([meanvalue;assignvalue(a,:)]));
                avals(a)=sim;
            end

            [mv,ix]=max(avals);
            if length(assigned)<2000
            thresh=mean(in(assigned))-1*std(in(assigned));
            end
            if thresh==1;
                thresh=othresh;
            end

            if mv<thresh
                break
            end

            toadd=avals>thresh;
            if ~sum(toadd)<mildnessfactor
                in(toassign(toadd))=avals(toadd);
            else

                [~,toadd]=sort(avals);
                try
                    in(toassign(toadd(end-mildnessfactor:end)))=avals(toadd(end-mildnessfactor:end));
                catch % if toadd has not mildnessfactor entries, just use all entries.
                    in(toassign(toadd))=avals(toadd);
                end
            end

            disp(['Threshold is ',num2str(thresh),', voxels included: ',num2str(sum(toadd)),'.']);

            %ea_dispercent(maxrun/mrun);
        end

        W=in(assigned);
        %W=W.^4;

        W=W/sum(W);

        meanvalue=zeros(1,length(srcs));
        for src=1:length(srcs)
            A=S{src}.img(assigned);
            meanvalue(src) = mean(W.'*A,2);
        end
        %ea_dispercent(1,'end');
    end
    %ea_dispercent(1,'end');

    labout=S{1};
    labout.img=in;
    labout.dt(1) = 16;
    labout.fname=[directory,label,'_firstlevel.nii'];
    spm_write_vol(labout,labout.img);

    % generate seed points for next part:
    [xx,yy,zz]=ind2sub(size(in),find(in>(ea_nanmean(in(:))+2*ea_nanstd(in(:)))));
    pts=[xx,yy,zz];

else
    nii=ea_load_nii(pts);
    [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>(ea_nanmean(nii.img(:))+2*ea_nanstd(nii.img(:)))));
    if isempty(xx)
            [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>(ea_nanmean(nii.img(:))+1*ea_nanstd(nii.img(:)))));

        if isempty(xx)
            [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>(ea_nanmean(nii.img(:)))));
        end
    end
    pts=[xx,yy,zz];
    if length(pts)>100000
        warning('Too many data points for algorithm submitted. Using 100000 random points from distribution.');
        pts=pts(round(linspace(1,length(pts),100000)),:);
    end
end

%% part two: iterate the whole brain by samples from this first estimate.
voxcnt=length(pts);
disp(['Included ',num2str(voxcnt),' to the ',label,'.']);
Xprob=zeros([size(S{1}.img)]);

%% determine covariance structure of given pointset
ixes=sub2ind(size(S{src}.img),pts(:,1),pts(:,2),pts(:,3));
for src=1:length(srcs)
    if useweights
        Aprob{src}=nan([size(S{1}.img)]);
    end
    profile(:,src)=S{src}.img(ixes);
end
up=mean(profile,1);
stp=std(profile,1);

% detemine V2. Delete if ~needed.
mask=S{1}.img;
mask(:)=0;
mask(ixes)=1;
mask=logical(mask);
unprofile=zeros(sum(~mask(:)),length(srcs));

for src=1:length(srcs)
    allup=mean(S{src}.img(:));
    allstd=std(S{src}.img(:));
    unprofile(:,src)=S{src}.img(~mask);
    vtwo(src)=mahal((up(src)-allup)/allstd,(unprofile(:,src)-allup)/allstd);
    vone(src)=(stp(src)/up(src)); % how constant are values in different acquisitions in the same tissue?
end
%end delete.
vtwo=vtwo./sum(vtwo);
vone=vone./sum(vone);
v=vone.*vtwo;
v=v./sum(v);

% equalize v
v=v.^(1/6);
v=v./sum(v); % sum of v is 1.

disp(['Mean values of profile for ',label,':']);
disp(num2str(up));
disp('Standard deviations:');
disp(num2str(stp));

metrics=[voxcnt,up,stp,vone,vtwo];
save([directory,'Metrics_',label,num2str(length(srcs))],'metrics');

for src=1:length(srcs);
    profile(:,src)=(profile(:,src)-up(src))/stp(src);
end

ea_dispercent(0,'Iterating voxels');
dimens=numel(Xprob(:));
chunk=5000000;
for ind=1:chunk:dimens

    if ind+chunk-1>dimens
       chunk=dimens-ind+1;
    end
    thisindprofile=zeros(chunk,length(srcs));
    for src=1:length(srcs)

        thisindprofile(:,src)=S{src}.img(ind:ind+chunk-1);
    end

    ea_dispercent(ind/dimens);

    %Xprob(ind:ind+slab-1)=1/exp(mahal(thisindprofile,profile));
    for src=1:length(srcs);
        thisindprofile(:,src)=(thisindprofile(:,src)-up(src))/stp(src);
        if useweights
            Aprob{src}(ind:ind+chunk-1)=mahal(thisindprofile(:,src),profile(:,src));
        end

    end
    if ~useweights

        Xprob(ind:ind+chunk-1)=mahal(thisindprofile,profile);

    end
end

labout=S{1};
if ~useweights
    labout.img=1/exp(0.01*Xprob);
else
    labout.img(:)=0;
    for src=1:length(srcs);
        labout.img=labout.img+((1/exp(0.1*Aprob{src}))*v(src));
    end
end
labout.dt(1) = 16;
labout.fname=[directory,label,'_secondlevel.nii'];
spm_write_vol(labout,labout.img);
matlabbatch{1}.spm.spatial.smooth.data = {labout.fname};
matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});


function c=ea_readcsv(pth)
fid=fopen(pth);
C=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %s %s %s','commentStyle', '#','delimiter', ',');
fclose(fid);
for coord=1:length(C{1})
   c(:,coord)=[C{2}(coord);C{3}(coord);C{4}(coord);1];
end
