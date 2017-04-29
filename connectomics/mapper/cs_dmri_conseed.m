function cs_dmri_conseed(dfold,cname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,space)
%# ea_load_nii




[sfile,roilist]=ea_handleseeds(sfile);


if isdeployed
    cbase=dfold;
else
    cbase=ea_getconnectomebase;
end
disp(['Command: ',cmd]);
switch cmd
    case 'seed'
        for s=1:length(sfile)
            map=ea_load_nii([cbase,'spacedefinitions',filesep,space]);
            cfile=[dfold,'dMRI',filesep,cname];
            
            if exist([cfile,filesep,'data.mat'],'file') % regular mat file
                if ~exist('fibers','var')
                    [fibers,fidx,voxmm,mat]=ea_loadfibertracts([cfile,filesep,'data.mat']);
                    if ~exist('fibers','var')
                        ea_error('Structural connectome file supplied in wrong format.');
                    end
                end
                
                redotree=0;
                ctype='mat';
            elseif exist([cfile,filesep,'data.fib.gz'],'file') % regular .fib.gz file
                
                ftr=track_seed_gqi([cfile,filesep,'data.fib.gz'],sfile{s});
                fibers=ftr.fibers;
                redotree=1;
                ctype='fibgz';
                
            else % connectome type not supported
                ea_error('Connectome file vanished or not supported!');
            end
            
            
            mapsz=size(map.img);
            
            seedfiles=sfile;
            
            map.img(:)=0;
            Vseed=ea_load_nii(seedfiles{s});
            
            maxdist=mean(abs(Vseed.voxsize))/2;
            
            Vseed.img(isnan(Vseed.img))=0;
            
            ixs=find(Vseed.img);
            % subtract nan values from these
            
            ixvals=Vseed.img(ixs);
            if sum(abs(ixvals-double(logical(ixvals))))<0.0001
                allbinary=1;
            else
                allbinary=0;
            end
            
            if ~allbinary || strcmp(ctype,'mat')
                
                [xx,yy,zz]=ind2sub(size(Vseed.img),ixs);
                XYZvx=[xx,yy,zz,ones(length(xx),1)]';
                clear ixs
                XYZmm=Vseed.mat*XYZvx;
                XYZmm=XYZmm(1:3,:)';
                clear Vseed
                if ~exist('tree','var') || redotree % only compute for first seed.
                    tree=KDTreeSearcher(fibers(:,1:3));
                end
                ids=rangesearch(tree,XYZmm,maxdist,'distance','chebychev');
                % select fibers for each ix
                ea_dispercent(0,'Iterating voxels');
                ixdim=length(ixvals);
                fiberstrength=zeros(size(fidx,1),1); % in this var we will store a mean value for each fiber (not fiber segment) traversing through seed
                fiberstrengthn=zeros(size(fidx,1),1); % counting variable to average strengths
                for ix=find(cellfun(@length,ids))'
                    % assign fibers on map with this weighted value.
                    %if ~isempty(ids{ix})
                        fibnos=unique(fibers(ids{ix},4)); % these fiber ids go through this particular voxel.
                        fiberstrength(fibnos)=fiberstrength(fibnos)+ixvals(ix);
                        fiberstrengthn(fibnos)=fiberstrengthn(fibnos)+1;
                    %end
                    ea_dispercent(ix/ixdim);
                    
                end
                nzz=~(fiberstrength==0);
                fiberstrength(nzz)=fiberstrength(nzz)./fiberstrengthn(nzz); % now each fiber has a strength mediated by the seed.
                ea_dispercent(1,'end');
                
                ea_dispercent(0,'Iterating fibers');
                cfibers=find(fiberstrength);
                cfibdim=length(cfibers);
                fcnt=1;
                for f=cfibers' % iterate through fibers that have assigned a nonzero value and paint to map.
                    allfibcs=fibers(fibers(:,4)==f,1:3);
                    
                    allfibcs=round(map.mat\[allfibcs,ones(size(allfibcs,1),1)]');
                    allfibcs(:,logical(sum(allfibcs<1,1)))=[];
                    topaint=sub2ind(mapsz,allfibcs(1,:),allfibcs(2,:),allfibcs(3,:));
                    map.img(topaint)=map.img(topaint)+fiberstrength(f);
                    ea_dispercent(fcnt/cfibdim);
                    fcnt=fcnt+1;
                end
                ea_dispercent(1,'end');
            else % if all is binary && using a .fib.gz file (i.e. all fibers go through seed already), can be much quicker.
                allfibcs=fibers(:,1:3);
                allfibcs=round(map.mat\[allfibcs,ones(size(allfibcs,1),1)]');
                topaint=sub2ind(mapsz,allfibcs(1,:),allfibcs(2,:),allfibcs(3,:));
                utopaint=unique(topaint);
                c=countmember(utopaint,topaint);
                map.img(utopaint)=c;
            end
            
            [~,fn]=fileparts(seedfiles{s});
            
            map.fname=fullfile(outputfolder,[fn,'_struc_',cmd,'.nii']);
            map.dt=[16,0];
            spm_write_vol(map,map.img);
            
        end
        
    case {'matrix','pmatrix'}
        
        for s=1:length(sfile)
            cfile=[dfold,'dMRI',filesep,cname];
            
            if exist([cfile,filesep,'data.mat'],'file') % regular mat file
                if ~exist('fibers','var')
                    [fibers,fidx,voxmm,mat]=ea_loadfibertracts([cfile,filesep,'data.mat']);
                    if ~exist('fibers','var')
                        ea_error('Structural connectome file supplied in wrong format.');
                    end
                end
                
                redotree=0;
                ctype='mat';
            elseif exist([cfile,filesep,'data.fib.gz'],'file') % regular .fib.gz file
                
                ftr=track_seed_gqi([cfile,filesep,'data.fib.gz'],sfile{s});
                fibers=ftr.fibers;
                redotree=1;
                ctype='fibgz';
                
            else % connectome type not supported
                ea_error('Connectome file vanished or not supported!');
            end
            
            
            
            seedfiles=sfile;
            
            Vseed{s}=ea_load_nii(seedfiles{s});
            
            maxdist{s}=mean(abs(Vseed{s}.voxsize))/2;
            
            Vseed{s}.img(isnan(Vseed{s}.img))=0;
            
            ixs{s}=find(Vseed{s}.img);
            % subtract nan values from these
            
            ixvals{s}=Vseed{s}.img(ixs{s});
            if sum(abs(ixvals{s}-double(logical(ixvals{s}))))<0.0001
                allbinary{s}=1;
            else
                allbinary{s}=0;
            end
            
            
            % now have all seeds and connectome - for each seed find fibers
            % connected to it:
            
            [xx,yy,zz]=ind2sub(size(Vseed{s}.img),ixs{s});
            XYZvx=[xx,yy,zz,ones(length(xx),1)]';
            XYZmm{s}=Vseed{s}.mat*XYZvx;
            XYZmm{s}=XYZmm{s}(1:3,:)';
            if ~exist('tree','var') || redotree % only compute for first seed.
                tree=KDTreeSearcher(fibers(:,1:3));
            end
            ids{s}=rangesearch(tree,XYZmm{s},maxdist{s},'distance','chebychev');
            
            % select fibers for each ix
            ea_dispercent(0,'Iterating voxels');
            ixdim=length(ixvals{s});
            fiberstrength{s}=zeros(size(fidx,1),1); % in this var we will store a mean value for each fiber (not fiber segment) traversing through seed
            fiberstrengthn{s}=zeros(size(fidx,1),1); % counting variable to average strengths
            for ix=find(cellfun(@length,ids{s}))' % iterate through voxels that found something
                % assign fibers on map with this weighted value.
                fibnos=unique(fibers(ids{s}{ix},4)); % these fiber ids go through this particular voxel.
                fiberstrength{s}(fibnos)=fiberstrength{s}(fibnos)+ixvals{s}(ix);
                fiberstrengthn{s}(fibnos)=fiberstrengthn{s}(fibnos)+1;
                ea_dispercent(ix/ixdim);
            end
            nzz=~(fiberstrength{s}==0);
            fiberstrength{s}(nzz)=fiberstrength{s}(nzz)./fiberstrengthn{s}(nzz); % now each fiber has a strength mediated by the seed.
            ea_dispercent(1,'end');
        end
        
        fiberstrength=cell2mat(fiberstrength);
        mat=zeros(length(sfile));
        
        for sxx=1:length(sfile)
            for syy=1:length(sfile)
                if sxx>syy
                    switch cmd
                        case 'matrix'
                            mask=logical((fiberstrength(:,sxx)>0).*(fiberstrength(:,syy)>0));
                        case 'pmatrix'
                            oix=1:length(sfile);
                            oix([sxx,syy])=[];
                            mask=fiberstrength(:,sxx)>0.*fiberstrength(:,syy)>0.*~(fiberstrength(:,oix)>0);
                    end
                    mat(sxx,syy)=sum(fiberstrength(mask,sxx)+fiberstrength(mask,syy))/2;

                elseif sxx==syy % also fill diagonal
                    mat(sxx,syy)=sum(fiberstrength(:,sxx));
                end
            end
        end
        mat=mat+mat'; % symmetrize matrix
        mat(logical(eye(length(sfile))))=mat(logical(eye(length(sfile))))/2;
        [~,fn]=fileparts(seedfiles{1});
        save(fullfile(outputfolder,[fn,'_struc_',cmd,'.mat']),'mat');
    otherwise
        warning('Structural connectivity only supported for seed / matrix / pmatrix commands.');
end




function C = countmember(A,B)
% COUNTMEMBER - count members
%
%   C = COUNTMEMBER(A,B) counts the number of times the elements of array A are
%   present in array B, so that C(k) equals the number of occurences of
%   A(k) in B. A may contain non-unique elements. C will have the same size as A.
%   A and B should be of the same type, and can be cell array of strings.
%
%   Examples:
%     countmember([1 2 1 3],[1 2 2 2 2])
%        -> 1     4     1     0
%     countmember({'a','b','c'},{'a','x','a'})
%        -> 2     0     0
%
%   See also ISMEMBER, UNIQUE, HISTC

% for Matlab R13 and up
% version 1.2 (dec 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0 (2005) created
% 1.1 (??): removed dum variable from [AU,dum,j] = unique(A(:)) to reduce
%    overhead
% 1.2 (dec 2008) - added comments, fixed some spelling and grammar
%    mistakes, after being selected as Pick of the Week (dec 2008)

error(nargchk(2,2,nargin)) ;

if ~isequal(class(A),class(B)),
    error('Both inputs should be the same class.') ;
end
if isempty(B),
    C = zeros(size(A)) ;
    return
elseif isempty(A),
    C = [] ;
    return
end

% which elements are unique in A,
% also store the position to re-order later on
[AU,j,j] = unique(A(:)) ;
% assign each element in B a number corresponding to the element of A
[L, L] = ismember(B,AU) ;
% count these numbers
N = histc(L(:),1:length(AU)) ;
% re-order according to A, and reshape
C = reshape(N(j),size(A)) ;

function ftr=track_seed_gqi(cfile,seedfile)


basedir = [ea_getearoot, 'ext_libs',filesep,'dsi_studio',filesep];
if ismac
    dsistudio = [basedir,'mac',filesep, 'dsi_studio.app',filesep,'Contents',filesep,'MacOS',filesep,'dsi_studio'];
elseif isunix
    ea_libs_helper([basedir, 'linux']);
    dsistudio = [basedir, 'linux',filesep,'dsi_studio'];
elseif ispc
    dsistudio = [basedir, 'win',filesep,'dsi_studio.exe'];
end


[pth,fn]=fileparts(seedfile);

cmd=[dsistudio,' --action=trk --source=',ea_path_helper(cfile),...
    ' --method=0',...
    ' --seed=',ea_path_helper(seedfile),...
    ' --seed_count=10000',...
    ' --output=',ea_path_helper([pth,filesep,'temp.mat'])];


err=ea_submitcmd(cmd);
if err
    ea_error(['Fibertracking with dsi_studio failed (error code=',num2str(err),').']);
end

% now store tract in lead-dbs format
ea_dispercent(0,'Converting fibers');
fibinfo=load([pth,filesep,'temp.mat']);
fibers=fibinfo.tracts;
idx=fibinfo.length';
clear fibinfo
fibers=fibers';


fibers(:,1)=78.0-fibers(:,1);
fibers(:,2)=76.0-fibers(:,2);
fibers(:,3)=-50.0+fibers(:,3);

clear length
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

ftr.fourindex=1;
ftr.ea_fibformat='1.0';
ftr.fibers=fibers;
ftr.idx=idx;
delete([pth,filesep,'temp.mat']);

