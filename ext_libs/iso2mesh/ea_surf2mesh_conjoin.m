function [node,elem,face,success]=ea_surf2mesh_conjoin(v,f,p0,p1,keepratio,maxvol,batchno,precision,options,stimname,side)
%
% [node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,regions,holes,forcebox)
%
% create quality volumetric mesh from isosurface patches
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/24
%
% input parameters:
%      v: input, isosurface node list, dimension (nn,3)
%         if v has 4 columns, the last column specifies mesh density near each node
%      f: input, isosurface face element list, dimension (be,3)
%      p0: input, coordinates of one corner of the bounding box, p0=[x0 y0 z0]
%      p1: input, coordinates of the other corner of the bounding box, p1=[x1 y1 z1]
%      keepratio: input, percentage of elements being kept after the simplification
%      maxvol: input, maximum tetrahedra element volume
%      regions: list of regions, specifying by an internal point for each region
%      holes: list of holes, similar to regions
%      forcebox: 1: add bounding box, 0: automatic
%
% outputs:
%      node: output, node coordinates of the tetrahedral mesh
%      elem: output, element list of the tetrahedral mesh
%      face: output, mesh surface element list of the tetrahedral mesh
%             the last column denotes the boundary ID
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fprintf(1,'generating tetrahedral mesh from closed surfaces ...\n');

exesuff=getexeext;
tMax = 300;
if(keepratio>1 || keepratio<0)
   warn(['The "keepratio" parameter is required to be between 0 and 1. '...
         'Your input is out of this range. surf2mesh will not perform '...
	 'simplification. Please double check to correct this.']);
end

% first, resample the surface mesh with cgal
if(keepratio<1-1e-9 && ~iscell(f))
	fprintf(1,'resampling surface mesh ...\n');
	[no,el]=meshresample(v(:,1:3),f(:,1:3),keepratio);
	el=unique(sort(el,2),'rows');
else
	no=v;
	el=f;
end

regions=[];
holes=[];

if(size(regions,2)>=4 && ~isempty(maxvol))
    warning('you specified both maxvol and the region based volume constraint,the maxvol setting will be ignored');
    maxvol=[];
end

dobbx=0;

% dump surface mesh to .poly file format
if(~iscell(el) && ~isempty(no) && ~isempty(el))
	saveoff(no(:,1:3),el(:,1:3),mwpath('post_vmesh.off'));
end

deletemeshfile(mwpath('post_vmesh.mtr'));

sweeptempdir;
savesurfpoly(no,el,holes,regions,p0,p1,mwpath('post_vmesh.poly'),dobbx);

moreopt='';
if(size(no,2)==4)
	moreopt=[moreopt ' -m '];
end
% call tetgen to create volumetric mesh
deletemeshfile(mwpath('post_vmesh.1.*'));
fprintf(1,'creating volumetric mesh from a surface mesh ...\n');

% write a protocol:
if ~exist([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname],'file')
    mkdir([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname]);
end
protocol_path=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'ea_genvat_horn_output_',num2str(side),'.txt'];
protocol=fopen(protocol_path,'a');
fprintf(protocol,'\n%s\n%s\n%s\n','================================',['BATCH #',num2str(batchno),', PRECISION = ',num2str(precision)],'================================');
fprintf(protocol,'%s\n%s\n%s\n','################################','################################','################################');

ea_delete(mwpath('post_vmesh.1.node'));
ea_kill('name', ['tetgen',exesuff]);

for p=1:2
    if p==1
        padd ='';
    elseif p==2
        padd='-p/0.02';
    end

    for tolerance=[8,10,5,4,12]
        fprintf(protocol,'\n\n%s\n%s\n',['ATTEMPT WITH P ''',padd,''' & TOLERANCE = 10^-',num2str(tolerance),' mm'],'-------------------------------------');

        % Start tetgen in background
        cmd = ['"' mcpath('tetgen') exesuff '" -A -T1e-',num2str(tolerance),' -pq1/0 ',padd,' -a -Y ' num2str(maxvol) ' ' moreopt ' "' mwpath('post_vmesh.poly') '" &'];
        system(cmd);

        if ea_isprocess('name', ['tetgen',exesuff]) % Check if tetgen started
            fprintf(protocol,'%s\n','Tetgen job started.');
            tStart=tic;

            % Break the 'tolerance' loop and toc if tetgen already finished so quickly
            if exist(mwpath('post_vmesh.1.node'),'file')
                if ~exist(mwpath('post_vmesh_skipped.node'),'file')
                    tEnd = toc(tStart);
                    break;
                else
                    ea_delete(mwpath('post_vmesh.1.node'));
                    ea_delete(mwpath('post_vmesh_skipped.node'));
                end
            else
                % Monitor if tetgen generated the mesh
                while ~exist(mwpath('post_vmesh.1.node'),'file')
                    % Break the while loop monitor if tetgen terminated itself
                    if ~ea_isprocess('name', ['tetgen',exesuff])
                        fprintf(protocol,'%s\n',['TERMINATED: Tetgen job terminated itself (because it failed) after ',num2str(tEnd),' seconds.']);
                        break
                    end

                    pause(2);
                    tEnd = toc(tStart);

                    % Break the while loop monitor if tetgen ran timeout
                    if tEnd > tMax
                        ea_kill('name', ['tetgen',exesuff]);
                        fprintf(protocol,'%s\n',['TERMINATED: Tetgen job has reached the time limit of ',num2str(tMax),' seconds and was killed by Lead-DBS.']);
                        break
                    end
                end

                % Break the 'tolerance' loop if tetgen already generated the mesh
                if exist(mwpath('post_vmesh.1.node'),'file')
                    if ~exist(mwpath('post_vmesh_skipped.node'),'file')
                        break;
                    else
                        ea_delete(mwpath('post_vmesh.1.node'));
                        ea_delete(mwpath('post_vmesh_skipped.node'));
                    end
                end
            end
        else
            fprintf(protocol,'%s\n','Tetgen job failed to start.');
            ea_kill('name', ['tetgen',exesuff]);
        end
    end

    % Break the 'ppad' loop if tetgen already generated the mesh
    if exist(mwpath('post_vmesh.1.node'),'file')
        if ~exist(mwpath('post_vmesh_skipped.node'),'file')
            break;
        else
            ea_delete(mwpath('post_vmesh.1.node'));
            ea_delete(mwpath('post_vmesh_skipped.node'));
        end
    end
end

if exist(mwpath('post_vmesh.1.node'),'file') && ~exist(mwpath('post_vmesh_skipped.node'),'file')
    % read in the generated mesh
    success=1;
    fprintf(protocol,'%s\n','TETGEN JOB SUCCESSFUL.');
    fprintf(protocol,'%s\n',['Tetgen job took ',num2str(tEnd),' seconds to complete.']);

    [node,elem,face]=readtetgen(mwpath('post_vmesh.1'));
else
    success=0;
end
fclose(protocol);


function sweeptempdir
 file=mwpath('post_vmesh.poly');
 pth=fileparts(file);
 warning('off');
 delete([pth,filesep,'*']);
 warning('on');
