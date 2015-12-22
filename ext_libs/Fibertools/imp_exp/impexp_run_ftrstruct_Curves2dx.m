function out = impexp_run_ftrstruct_Curves2dx(job,cmd)
% Query an ftr_struct and extract some fiber sets in a format suitable for
% import by fancy_render.
% Structure of job:
% job.ftr.ftrname{1} - file name
%     or .ftrstruct  - ftr struct
% job.fibers.all     - all fibers
%        or .index   - vector of fiber indices
%        or .names   - cell string of fiber names
% Optional argument cmd:
% 'vout' - do not run job, return virtual outputs
% Output:
% out.dxfiles - cell string of output file names

switch cmd,
    case 'vout',
        % don't explicitly check cmd
        out            = cfg_dep;
        out.sname      = 'Exported fibers';
        out.src_output = substruct('.','dxfiles');
        out.tgt_spec   = cfg_findspec({{'filter', 'any', 'strtype', 'e'}});
    case 'run'
        % get and check ftrstruct
        try
            if isfield(job.ftr, 'ftrname')
                job.ftr.ftrstruct = load(job.ftr.ftrname{1});
            end
            [ver, errstr] = ftrstruct_query(job.ftr.ftrstruct, 'getVer');
            if ~isempty(errstr)
                error('ftrstruct_Curves2dx:load_ftrstruct', errstr);
            end
        catch le
            rethrow(le);
        end

        % get and check fiberNames
        [names errstr] = ftrstruct_query(job.ftr.ftrstruct, 'fiberNames');
        if ~isempty(errstr)
            error('ftrstruct_Curves2dx:get_fiberNames', errstr);
        end
        if isfield(job.fibers, 'all')
            job.fibers.names = names;
        elseif isfield(job.fibers, 'index')
            job.fibers.names = names(job.fibers.index(job.fibers.index <= numel(names)));
        else
            job.fibers.names = intersect(job.fibers.names, names);
        end
        if isempty(job.fibers.names)
            error('ftrstruct_Curves2dx:get_fiberNames', 'No matching fiber (names) found.');
        end

        % for each fiberName, create a separate file
        if isfield(job.ftr, 'ftrname')
            [pth nam] = fileparts(job.ftr.ftrname{1});
        else
            pth = pwd;
            nam = '';
        end
        out.dxfiles = {};
        for k = 1:numel(job.fibers.names)
            [Curves errstr] = ftrstruct_query(job.ftr.ftrstruct, 'getCurve', job.fibers.names{k});
            if isempty(errstr)
                out.dxfiles{end+1} = fullfile(pth, sprintf('%s.dx', genvarname(sprintf('%s_%s', nam, job.fibers.names{k}))));
                write_Curves2dx(out.dxfiles{end}, Curves, job.ftr.ftrstruct.hMatrix);
            else
                warning('ftrstruct_Curves2dx:getCurve', errstr)
            end
        end
end

function write_Curves2dx(dxfname, curves, hMatrix, data)
% dxfname - file name
% curves  - cell array with Nx3 voxel coordinates in each cell
% hMatrix - matrix to transform from voxel to edges space - will be used to
%           compute transformation from voxel to NIFTI space
% data    - optional, cell array with Nx1 data values (one per point on
%           line) or 1x1 (one per line) in each cell
if nargin == 3
    data = num2cell(ones(1, numel(curves)));
end
% transformation from index into NIFTI world coordinates - see
% mrstruct_to_nifti for a history of this
% wegen transponieren bei mrstruct
transMx= [[0 1 0 0]; [1 0 0 0]; [0 0 1 0]; [0 0 0 1]];
% re-conversion according to spm_dicom_convert
patient_to_tal = diag([-1 -1 1 1]);
analyze_to_dicom = [diag([1 -1 1]) [0 (hdrStrc.dim(2) + 1) 0]'; 0 0 0 1];  %%% 28.1.2008 BWK Teil von MVM korregiert
corrMy= diag(ones(1, 4)); corrMy(1:3, 4)= -1;        %%% bei neuer norm muss das rein EDGES_NEW
vox2nifti = patient_to_tal*hMatrix*corrMy/transMx*analyze_to_dicom;      %%% bei neuer norm muss das rein EDGES_NEW
fid = fopen(dxfname, 'w');
for k = 1:numel(curves)
    % positions component
    fprintf(fid, '# curve %d\n', k);
    fprintf(fid, 'object %d class array type float rank 1 shape 3 items %d data follows\n', (k-1)*4+2, size(curves{k},1));
    newcurves = vox2nifti*[curves{k} ones(size(curves{k},1),1)];
    str = cellstr(num2str(newcurves(:,1:3)));
    fprintf(fid, '%s\n', str{:});
    fprintf(fid, 'attribute "dep" string "positions"\n#\n');
    % connections component
    fprintf(fid, 'object %d class gridconnections counts %d\n', (k-1)*4+3, size(curves{k},1));
    fprintf(fid, 'attribute "element type" string "lines"\n');
    fprintf(fid, 'attribute "dep" string "connections"\n');
    fprintf(fid, 'attribute "ref" string "positions"\n#\n');
    % data component - can be one item or numel(curves{k})
    if numel(data{k}) == 1
        cl = 'constantarray';
    else
        cl = 'array';
    end
    fprintf(fid, 'object %d class %s type float rank 0 items %d data follows\n', (k-1)*4+4, cl, size(curves{k},1));
    fprintf(fid, '%f ', data{k}(:)');
    fprintf(fid, '\nattribute "dep" string "positions"\n#\n');
    % construct field object
    fprintf(fid, 'object %d class field\n', (k-1)*4+1);
    fprintf(fid, 'component "positions" value %d\n', (k-1)*4+2);
    fprintf(fid, 'component "connections" value %d\n', (k-1)*4+3);
    fprintf(fid, 'component "data" value %d\n', (k-1)*4+4);
    fprintf(fid, 'attribute "name" string "default"\n#\n');
end
% collect fields into group
fprintf(fid, 'object "default" class group\n');
for k = 0:(numel(curves)-1)
    fprintf(fid, 'member %d value %d\n', k, k*4+1);
end
fprintf(fid, '#\nend\n');
fclose(fid);
