function fid = ea_native2acpc(cfg,leadfig,automan)
% leadfig can be either a handle to the lead figure with selected patient
% directories or a cell of patient directories.
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% fidpoints_mm = [-0.4,1.53,-2.553       % AC
%     -0.24,-26.314,-4.393            % PC
%     -0.4,1.53,20];              % Midsag

if iscell(leadfig)
    uidir = leadfig;
else
    uidir = getappdata(leadfig,'uipatdir');
end

if ~exist('automan','var')
    automan = 'auto';
end

fidpoints_mm = [0.25,1.298,-5.003       % AC
    -0.188,-24.756,-2.376            % PC
    0.25,1.298,55];              % Midsag

native = [cfg.xmm,cfg.ymm,cfg.zmm];

if isempty(uidir)
    ea_error('Please choose and normalize patients first.');
end

% ea_dispercent(0, 'Iterating through patients');
for pt = 1:length(uidir)
    % ea_dispercent(pt/length(uidir));
    options = ea_getptopts(uidir{pt});
    json = loadjson(options.subj.norm.log.method);
    if contains(json.method, 'SPM')
        normTemplate = [ea_space, 'TPM.nii'];
    else
        spacedef = ea_getspacedef;
        normTemplate = [ea_space, spacedef.templates{1}, '.nii'];
    end

    switch automan

        case 'auto'
            fidpoints_vox = ea_getfidpoints(fidpoints_mm,normTemplate);

            % warp into patient space
            [fpinsub_mm] = ea_map_coords(fidpoints_vox', normTemplate, [options.subj.subjDir,filesep,'forwardTransform'], '');
            fpinsub_mm = fpinsub_mm';

            fid(pt).AC = fpinsub_mm(1,:);
            fid(pt).PC = fpinsub_mm(2,:);
            fid(pt).MSP = fpinsub_mm(3,:);

        case 'manual'
            copyfile([options.subj.subjDir,filesep,'F.fcsv'],[options.subj.subjDir,filesep,'F.dat'])
            Ct = readtable([options.subj.subjDir,filesep,'F.dat']);

            % AC
            cnt = 1;
            fid(pt).AC = zeros(1,3);
            for i = 17:19
                thisval = table2array(Ct(i,1));
                fid(pt).AC(cnt) = str2double(thisval{1});
                cnt = cnt+1;
            end

            % PC
            cnt = 1;
            fid(pt).PC = zeros(1,3);
            for i = 31:33
                thisval = table2array(Ct(i,1));
                fid(pt).PC(cnt) = str2double(thisval{1});
                cnt = cnt+1;
            end

            % MSP
            cnt = 1;
            fid(pt).MSP = zeros(1,3);
            for i = 45:47
                thisval = table2array(Ct(i,1));
                fid(pt).MSP(cnt) = str2double(thisval{1});
                cnt = cnt+1;
            end

            fpinsub_mm(1,:) = fid(pt).AC;
            fpinsub_mm(2,:) = fid(pt).PC;
            fpinsub_mm(3,:) = fid(pt).MSP;
    end

    % x-dimension
    A = fpinsub_mm(3,:)-fpinsub_mm(1,:);
    B = fpinsub_mm(2,:)-fpinsub_mm(1,:);
    xvec = cross(A,B); % normal to given plane

    xvec = xvec/norm(xvec);
    % y-dimension (just move from ac to pc and scale by y dimension):
    yvec = (fpinsub_mm(2,:)-fpinsub_mm(1,:));
    yvec = yvec/norm(yvec);
    yvec = -yvec;

    zvec = cross(xvec,yvec);
    zvec = zvec/norm(zvec);

    switch cfg.acmcpc
        case 1 % relative to AC:
            warpcoord_mm = mldivide([xvec',yvec',zvec'],native'-fpinsub_mm(1,:)');
        case 2 % relative to midcommissural point:
            warpcoord_mm = mldivide([xvec',yvec',zvec'],native'-mean([fpinsub_mm(1,:);fpinsub_mm(2,:)],1)');
        case 3 % relative to PC:
            warpcoord_mm = mldivide([xvec',yvec',zvec'],native'-fpinsub_mm(2,:)');
    end

    fid(pt).WarpedPointACPC = [warpcoord_mm(1),warpcoord_mm(2),warpcoord_mm(3)];
    fid(pt).WarpedPointNative = native;
end
% ea_dispercent(1,'end');

assignin('base','fid',fid);


function fidpoints_vox = ea_getfidpoints(fidpoints_mm,tempfile)

V = spm_vol(tempfile);
fidpoints_vox = V(1).mat\[fidpoints_mm,ones(size(fidpoints_mm,1),1)]';
fidpoints_vox = fidpoints_vox(1:3,:)';
