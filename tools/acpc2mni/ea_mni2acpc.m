function fid=ea_mni2acpc(cfg,leadfig)
% leadfig can be either a handle to the lead figure with selected patient
% directories or a cell of patient directories.
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% fidpoints_mm=[-0.4,1.53,-2.553       % AC
%     -0.24,-26.314,-4.393            % PC
%     -0.4,1.53,20];              % Midsag
fidpoints_mm=[ 0.25,   1.298,   -5.003    % AC
              -0.188, -24.756,  -2.376    % PC
               0.25,   1.298,    55];     % Midsag

if iscell(leadfig)
    uidir=leadfig;
else
    uidir=getappdata(leadfig,'uipatdir');
end

% prompt for MNI-coordinates:


mni=[cfg.xmm,cfg.ymm,cfg.zmm];
fidpoints_mm=[fidpoints_mm;mni];

if isempty(uidir)
    ea_error('Please choose and normalize patients first.');
end

disp('*** Converting MNI-coordinates to ACPC based on normalizations in selected patients.');
%ea_dispercent(0,'Iterating through patients');
for pt=1:length(uidir)
 %   ea_dispercent(pt/length(uidir));
    directory=[uidir{pt},filesep];
    options = ea_getptopts(directory);

    [~,template]=ea_whichnormmethod(directory);

    fidpoints_vox=ea_getfidpoints(fidpoints_mm,template);

    % warp into patient space:

    [fpinsub_mm] = ea_map_coords(fidpoints_vox', template, [options.subj.subjDir,filesep,'forwardTransform'], '');
    fpinsub_mm=fpinsub_mm';


    fid(pt).AC=fpinsub_mm(1,:);
    fid(pt).PC=fpinsub_mm(2,:);
    fid(pt).MSP=fpinsub_mm(3,:);

    % x-dimension
    A=fpinsub_mm(3,:)-fpinsub_mm(1,:);
    B=fpinsub_mm(2,:)-fpinsub_mm(1,:);
    xvec=cross(A,B); %normal to given plane
    xvec=xvec/norm(xvec);
    % y-dimension (just move from ac to pc and scale by y dimension):
    yvec=(fpinsub_mm(2,:)-fpinsub_mm(1,:));
    yvec=yvec/norm(yvec);
    yvec=-yvec; % this vector points to anterior of the AC!! (but acpc y coordinates are given with - sign).

    zvec=cross(xvec,yvec);
    zvec=zvec/norm(zvec);

    switch cfg.acmcpc
        case 1 % relative to AC:
            warpcoord_mm=mldivide([xvec',yvec',zvec'],fpinsub_mm(4,:)'-fpinsub_mm(1,:)');
        case 2 % relative to midcommissural point:
            warpcoord_mm=mldivide([xvec',yvec',zvec'],fpinsub_mm(4,:)'-mean([fpinsub_mm(1,:);fpinsub_mm(2,:)],1)');
        case 3 % relative to PC:
            warpcoord_mm=mldivide([xvec',yvec',zvec'],fpinsub_mm(4,:)'-fpinsub_mm(2,:)');
    end

    fid(pt).WarpedPointACPC=[warpcoord_mm(1),warpcoord_mm(2),warpcoord_mm(3)];
    fid(pt).WarpedPointNative=fpinsub_mm(4,:);
end
%ea_dispercent(1,'end');

assignin('base','fid',fid);





function fidpoints_vox=ea_getfidpoints(fidpoints_mm,tempfile)

V=spm_vol(tempfile);
fidpoints_vox=V(1).mat\[fidpoints_mm,ones(size(fidpoints_mm,1),1)]';
fidpoints_vox=fidpoints_vox(1:3,:)';





function o=cell2acpc(acpc)

acpc=ea_strsplit(acpc{1},' ');
if length(acpc)~=3
    acpc=ea_strsplit(acpc{1},',');
    if length(acpc)~=3
        ea_error('Please enter 3 values separated by spaces or commas.');
    end
end
for dim=1:3
o(dim,1)=str2double(acpc{dim});
end
