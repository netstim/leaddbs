% interface program for batch program
%
% Author: Susanne Schnell
% PC 25.07.2008


function out = dti_tracking_mori_ui(P)
%load DTD
[path,nam,ext] = fileparts(P.filename{1});
dtdStruct = dtdstruct_read(P.filename{1});

% create or load start mask
if isfield(P.start,'startdef')
    %load mask as mrstruct
    [startMask,errStr] = maskstruct_read(P.start.startdef.startfile{1});
    if isequal(size(dtdStruct.b0_image_struc.dataAy), maskstruct_query(startMask, 'sizeAy'))
        if isfield(P.start.startdef.startmask,'startnumber')
            startMask = mrstruct_init('volume',...
                double(maskstruct_query(startMask,'getMaskVc',P.start.startdef.startmask.startnumber)),dtdStruct.b0_image_struc);
        elseif isfield(P.start.startdef.startmask,'startname')
            startMask = mrstruct_init('volume',...
                double(maskstruct_query(startMask,'getMask',P.start.startdef.startmask.startname)),dtdStruct.b0_image_struc);
        end
        startMask.user.upperLim_FA= 'undef';   startMask.user.lowerLim_FA= 'undef';
        startMask.user.upperLim_TrD= 'undef';   startMask.user.lowerLim_TrD= 'undef';
    else
        error('error: invalid roidata of startmask');
    end
    
    if ~isempty(errStr)
        error('Selected file for stop mask is not of type "mrstruct".')
    end
else    
    mrStruct= dtdstruct_query(dtdStruct, 'getFA');
    startMask= mrstruct_init('volume', ...
        double(mrStruct.dataAy > P.start.startthreshold.startfaLim),mrStruct);
    mrStruct= dtdstruct_query(dtdStruct, 'getTrace');
    startMask.dataAy= double(startMask.dataAy & (mrStruct.dataAy < P.start.startthreshold.startTrLim));
    startMask.user.upperLim_FA= [];
    startMask.user.lowerLim_FA= [P.start.startthreshold.startfaLim];
    startMask.user.upperLim_TrD= [P.start.startthreshold.startTrLim];
    startMask.user.lowerLim_TrD= [];
end

% create or load stop mask
if isfield(P.stop,'stopdef')
    %load mask as mrstruct
    [stopMask,errStr] = maskstruct_read(P.stop.stopdef.stopfile{1});
    if isequal(size(dtdStruct.b0_image_struc.dataAy), maskstruct_query(stopMask, 'sizeAy'))
        if isfield(P.stop.stopdef.stopmask,'stopnumber')
            stopMask = mrstruct_init('volume',...
                double(maskstruct_query(stopMask,'getMaskVc',P.stop.stopdef.stopmask.stopnumber)),dtdStruct.b0_image_struc);
        elseif isfield(P.stop.stopdef.stopmask,'stopname')
            stopMask = mrstruct_init('volume',...
                double(maskstruct_query(stopMask,'getMask',P.stop.stopdef.stopmask.stopname)),dtdStruct.b0_image_struc);
        end
        stopMask.user.upperLim_FA= 'undef';
        stopMask.user.lowerLim_FA= 'undef';
        stopMask.user.upperLim_TrD= 'undef';
        stopMask.user.lowerLim_TrD= 'undef';
    else
        error('error: invalid roidata of startmask');
    end
    
    if ~isempty(errStr)
        error('Selected file for stop mask is not of type "mrstruct".')
    end
else
    mrStruct= dtdstruct_query(dtdStruct, 'getFA');
    stopMask= mrstruct_init('volume', double(mrStruct.dataAy > P.stop.stopthreshold.stopfaLim), mrStruct);
    mrStruct= dtdstruct_query(dtdStruct, 'getTrace');
    stopMask.dataAy= double(stopMask.dataAy & (mrStruct.dataAy < P.stop.stopthreshold.stopTrLim));
    stopMask.user.upperLim_FA= [];
    stopMask.user.lowerLim_FA= [P.stop.stopthreshold.stopfaLim];
    stopMask.user.upperLim_TrD= [P.stop.stopthreshold.stopTrLim];
    stopMask.user.lowerLim_TrD= [];
end

% random sampling?
if P.randsampler == 1
    [ftrStruct, errStr]= ftrack_mori(dtdStruct, startMask, stopMask, abs(cos(pi*P.maxangle/180)), P.minvox, P.randsampno);
else
    [ftrStruct, errStr]= ftrack_mori(dtdStruct, startMask, stopMask, abs(cos(pi*P.maxangle/180)), P.minvox, []);
end

% save file
if isempty(errStr)
    [path,nam,ext] = fileparts(P.filename{1});
    if isfield(P.newmorifile,'auto')
        if length(nam) > 4 && strcmp(nam(end-3:end), '_DTD')
            nam = nam(1:end-4);
        end
        out.files{1} = fullfile(path, [nam,'_FTR.mat']);
    else
        out.files{1} = fullfile(P.newmorifile.out.dir{1}, P.newmorifile.out.fname);
    end
    ftrstruct_write(ftrStruct,out.files{1});
else
    error(errStr);
end
