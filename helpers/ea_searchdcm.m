function dcm = ea_searchdcm(src)
% Walk through DICOM folder to find DCMs to be converted
%
% Return a N*1 cell, N is the number of detected series. Each element of
% the cell is a cell of DCM file paths to be converted to NIfTI.
%
% Adopted from dicm2nii by Xiangrui Li
% See: https://github.com/xiangruili/dicm2nii


%% Retrieve all files in the folder
fnames = ea_regexpdir(src, '.*', 1, 'file');

nFile = numel(fnames);
if nFile<1
    error(' No files found in the data source.');
end

%% Check each file, store partial header in cell array hh
% first 3 fields are must. First 10 indexed in code
flds = {'Columns' 'Rows' 'BitsAllocated' 'SeriesInstanceUID' 'SeriesNumber' ...
    'ImageOrientationPatient' 'ImagePositionPatient' 'PixelSpacing' ...
    'SliceThickness' 'SpacingBetweenSlices' ... % these 10 indexed in code
    'PixelRepresentation' 'BitsStored' 'HighBit' 'SamplesPerPixel' ...
    'PlanarConfiguration' 'EchoTime' 'RescaleIntercept' 'RescaleSlope' ...
    'InstanceNumber' 'NumberOfFrames' 'B_value' 'DiffusionGradientDirection' ...
    'RTIA_timer' 'RBMoCoTrans' 'RBMoCoRot' 'AcquisitionNumber'};
dict = dicm_dict('SIEMENS', flds); % dicm_hdr will update vendor if needed

% read header for all files, use parpool if available and worthy
fprintf('Validating %g files ...\n', nFile);
hh = cell(1, nFile); errStr = cell(1, nFile);
useParTool = 0;
doParFor = nFile>2000 && useParTool;
for k = 1:nFile
    [hh{k}, errStr{k}, dict] = dicm_hdr(fnames{k}, dict);
    if doParFor && ~isempty(hh{k}) % parfor wont allow updating dict
        parfor i = k+1:nFile
            [hh{i}, errStr{i}] = dicm_hdr(fnames{i}, dict);
        end
        break;
    end
end

%% sort headers into cell h by SeriesInstanceUID, EchoTime and InstanceNumber
h = {}; % in case of no dicom files at all
errInfo = '';
seriesUIDs = {}; ETs = {};
for k = 1:nFile
    s = hh{k};
    if isempty(s) || any(~isfield(s, flds(1:3))) || ~isfield(s, 'PixelData') ...
            || (isstruct(s.PixelData) && s.PixelData.Bytes<1)
        if ~isempty(errStr{k}) % && ~contains(errInfo, errStr{k})
            errInfo = sprintf('%s\n%s\n', errInfo, errStr{k});
        end
        continue; % skip the file
    end

    if isfield(s, flds{4})
        sUID = s.SeriesInstanceUID;
    else
        if isfield(s, 'SeriesNumber'), sN = s.SeriesNumber;
        else, sN = fix(toc*1e6);
        end
        sUID = num2str(sN); % make up UID
        if isfield(s, 'SeriesDescription')
            sUID = [s.SeriesDescription sUID];
        end
    end

    m = find(strcmp(sUID, seriesUIDs));
    if isempty(m)
        m = numel(seriesUIDs)+1;
        seriesUIDs{m} = sUID;
        ETs{m} = [];
    end

    % EchoTime is needed for Siemens fieldmap mag series
    et = tryGetField(s, 'EchoTime');
    if isempty(et), i = 1;
    else
        i = find(et == ETs{m}); % strict equal?
        if isempty(i)
            i = numel(ETs{m}) + 1;
            ETs{m}(i) = et;
            if i>1
                [ETs{m}, ind] = sort(ETs{m});
                i = find(et == ETs{m});
                h{m}{end+1}{1} = [];
                h{m} = h{m}(ind);
            end
        end
    end
    j = tryGetField(s, 'InstanceNumber');
    if isempty(j) || j<1
        try j = numel(h{m}{i}) + 1;
        catch, j = 1;
        end
    end
    h{m}{i}{j} = s; % sort partial header
end
clear hh errStr;

%% Check headers: remove dim-inconsistent series
nRun = numel(h);
if nRun<1 % no valid series
    errorLog(sprintf('No valid files found:\n%s.', errInfo));
    return;
end
keep = true(1, nRun); % true for useful runs
subjs = cell(1, nRun); vendor = cell(1, nRun);
sNs = ones(1, nRun); studyIDs = cell(1, nRun);
fldsCk = {'ImageOrientationPatient' 'NumberOfFrames' 'Columns' 'Rows' ...
          'PixelSpacing' 'RescaleIntercept' 'RescaleSlope' 'SamplesPerPixel' ...
          'SpacingBetweenSlices' 'SliceThickness'}; % last for thickness
for i = 1:nRun
    h{i} = [h{i}{:}]; % concatenate different EchoTime
    ind = cellfun(@isempty, h{i});
    h{i}(ind) = []; % remove all empty cell for all vendors

    s = h{i}{1};
    if ~isfield(s, 'LastFile') % avoid re-read for PAR/HEAD/BV file
        s = dicm_hdr(s.Filename); % full header for 1st file
    end
    if ~isfield(s, 'Manufacturer'), s.Manufacturer = 'Unknown'; end
    subjs{i} = PatientName(s);
    vendor{i} = s.Manufacturer;
    if isfield(s, 'SeriesNumber'), sNs(i) = s.SeriesNumber;
    else, sNs(i) = fix(toc*1e6);
    end
    studyIDs{i} = tryGetField(s, 'StudyID', '1');
    series = sprintf('Subject %s, %s (Series %g)', subjs{i}, ProtocolName(s), sNs(i));
    s = multiFrameFields(s); % no-op if non multi-frame
    if isempty(s), keep(i) = 0; continue; end % invalid multiframe series
    s.isDTI = isDTI(s);
    if ~isfield(s, 'AcquisitionDateTime') % assumption: 1st instance is earliest
        try s.AcquisitionDateTime = [s.AcquisitionDate s.AcquisitionTime]; end
    end

    h{i}{1} = s; % update record in case of full hdr or multiframe

    nFile = numel(h{i});
    if nFile>1 && tryGetField(s, 'NumberOfFrames', 1) > 1 % seen in vida
        for k = 2:nFile % this can be slow
            h{i}{k} = dicm_hdr(h{i}{k}.Filename); % full header
            h{i}{k} = multiFrameFields(h{i}{k});
        end
        if ~isfield(s, 'EchoTimes') && isfield(s, 'EchoTime')
            h{i}{1}.EchoTimes = nan(1, nFile);
            for k = 1:nFile, h{i}{1}.EchoTimes(k) = h{i}{k}.EchoTime; end
        end
    end

    % check consistency in 'fldsCk'
    nFlds = numel(fldsCk);
    if isfield(s, 'SpacingBetweenSlices'), nFlds = nFlds - 1; end % check 1 of 2
    for k = 1:nFlds*(nFile>1)
        if isfield(s, fldsCk{k}), val = s.(fldsCk{k}); else, continue; end
        val = repmat(double(val), [1 nFile]);
        for j = 2:nFile
            if isfield(h{i}{j}, fldsCk{k}), val(:,j) = h{i}{j}.(fldsCk{k});
            else, keep(i) = 0; break;
            end
        end
        if ~keep(i), break; end % skip silently
        ind = any(abs(bsxfun(@minus, val, val(:,1))) > 1e-4, 1);
        if sum(ind)>1 % try 2nd, in case only 1st is inconsistent
            ind = any(abs(bsxfun(@minus, val, val(:,2))) > 1e-4, 1);
        end
        if ~any(ind), continue; end % good
        if any(strcmp(fldsCk{k}, {'RescaleIntercept' 'RescaleSlope'}))
            h{i}{1}.ApplyRescale = true;
            continue;
        end
        if numel(ind)>2 && sum(ind)==1 % 2+ files but only 1 inconsistent
            h{i}(ind) = []; % remove first or last, but keep the series
            nFile = nFile - 1;
            if ind(1) % re-do full header for new 1st file
                s = dicm_hdr(h{i}{1}.Filename);
                s.isDTI = isDTI(s);
                h{i}{1} = s;
            end
        else
            errorLog(['Inconsistent ''' fldsCk{k} ''' for ' series '. Series skipped.']);
            keep(i) = 0; break;
        end
    end

    nSL = nMosaic(s); % nSL>1 for mosaic
    if ~isempty(nSL) && nSL>1
        h{i}{1}.isMos = true;
        h{i}{1}.LocationsInAcquisition = nSL;
        if s.isDTI, continue; end % allow missing directions for DTI
        a = zeros(1, nFile);
        for j = 1:nFile, a(j) = tryGetField(h{i}{j}, 'InstanceNumber', 1); end
        if any(diff(a) ~= 1) % like CMRR ISSS seq or multi echo
            errorLog(['InstanceNumber discontinuity detected for ' series '.' ...
                'See VolumeTiming in NIfTI ext or dcmHeaders.mat.']);
            dict = dicm_dict('', { 'AcquisitionDate' 'AcquisitionTime'});
            vTime = nan(1, nFile);
            for j = 1:nFile
                s2 = dicm_hdr(h{i}{j}.Filename, dict);
                dt = [s2.AcquisitionDate s2.AcquisitionTime];
                vTime(j) = datenum(dt, 'yyyymmddHHMMSS.fff');
            end
            vTime = vTime - min(vTime);
            h{i}{1}.VolumeTiming = vTime*24*3600; % day to seconds
        end
        continue; % no other check for mosaic
    end

    if ~keep(i) || nFile<2 || ~isfield(s, 'ImagePositionPatient'), continue; end
    if tryGetField(s, 'NumberOfFrames', 1) > 1, continue; end % Siemens Vida

    ipp = zeros(nFile, 1);
    iSL = xform_mat(s); iSL = iSL(3);
    for j = 1:nFile, ipp(j,:) = h{i}{j}.ImagePositionPatient(iSL); end
    gantryTilt = abs(tryGetField(s, 'GantryDetectorTilt', 0)) > 0.1;
    [err, nSL, sliceN, isTZ] = checkImagePosition(ipp, gantryTilt);
    if ~isempty(err)
        errorLog([err ' for ' series '. Series skipped.']);
        keep(i) = 0; continue; % skip
    end
    h{i}{1}.LocationsInAcquisition = uint16(nSL); % best way for nSL?

    nVol = nFile / nSL;
    if isTZ % Philips
        ind = reshape(1:nFile, [nVol nSL])';
        h{i} = h{i}(ind(:));
    end

    % re-order slices within vol. No SliceNumber since files are organized
    if all(diff(sliceN, 2) == 0), continue; end % either 1:nSL or nSL:-1:1
    if sliceN(end) == 1, sliceN = sliceN(nSL:-1:1); end % not important
    inc = repmat((0:nVol-1)*nSL, nSL, 1);
    ind = repmat(sliceN(:), nVol, 1) + inc(:);
    h{i} = h{i}(ind); % sorted by slice locations

    if sliceN(1) == 1, continue; end % first file kept: following update h{i}{1}
    h{i}{1} = dicm_hdr(h{i}{1}.Filename); % read full hdr
    s = h{i}{sliceN==1}; % original first file
    fldsCp = {'AcquisitionDateTime' 'isDTI' 'LocationsInAcquisition'};
    for j = 1:numel(fldsCp)
        if isfield(h{i}{1}, fldsCk{k}), h{i}{1}.(fldsCp{j}) = s.(fldsCp{j}); end
    end
end
h = h(keep);

dcm = cell(length(h), 1);
for i=1:length(h)
    dcm{i} = cellfun(@(X) X.Filename, h{i}, 'Uni', 0)';
end


%% subfuction: check whether parpool is available
% Return true if it is already open, or open it if available
function doParal = useParTool
doParal = usejava('jvm');
if ~doParal, return; end

if isempty(which('parpool')) % for early matlab versions
    try
        if matlabpool('size')<1 %#ok<*DPOOL>
            try
                matlabpool;
            catch me
                fprintf(2, '%s\n', me.message);
                doParal = false;
            end
        end
    catch
        doParal = false;
    end
    return;
end

% Following for later matlab with parpool
try
    if isempty(gcp('nocreate'))
        try
            parpool;
        catch me
            fprintf(2, '%s\n', me.message);
            doParal = false;
        end
    end
catch
    doParal = false;
end


%% Subfunction: get field if exist, return default value otherwise
function val = tryGetField(s, field, dftVal)
if isfield(s, field), val = s.(field);
elseif nargin>2, val = dftVal;
else, val = [];
end


%% Subfunction: return PatientName
function subj = PatientName(s)
subj = tryGetField(s, 'PatientName');
if isempty(subj), subj = tryGetField(s, 'PatientID', 'Anonymous'); end


%% Subfunction: return SeriesDescription
function name = ProtocolName(s)
name = tryGetField(s, 'SeriesDescription');
if isempty(name) || (strncmp(s.Manufacturer, 'SIEMENS', 7) && any(regexp(name, 'MoCoSeries$')))
    name = tryGetField(s, 'ProtocolName');
end
if isempty(name), [~, name] = fileparts(s.Filename); end


%% subfunction: extract useful fields for multiframe dicom
function s = multiFrameFields(s)
pffgs = 'PerFrameFunctionalGroupsSequence';
sfgs = 'SharedFunctionalGroupsSequence';
if any(~isfield(s, {sfgs pffgs})), return; end
try nFrame = s.NumberOfFrames; catch, nFrame = numel(s.(pffgs).FrameStart); end

% check slice ordering (Philips often needs SortFrames)
n = numel(MF_val('DimensionIndexValues', s, 1));
s2 = struct('DimensionIndexValues', nan(n, nFrame), 'B_value', zeros(1, nFrame));
s2 = dicm_hdr(s, s2, 1:nFrame); a = s2.DimensionIndexValues';
[ind, nSL] = sort_frames([a(:,2) s2.B_value'], a(:, [3:end 1]));
if ~isequal(ind, 1:nFrame)
    if ind(1) ~= 1 || ind(end) ~= nFrame
        s = dicm_hdr(s.Filename, [], ind([1 end])); % re-read new frames [1 end]
    end
    s.SortFrames = ind; % will use to sort img and get iVol/iSL for PerFrameSQ
end
if ~isfield(s, 'LocationsInAcquisition'), s.LocationsInAcquisition = nSL; end

% copy important fields into s
flds = {'EchoTime' 'PixelSpacing' 'SpacingBetweenSlices' 'SliceThickness' ...
        'RepetitionTime' 'FlipAngle' 'RescaleIntercept' 'RescaleSlope' ...
        'ImageOrientationPatient' 'ImagePositionPatient' ...
        'InPlanePhaseEncodingDirection' 'MRScaleSlope' 'CardiacTriggerDelayTime'};
iF = 1; if isfield(s, 'SortFrames'), iF = s.SortFrames(1); end
for i = 1:numel(flds)
    if isfield(s, flds{i}), continue; end
    a = MF_val(flds{i}, s, iF);
    if ~isempty(a), s.(flds{i}) = a; end
end

if ~isfield(s, 'EchoTime')
    a = MF_val('EffectiveEchoTime', s, iF);
    if ~isempty(a), s.EchoTime = a;
    else, try s.EchoTime = str2double(s.EchoTimeDisplay); end
    end
end

% for Siemens: the redundant copy makes non-Siemens code faster
if isfield(s.(sfgs).Item_1, 'CSASeriesHeaderInfo')
    s.CSASeriesHeaderInfo = s.(sfgs).Item_1.CSASeriesHeaderInfo.Item_1;
end
fld = 'CSAImageHeaderInfo';
if isfield(s.(pffgs).Item_1, fld)
    s.(fld) = s.(pffgs).(sprintf('Item_%g', iF)).(fld).Item_1;
end

% check ImageOrientationPatient consistency for 1st and last frame only
iF = nFrame; if isfield(s, 'SortFrames'), iF = s.SortFrames(iF); end
a = MF_val('ImagePositionPatient', s, iF);
if ~isempty(a), s.LastFile.ImagePositionPatient = a; end
fld = 'ImageOrientationPatient';
val = MF_val(fld, s, iF);
if ~isempty(val) && isfield(s, fld) && any(abs(val-s.(fld))>1e-4)
    s = []; return; % inconsistent orientation, skip
end


%% Subfunction: return true if series is DTI
function tf = isDTI(s)
tf = isType(s, '\DIFFUSION'); % Siemens, Philips
if tf, return; end
if isfield(s, 'ProtocolDataBlock') % GE, not labeled as \DIFFISION
    IOPT = tryGetField(s.ProtocolDataBlock, 'IOPT');
    if isempty(IOPT), tf = tryGetField(s, 'DiffusionDirection', 0)>0;
    else, tf = ~isempty(regexp(IOPT, 'DIFF', 'once'));
    end
elseif strncmpi(s.Manufacturer, 'Philips', 7)
    tf = strcmp(tryGetField(s, 'MRSeriesDiffusion', 'N'), 'Y');
else % Some Siemens DTI are not labeled as \DIFFUSION
    tf = ~isempty(csa_header(s, 'B_value'));
end


%% Subfunction: return true if keyword is in s.ImageType
function tf = isType(s, keyword)
typ = tryGetField(s, 'ImageType', '');
tf = contains(typ, keyword); %#ok<*STREMP>


%% Write error info to a file in case user ignores Command Window output
function errorLog(errInfo)
fprintf(2, ' %s\n', errInfo); % red text in Command Window


%% Subfunction: return NumberOfImagesInMosaic if Siemens mosaic, or [] otherwise.
% If NumberOfImagesInMosaic in CSA is >1, it is mosaic, and we are done.
% If not exists, it may still be mosaic due to Siemens bug seen in syngo MR
% 2004A 4VA25A phase image. Then we check EchoColumnPosition in CSA, and if it
% is smaller than half of the slice dim, sSliceArray.lSize is used as nMos. If
% no CSA at all, the better way may be to peek into img to get nMos. Then the
% first attempt is to check whether there are padded zeros. If so we count zeros
% either at top or bottom of the img to decide real slice dim. In case there is
% no padded zeros, we use the single zero lines along row or col seen in most
% (not all, for example some phase img, derived data like moco series or tmap
% etc) mosaic. If the lines are equally spaced, and nMos is divisible by mosaic
% dim, we accept nMos. Otherwise, we fall back to NumberOfPhaseEncodingSteps,
% which is used by dcm2nii, but is not reliable for most mosaic due to partial
% fourier or less 100% phase fov.
function nMos = nMosaic(s)
nMos = csa_header(s, 'NumberOfImagesInMosaic'); % healthy mosaic dicom
if ~isempty(nMos), return; end % seen 0 for GLM Design file and others

% The next fix detects mosaic which is not labeled as MOSAIC in ImageType, nor
% NumberOfImagesInMosaic exists, seen in syngo MR 2004A 4VA25A phase image.
res = csa_header(s, 'EchoColumnPosition'); % half or full of slice dim
if ~isempty(res)
    dim = max([s.Columns s.Rows]);
    interp = asc_header(s, 'sKSpace.uc2DInterpolation');
    if ~isempty(interp) && interp, dim = dim / 2; end
    if dim/res/2 >= 2 % nTiles>=2
        nMos = asc_header(s, 'sSliceArray.lSize'); % mprage lSize=1
    end
    return; % Siemens non-mosaic returns here
end

% The fix below is for dicom labeled as \MOSAIC in ImageType, but no CSA.
if ~isType(s, '\MOSAIC'), return; end % non-Siemens returns here
try nMos = s.LocationsInAcquisition; return; end % try Siemens private tag

dim = double([s.Columns s.Rows]); % slice or mosaic dim
img = dicm_img(s, 0) ~= 0; % peek into img to figure out nMos
nP = tryGetField(s, 'NumberOfPhaseEncodingSteps', 4); % sliceDim >= phase steps
c = img(dim(1)-nP:end, dim(2)-nP:end); % corner at bottom-right
done = false;
if all(~c(:)) % at least 1 padded slice: not 100% safe
    c = img(1:nP+1, dim(2)-nP:end); % top-right
    if all(~c(:)) % all right tiles padded: use all to determine
        ln = sum(img);
    else % use several rows at bottom to determine: not as safe as all
        ln = sum(img(dim(1)-nP:end, :));
    end
    z = find(ln~=0, 1, 'last');
    nMos = dim(2) / (dim(2) - z);
    done = mod(nMos,1)==0 && mod(dim(1),nMos)==0;
end
if ~done % this relies on zeros along row or col seen in most mosaic
    ln = sum(img, 2) == 0;
    if sum(ln)<2
        ln = sum(img) == 0; % likely PhaseEncodingDirectionPositive=0
        i = find(~ln, 1, 'last'); % last non-zero column in img
        ln(i+2:end) = []; % leave only 1 true for padded zeros
    end
    nMos = sum(ln);
    done = nMos>1 && all(mod(dim,nMos)==0) && all(diff(find(ln),2)==0);
end
if ~done && isfield(s, 'NumberOfPhaseEncodingSteps')
    nMos = min(dim) / nP;
    done = nMos>1 && mod(nMos,1)==0 && all(mod(dim,nMos)==0);
end

if ~done
    errorLog([ProtocolName(s) ': NumberOfImagesInMosaic not available.']);
    nMos = []; % keep mosaic as it is
    return;
end

nMos = nMos * nMos;
img = mos2vol(uint8(img), nMos); % find padded slices: useful for STC
while 1
    a = img(:,:,nMos);
    if any(a(:)), break; end
    nMos = nMos - 1;
end


%% Subfunction, return a parameter from CSA Image/Series header
function val = csa_header(s, key)
val = [];
fld = 'CSAImageHeaderInfo';
if isfield(s, fld) && isfield(s.(fld), key), val = s.(fld).(key); return; end
fld = 'CSASeriesHeaderInfo';
if isfield(s, fld) && isfield(s.(fld), key), val = s.(fld).(key); return; end


%% Subfunction: get dicom xform matrix and related info
function [ixyz, R, pixdim, xyz_unit] = xform_mat(s, dim)
haveIOP = isfield(s, 'ImageOrientationPatient');
if haveIOP, R = reshape(s.ImageOrientationPatient, 3, 2);
else, R = [1 0 0; 0 1 0]';
end
R(:,3) = cross(R(:,1), R(:,2)); % right handed, but sign may be wrong
foo = abs(R);
[~, ixyz] = max(foo); % orientation info: perm of 1:3
if ixyz(2) == ixyz(1), foo(ixyz(2),2) = 0; [~, ixyz(2)] = max(foo(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
if nargout<2, return; end
iSL = ixyz(3); % 1/2/3 for Sag/Cor/Tra slice
signSL = sign(R(iSL, 3));

try
    pixdim = s.PixelSpacing;
    xyz_unit = 2; % mm
catch
    pixdim = [1 1]'; % fake
    xyz_unit = 0; % no unit information
end
thk = tryGetField(s, 'SpacingBetweenSlices');
if isempty(thk), thk = tryGetField(s, 'SliceThickness', pixdim(1)); end
pixdim = [pixdim; thk];
haveIPP = isfield(s, 'ImagePositionPatient');
if haveIPP, ipp = s.ImagePositionPatient; else, ipp = -(dim'.* pixdim)/2; end
% Next is almost dicom xform matrix, except mosaic trans and unsure slice_dir
R = [R * diag(pixdim) ipp];

% rest are former: R = verify_slice_dir(R, s, dim, iSL)
if dim(3)<2, return; end % don't care direction for single slice

if s.Columns > dim(1) % Siemens mosaic: use dim(1) since no transpose to img
    R(:,4) = R * [ceil(sqrt(dim(3))-1)*dim(1:2)/2 0 1]'; % real slice location
    vec = csa_header(s, 'SliceNormalVector'); % mosaic has this
    if ~isempty(vec) % exist for all tested data
        if sign(vec(iSL)) ~= signSL, R(:,3) = -R(:,3); end
        return;
    end
elseif isfield(s, 'LastFile') && isfield(s.LastFile, 'ImagePositionPatient')
    R(:, 3) = (s.LastFile.ImagePositionPatient - R(:,4)) / (dim(3)-1);
    thk = sqrt(sum(R(:,3).^2)); % override slice thickness if it is off
    if abs(pixdim(3)-thk)/thk > 0.01, pixdim(3) = thk; end
    return; % almost all non-mosaic images return from here
end

% Rest of the code is almost unreachable
if isfield(s, 'CSASeriesHeaderInfo') % Siemens both mosaic and regular
    ori = {'Sag' 'Cor' 'Tra'}; ori = ori{iSL};
    sNormal = asc_header(s, ['sSliceArray.asSlice[0].sNormal.d' ori]);
    if asc_header(s, ['sSliceArray.ucImageNumb' ori]), sNormal = -sNormal; end
    if sign(sNormal) ~= signSL, R(:,3) = -R(:,3); end
    if ~isempty(sNormal), return; end
end

pos = []; % volume center we try to retrieve
if isfield(s, 'LastScanLoc') && isfield(s, 'FirstScanLocation') % GE
    pos = (s.LastScanLoc + s.FirstScanLocation) / 2; % mid-slice center
    if iSL<3, pos = -pos; end % RAS convention!
    pos = pos - R(iSL, 1:2) * (dim(1:2)'-1)/2; % mid-slice location
end

if isempty(pos) && isfield(s, 'Stack') % Philips
    ori = {'RL' 'AP' 'FH'}; ori = ori{iSL};
    pos = tryGetField(s.Stack.Item_1, ['MRStackOffcentre' ori]);
    pos = pos - R(iSL, 1:2) * dim(1:2)'/2; % mid-slice location
end

if isempty(pos) % keep right-handed, and warn user
    if haveIPP && haveIOP
        errorLog(['Please check whether slices are flipped: ' s.NiftiName]);
    else
        errorLog(['No orientation/location information found for ' s.NiftiName]);
    end
elseif sign(pos-R(iSL,4)) ~= signSL % same direction?
    R(:,3) = -R(:,3);
end


%% subfunction: check ImagePostionPatient from multiple slices/volumes
function [err, nSL, sliceN, isTZ] = checkImagePosition(ipp, gantryTilt)
a = diff(sort(ipp));
tol = max(a)/100; % max(a) close to SliceThichness. 1% arbituary
if nargin>1 && gantryTilt, tol = tol * 10; end % arbituary
nSL = sum(a > tol) + 1;
err = ''; sliceN = []; isTZ = false;
nVol = numel(ipp) / nSL;
if mod(nVol,1), err = 'Missing file(s) detected'; return; end
if nSL<2, return; end

isTZ = nVol>1 && all(abs(diff(ipp(1:nVol))) < tol);
if isTZ % Philips XYTZ
    a = ipp(1:nVol:end);
    b = reshape(ipp, nVol, nSL);
else
    a = ipp(1:nSL);
    b = reshape(ipp, nSL, nVol)';
end
[~, sliceN] = sort(a); % no descend since wrong for PAR/singleDicom
if any(abs(diff(a,2))>tol), err = 'Inconsistent slice spacing'; return; end
if nVol>1
    b = diff(b);
    if any(abs(b(:))>tol), err = 'Irregular slice order'; return; end
end
