%% extractElectrodePointclouds - Preprocess CT Data for Electrode Artifacts
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu

function [elecsPointcloudStruct, brainMask] = extractElectrodePointclouds(niiCT, varargin)
    disp(['Voxel size in elecsPointcloudStruct: ' num2str(niiCT.voxsize')]);

    % CONSTANTS
    LAMBDA_1 = 25;  % elec latent space length [mm]

    %% Optional Arguments and Default Values
    argParser = inputParser();
    argParser.KeepUnmatched = true;

    argParser.addParameter('noMask', false); % for phantom studies where no brain is present in data
    argParser.addParameter('brainMask', ''); % for manually providing brain mask (binary segmentation image file path)
    argParser.addParameter('medtronicXMLPlan', '', @(x)(ischar(x)));
    argParser.addParameter('metalThreshold', 800, @(x)(isnumeric(x)));
    
    argParser.parse(varargin{:});
    args = argParser.Results;

    % Check if the CT is in "standard" range [-1024 4096], if not
    % assume a 1024 offset was added to make it strictly positive and
    % handle this here
    if(min(niiCT.img(:)) >= 0)
        METAL_THRESHOLD = args.metalThreshold + 1024;
    else
        METAL_THRESHOLD = args.metalThreshold;
    end
    
    %% determine brainMask
    if(args.noMask)
        disp('Using NO brain mask as "noMask" parameter was set...');
        brainMask = logical(niiCT.voxdim); % for Phantom etc. (process whole image)
    elseif(~isempty(args.brainMask))
        disp('Using brain mask provied by parameter "brainMask"...');
        niiBrainMask = NiftiSeg(args.brainMask);
        brainMask = niiBrainMask.img;
    else
        disp('Extracting convex hull brain mask...');
        [brainMask, ~] = extractBrainConvHull(niiCT);
    end
    
    %% detect metal artifacts inside the brain (hopefully representing electrodes)
    disp(['Thresholding ' niiCT.filepath  ' for metal with METAL_THRESHOLD = ' num2str(METAL_THRESHOLD) '...']);
    maskedImg = niiCT.img;

   %  [xx,yy,zz] = ndgrid(-10:10);
    %structEle = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2.5 / sqrt(max(niiCT.voxsize));
    %brainMask = imerode(imerode(brainMask,structEle),structEle);

    maskedImg(~(brainMask)) = NaN;
    threImg = (maskedImg > METAL_THRESHOLD);
    %% largest connectet components of metal inside the brain represent electrodes
    cc = bwconncomp(threImg,26);
    disp([num2str(cc.NumObjects) ' potential metal components detected within brain.']);
    
    ccProps = regionprops(cc, niiCT.img, 'Area', 'PixelIdxList', 'PixelList', 'PixelValues', 'BoundingBox'); % TODO not really needed as that info is already in cc
    [areas, idxs] = sort([ccProps.Area], 'descend'); % sort by size
    
    %%  try to guess the number and idxs of electrodes in the image
    elecIdxs = [];
    minVoxelNumber =  (1.2 * (1.27/2))^2 * pi * 40 / prod(niiCT.voxsize); % assumin at least 40mm in brain and 20% partial voluming
    maxVoxelNumber =  (3 * (1.27/2))^2 * pi * 80 / prod(niiCT.voxsize);  % assumin 80mm in brain and 300% partial voluming 
   % maxVoxelNumber = Inf; % FIXME
    % DEBUG: figure, scatterMatrix3(ccProps(1).PixelList)
    largeComponents = areas(areas >= minVoxelNumber & areas <= maxVoxelNumber); % Voxels
    componentIdxs = idxs(areas >= minVoxelNumber & areas <= maxVoxelNumber);
    
    for i = 1:length(largeComponents)
        [~,~,latent] = pca(ccProps(componentIdxs(i)).PixelList .* repmat(niiCT.voxsize', length(ccProps(componentIdxs(i)).PixelList) ,1));
        if(length(latent) < 3)
            continue
        end
        latent = sqrt(latent) * 2; % axes length == variance, * 2 for full (instead half) axes, note that the variance is not the extrem value!
       % an electrode has large variance in one direction and about the same in the
       % two others FIXME magic constants! 
      % pointCloudExtend = norm(max(ccProps(componentIdxs(i)).PixelList) - min(ccProps(componentIdxs(i)).PixelList))
       lowerAxesLength = sort(latent(2:3));
       if(latent(1) > LAMBDA_1 && latent(1) / mean(latent(2:3)) > 10 && lowerAxesLength(2) / (lowerAxesLength(1)+0.001) < 8) % 
           elecIdxs(end+1) = componentIdxs(i);  %#ok<AGROW>
       end
    end
    nElecs = length(elecIdxs);
    disp(['Guessing that ' num2str(nElecs) ' of them are Electrodes...']);
    
    if(nElecs == 0)
        if(METAL_THRESHOLD < 3000)
            disp('Something is weird with your CT data...  Trying again with higher metal threshold. ')
            [elecsPointcloudStruct, brainMask] = extractElectrodePointclouds(niiCT, 'brainMask', args.brainMask, 'metalThreshold', METAL_THRESHOLD * 1.2, 'medtronicXMLPlan', args.medtronicXMLPlan);
            return;
        else
            %% We tried hard but  didn't find an object that looks like an electrode in a reasonalbe HU range, notify the user and quit
            error(['NO electrode artifact found within brain mask. Did you supply a post-op brain CT image? \n Try the ''no mask'' parameter in case of phantom scans without brain. \n Try providing an externally created brain mask using the "brainMask" parameter in other cases.']); %#ok<NBRAK>
        end
    end
    %% if we have an xmlElectrodeDefinition, try to find the electrodes
    % specified there
    if(exist((args.medtronicXMLPlan), 'file'))
        reportedElecs = readMedtronicXMLTrajectory((args.medtronicXMLPlan));
        XMLDefintionPresent = true;
        nReportedElecs = length(reportedElecs.trajects);
        disp(['xmlElectrodeDefinition ist reporting ' num2str(nReportedElecs) ' electrodes defined']);
        if(nReportedElecs ~= nElecs)
            warning('extractElectrodePointclouds::Number of Electrodes specified in XML defintion does not match number of found electrodes!');
        end
    else
        XMLDefintionPresent = false;
        disp(['No xmlElectrodeDefition given. Guessing that there are ' num2str(nElecs) ' electrodes in the image']);
    end
    
    %% create output struct and try to associate electrodes from xml defitions
    % if given
    elecsPointcloudStruct = struct();
    
    for i=1:nElecs
        elecsPointcloudStruct(i).pixelIdxs =ccProps(elecIdxs(i)).PixelIdxList;
        pixelList = ccProps(elecIdxs(i)).PixelList;
        pixelList(:,[1 2 3]) = pixelList(:,[2 1 3]); % Swap i,j to X,Y FIXME check this!
        elecsPointcloudStruct(i).pixelValues =ccProps(elecIdxs(i)).PixelValues;

        elecsPointcloudStruct(i).pointCloudMm =(pixelList-1) * abs(niiCT.transformationMatrix(1:3,1:3));  % minus 1 is done in the get funtions but manually here
        elecsPointcloudStruct(i).pointCloudWorld = niiCT.getNiftiWorldCoordinatesFromMatlabIdx(pixelList')';%bsxfun(@plus,(pixelList-1) * niiCT.transformationMatrix(1:3,1:3), niiCT.transformationMatrix(1:3,4)'); 
        
     %   elecsPointcloudStruct(i).surroundingPoints = setdiff(bbPointCloud, pixelList, 'rows')* abs(niiCT.transformationMatrix(1:3,1:3));
        
        elecMask = false(size(maskedImg)); % FIXME check this swaps!!
        elecMask(elecsPointcloudStruct(i).pixelIdxs ) = true; % TODO: make sure we don't have to Swap i,j to X,Y here!
        elecsPointcloudStruct(i).binaryMaskImage = elecMask;
        
        if(XMLDefintionPresent)
            %try to find associated electrode defintion 
            trajNo = getMostPropableAssociatedXmlDefinition(elecsPointcloudStruct(i).pointCloudMm);
            elecsPointcloudStruct(i).associatedXmlDefinition = reportedElecs.trajects(trajNo);
        end
    end
    
    
    %% Nested helper functions
    function trajNo = getMostPropableAssociatedXmlDefinition(pointCloud)
        noSteps = 50; % no. points to sample on a trajcet
        dists = NaN(nReportedElecs,50);
        for j=1:nReportedElecs
            % this could be a more robust by using more points on the
            % trajectory
            entry = convertMedtronicCoordToLPI(reportedElecs.trajects(j).entry, niiCT)';
            target = convertMedtronicCoordToLPI(reportedElecs.trajects(j).target, niiCT)';
            direct = (target - entry);
            len = norm(target - entry);
            direct = direct / len;
            
            %TODO distnace for coloring :-)
            steps = 1:noSteps;
            stepLen = len / noSteps;
            sampleOn = bsxfun(@plus,((stepLen .* steps)' * direct), entry);
            
            [~, dist] = dsearchn(pointCloud, ...
                sampleOn);
            
            dists(j,:) = dist;
        end
        [~, idx] = sort(mean(dists.^2,2));
        trajNo = idx(1);
    end
end