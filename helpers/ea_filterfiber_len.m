function fiberFiltered = ea_filterfiber_len(ftr, fiblen, coords, S)
% Filter fibers based on the fiber length threshold and active contacts
% coordinates

% Load stimparameters in case needed
if ischar(S) && isfile(S)
    load(S, 'S');
end

disp('Collecting stimulation parameters...')

% Active contacts indices
activeContacts = cell(size(S.activecontacts));
for i=1:length(activeContacts)
    activeContacts{i} = find(S.activecontacts{i});
end

% Active contacts coordinates
stimCoords = cell(size(coords));
for i=1:length(stimCoords)
    stimCoords{i} = mean(coords{i}(activeContacts{i},:),1);
end

% Load fiber connectome
if ischar(ftr) && isfile(ftr)
    disp('Loading fibers...');
    ftr = load(ftr);
end

if ~iscell(ftr)
    ftr = {ftr};
end

fiberFiltered = cell(size(ftr));
for i=1:length(fiberFiltered)
    fiberFiltered{i}.fibers = [];
    fiberFiltered{i}.idx = [];
end

% Inline function to calculate the mean neighbouring distance
meandist = @(x) mean(sqrt(sum(diff(x).^2,2)));

% Check if fibers pass through the ROI
for i=1:length(ftr)
    if ~isempty(ftr{i})
        % Prepare stard and end index of each fiber
        endIndex = cumsum(ftr{i}.idx);
        startIndex = [1;endIndex+1];
        startIndex = startIndex(1:end-1);
        fibIndex = arrayfun(@(x,y) [x:y]', startIndex, endIndex, 'Uni', 0);
        
        % Calculate numbers of points to keep
        voxsize = mean(cellfun(@(x) meandist(ftr{i}.fibers(x,1:3)), fibIndex));
        numPoints = ceil(fiblen/voxsize);
        for f=1:length(fibIndex)
            if length(fibIndex{f})>numPoints
                % Put the center of the filtered fiber at the stim target
                closestInd = dsearchn(ftr{i}.fibers(fibIndex{f},1:3), stimCoords{i});
                if closestInd<floor(numPoints/2)
                    selectInd = 1:numPoints;
                elseif (length(fibIndex{f})-closestInd)<ceil(numPoints/2)
                    selectInd = length(fibIndex{f})-numPoints+1:length(fibIndex{f});
                elseif mod(numPoints, 2)
                    selectInd = (closestInd-floor(numPoints/2)):(closestInd+floor(numPoints/2));
                else
                    selectInd = (closestInd-numPoints/2+1):(closestInd+numPoints/2);
                end
                % Select fibers
                fiberFiltered{i}.fibers = [fiberFiltered{i}.fibers; ftr{i}.fibers(fibIndex{f}(selectInd),:)];
                fiberFiltered{i}.idx = [fiberFiltered{i}.idx; length(selectInd)];
            end
        end
    end
end

disp('Finished!');


function trim