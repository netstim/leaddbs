function fibhandle=ea_showfiber(fibers,fibidx,color,fiberalpha)

if ~exist('color','var')
    color=nan;
end

if ~exist('fiberalpha','var')
    fiberalpha=0.2;
end

% Adapt tube width and fv reduce factor for rat (and other) space[s]
if ismember(ea_getspace, {'Waxholm_Space_Atlas_SD_Rat_Brain'})
    sampleFactor = 1; % Do not reduce the points along the fiber
    tubeWidth = 0.02; % Smaller tube width
    reduceFactor = 0; % No patch reduce
else
    sampleFactor = 5; % Reduce the points along the fiber by a factor of 5
    prefs = ea_prefs;
    tubeWidth = prefs.d3.fiberwidth; % Larger tube width
    reduceFactor = 0.1; % Set patch reduce factor
end

if size(fibers,1) < size(fibers,2)
    fibers = fibers';
end

fibersnew = mat2cell(fibers(:, 1:3), fibidx);

fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/sampleFactor))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);

k = 1;
while k <= length(fibersnew)
    if length(fibersnew{k}(:,1)) <= 1
        fibersnew(k) = [];
    else
        k = k+1;
    end
end
if isnan(color)
    ea_dispercent(0,'Determine Directions');
    k = 1;
    for k = 1:length(fibersnew)
        ea_dispercent(k/length(fibersnew))
            fibersdiff{k} = abs(diff(fibersnew{k}));
            fibersdiff{k} = vertcat(fibersdiff{k},fibersdiff{k}(end,:));
            for l = 1:length(fibersdiff{k}(:,1))
                fibersdiff{k}(l,:) = fibersdiff{k}(l,:)/norm(fibersdiff{k}(l,:));
            end
    end
    ea_dispercent(1,'end');
end

% Downsample fibers. TODO: Need further validation!
numFiberThreshold = 500;
if length(fibersnew) < numFiberThreshold
    idx = 1:length(fibersnew);
else
    idx = ceil(linspace(1, length(fibersnew), numFiberThreshold));
end

try
    fibhandle = streamtube(fibersnew(idx), tubeWidth);
catch
    keyboard
end
set(fibhandle(:),'CDataMapping','direct')

fprintf('\n');
if isnan(color)    
    ea_dispercent(0,'Adding color information');
    for k = 1:length(fibhandle)
        thiscolor = fibersdiff{k};
        thiscolor = repmat(thiscolor,1,1,length(fibhandle(k).ZData(1,:,1)));
        thiscolor = permute(thiscolor,[1 3 2]);
        set(fibhandle(k),'CData',thiscolor)
        ea_dispercent(k/length(fibhandle))
    end
    ea_dispercent(1,'end');
else
    ea_dispercent(0,'Adding color information');
    for k = 1:length(fibhandle)
        thiscolor = repmat(color,length(fibhandle(k).ZData(:,1)),1);
        thiscolor = repmat(thiscolor,1,1,length(fibhandle(k).ZData(1,:,1)));
        thiscolor = permute(thiscolor,[1 3 2]);
        set(fibhandle(k),'CData',thiscolor)
        ea_dispercent(k/length(fibhandle))
    end
    ea_dispercent(1,'end');
end

% we could be done here - but now lets concatenate the tracts for faster
% visualization
afv = ea_concatfv(fibhandle, 0, reduceFactor);
delete(fibhandle);

fibhandle=patch('Faces',afv.faces,'Vertices',afv.vertices,'FaceVertexCData',afv.facevertexcdata,'EdgeColor','none','FaceAlpha',fiberalpha,'CDataMapping','direct','FaceColor','flat');
