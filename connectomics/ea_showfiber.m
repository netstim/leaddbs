function fibhandle=ea_showfiber(fibers,fibidx,col,fiberalpha)

if ~exist('col','var')
    col=nan;
end

if ~exist('fiberalpha','var')
    fiberalpha=0.2;
end

if ~(size(fibers,1)==4)
    fibers=fibers';
end
fibersnew=mat2cell(fibers(1:3,:)',fibidx);

fibersnew = cellfun(@(f,len) f(round(linspace(1,len,round(len/5))),:), fibersnew, num2cell(cellfun(@(p) size(p,1), fibersnew)), 'UniformOutput', 0);
%fibersnew=cellfun(@downsample,fibersnew,repmat({5},length(fibersnew),1),'UniformOutput',0);

k = 1;
while k <= length(fibersnew)
    if length(fibersnew{k}(:,1)) <= 1
        fibersnew(k) = [];
    else
        k = k+1;
    end
end
if isnan(col)
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
fibhandle = streamtube(fibersnew,0.25);
set(fibhandle(:),'CDataMapping','direct')

fprintf('\n');
if isnan(col)    
    ea_dispercent(0,'Adding color information');
    for k = 1:length(fibhandle)
        thiscol = fibersdiff{k};
        thiscol = repmat(thiscol,1,1,length(fibhandle(k).ZData(1,:,1)));
        thiscol = permute(thiscol,[1 3 2]);
        set(fibhandle(k),'CData',thiscol)
        ea_dispercent(k/length(fibhandle))
    end
    ea_dispercent(1,'end');
else
    ea_dispercent(0,'Adding color information');
    for k = 1:length(fibhandle)
        thiscol = repmat(col,length(fibhandle(k).ZData(:,1)),1);
        thiscol = repmat(thiscol,1,1,length(fibhandle(k).ZData(1,:,1)));
        thiscol = permute(thiscol,[1 3 2]);
        set(fibhandle(k),'CData',thiscol)
        ea_dispercent(k/length(fibhandle))
    end
    ea_dispercent(1,'end');
end
% we could be done here - but now lets concatenate the tracts for faster
% visualization:

afv=ea_concatfv(fibhandle,0,0.5);
delete(fibhandle);

%dafv=reducepatch(afv,0.2,'fast');

fibhandle=patch('Faces',afv.faces,'Vertices',afv.vertices,'FaceVertexCData',afv.facevertexcdata,'EdgeColor','none','FaceAlpha',fiberalpha,'CDataMapping','direct','FaceColor','flat');
