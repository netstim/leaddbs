function fibhandle=ea_showfiber(fibers,fibidx)

fibersnew=mat2cell(fibers(1:3,:)',fibidx);

fibersnew=cellfun(@downsample,fibersnew,repmat({5},length(fibersnew),1),'UniformOutput',0);

ea_dispercent(0,'Determine Directions');
k = 1;
while k <= length(fibersnew)
    ea_dispercent(k/length(fibersnew))
    if length(fibersnew{k}(:,1)) > 1
        fibersdiff{k} = abs(diff(fibersnew{k}));
        fibersdiff{k} = vertcat(fibersdiff{k},fibersdiff{k}(end,:));
        for l = 1:length(fibersdiff{k}(:,1))
            fibersdiff{k}(l,:) = fibersdiff{k}(l,:)/norm(fibersdiff{k}(l,:));
        end
        k = k+1;
    else
        fibersnew(k) = [];
    end
end
ea_dispercent(1,'end');

fibhandle = streamtube(fibersnew,0.25);
set(fibhandle(:),'EdgeAlpha',0)
set(fibhandle(:),'FaceAlpha',0.25)
set(fibhandle(:),'CDataMapping','direct')
set(fibhandle(:),'EdgeAlpha',0);

ea_dispercent(0,'Adding color information');
for k = 1:length(fibhandle)
    thiscol = fibersdiff{k};
    thiscol = repmat(thiscol,1,1,length(fibhandle(k).ZData(1,:,1)));
    thiscol = permute(thiscol,[1 3 2]);
    set(fibhandle(k),'CData',thiscol)
    ea_dispercent(k/length(fibhandle))
end
ea_dispercent(1,'end');

% we could be done here - but now lets concatenate the tracts for faster
% visualization:
afv=ea_concatfv(fibhandle,0,0.1);
delete(fibhandle);

%dafv=reducepatch(afv,0.2,'fast');

fibhandle=patch('Faces',afv.faces,'Vertices',afv.vertices,'FaceVertexCData',afv.facevertexcdata,'EdgeColor','none','FaceAlpha',0.25,'CDataMapping','direct','FaceColor','flat');
keyboard

ea_dispercent(1,'end');