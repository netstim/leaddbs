function mycols = ea_explorer_isocolors(img,vertices)

%% use this if you have RAM issues
fprintf('Sampling exact isocolors for volume.\nIf taking too long or too RAM-heavy please downsample mesh using the hullsimplify option.\n')
tic
vertices=single(vertices);
[x,y,z] = ind2sub(size(img),[1:numel(img(:))]);
coords = vertcat(x,y,z)';
allvals = single(img(:));
nanidx = find(isnan(allvals));
allvals(nanidx) = [];
coords(nanidx,:) =[];
allvals = repmat(allvals,1,size(vertices,1));
clear x y z nanidx

xdiff = coords(:,1)-vertices(:,1)';
ydiff = coords(:,2)-vertices(:,2)';
zdiff = coords(:,3)-vertices(:,3)';
norms = sqrt(xdiff.^2 + ydiff.^2 + zdiff.^2);
clear xdiff ydiff zdiff
minnorms = min(norms,[],1);

% if ~isempty(find(minnorms > 1))
%     keyboard
% end

minnorms = num2cell(minnorms);
norms = num2cell(norms,1);
minnormidx = cellfun(@ismember,norms,minnorms,'UniformOutput',false);
clear norms minnorms
minnormidx = horzcat(minnormidx{:});
selectedvals = allvals.*minnormidx;
mycols = (sum(selectedvals,1) ./ sum(minnormidx,1))';
toc
end