function out = ea_near_line(X,Y)
% This function finds the nearest start and end points of Y in X. Then gets
% the portion of X between those points that is nearest to Y.
%
% Each row of out contains the index of X corresponding to the row in Y


% generate indexes for each seperate contour
dX = diff((sqrt(sum(diff(X).^2,2))));
jumpind = [0 0 dX'>1];
blocks = cumsum(jumpind);

% variable to save the minor distance from drawing to contour
Dglobal = inf;

for i=0:blocks(end) % iterate over contours
    
    subX = X(blocks==i,:);
    
    idx0 = knnsearch(subX, Y([1,end],:)); % find first and last points of the drawing in the contour
    
    % two possible ways to join the two points
    idx{1} = min(idx0):max(idx0);
    idx{2} = [1:min(idx0) max(idx0):length(subX)];
    
    % calculate the distance between drawing and contour
    [~, D] = knnsearch(subX(idx{1},:),Y);
    d(1) = sum(D);
    [~, D] = knnsearch(subX(idx{2},:),Y);
    d(2) = sum(D);
    
    [Dlocal, idmin] = min(d); % get minimum
    
    if Dlocal < Dglobal
        Dglobal = Dlocal; % save new minimum
        out = knnsearch(X,subX(idx{idmin},:)); % get the indexes in original X array
        mid = min(idx0); % save mid index for later
        aux = idmin; % save array index
    end
    
end

% correct the index
f = find(out == knnsearch(X,Y(1,:)));

if f == 1
    % do nothing
elseif f == mid
    out = [flip(out(1:mid)); out(mid+1:end)];
elseif f == mid+1
    out = [out(mid+1:end); out(1:mid)];
else
    out = flip(out);
end

if aux == 2 && out(end) ~= knnsearch(X,Y(end,:))
    out = [out(1:mid); flip(out(mid+1:end))];
end

% resample to get same amount of points as Y
out = out(round(linspace(1,length(out),length(Y))));

end