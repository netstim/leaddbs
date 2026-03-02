function [structures] = ea_updatecortex(options,resultfig,sides,structures,labelidx,alpha,graycolor)

cortex = getappdata(resultfig,'cortex');
annot = getappdata(resultfig,'annot');

if ~exist('alpha','var') || isempty(alpha)
    alpha = options.prefs.d3.cortexalpha;
end

if ~exist('graycolor','var')
    graycolor = 0;
end

if iscell(labelidx{1}) % both sides
   labelstruct = labelidx;
end
for s = sides
    if exist('labelstruct','var')
        labelidx = labelstruct{s};
    end
    vertices=cortex{s}.Vertices;
    annot(s).adat=ones(size(vertices,1),1)*alpha;

    % Choose gyri
    % labels = cell2mat(arrayfun(@(x) find(x==annot(s).label),annot(s).colortable.table(:,5),'uni',0));
    structidx = arrayfun(@(x) find(annot(s).label==annot(s).colortable.table(x,5)),[labelidx{:}],'uni',0);
    for i=1:length(structidx)
        colorindex = structidx{i};
        %annot(s).cdat(colorindex,:) = repmat(annot(s).colortable.table(i,1:3)/256,[length(colorindex),1]);
        annot(s).adat(colorindex,:) = alpha;
    end

    labelsoff = setdiff(1:length(annot(s).colortable.table),cell2mat(labelidx));
    invisidx = arrayfun(@(x) find(annot(s).label==annot(s).colortable.table(x,5)),labelsoff,'uni',0);

    % remove unlabeled
    unl = annot(s).cdat==0.65;
    invisidx{1} = find(unl(:,1));

    for i=1:length(invisidx)
        %annot(s).cdat(invisidx{i},:) = repmat(annot(s).colortable.table(i,1:3)/256,[length(colorindex),1]);
        annot(s).adat(invisidx{i},:) = 0;
    end

    %annot(side).adat(annot(side).adat==0.1111)=0;
    if ~graycolor
        set(cortex{s},'FaceVertexCData',annot(s).cdat);
    else
        set(cortex{s},'FaceVertexCData',repmat([0.65, 0.65, 0.65], length(cortex{s}.Vertices), 1));
    end
    set(cortex{s},'FaceVertexAlphaData',annot(s).adat,'FaceAlpha','interp','AlphaDataMapping','none');

end
    setappdata(resultfig,'cortex',cortex);
    setappdata(resultfig,'annot',annot);
end