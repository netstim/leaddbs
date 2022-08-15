function afv=ea_concatfv(fv,usecork,reduce)

if ~exist('reduce','var')
    reduce=0;
end
if ~exist('usecork','var')
    usecork=0;
end
if strcmp(class(fv),'matlab.graphics.primitive.Surface')
    for f=1:length(fv)
        nfv(f)=surf2patch(fv(f),'triangles');
    end
    fv=nfv;
    clear nfv
end
if reduce
    ea_dispercent(0,'Reducing patch');
    for f=1:length(fv)
        kfv=reducepatch(fv(f),reduce,'fast');
        [~,ic]=ismember(kfv.vertices,fv(f).vertices,'rows');
        kfv.facevertexcdata=fv(f).facevertexcdata(ic,:);
        fv(f)=kfv;
        ea_dispercent(f/length(fv));
    end
    ea_dispercent(1,'end');
end

if ~usecork
    vertlen=cellfun(@length,{fv.vertices});
    facelen=cellfun(@length,{fv.faces});
    if isfield(fv,'facevertexcdata')
        facecollen=cellfun(@length,{fv.facevertexcdata});
    end
    afv.faces=zeros(sum(facelen),3);
    afv.vertices=zeros(sum(vertlen),3);
    foffset=1;
    voffset=1;
    fcoffset=1;

    fprintf('\n\n');
    ea_dispercent(0,'Concatenating patch');
    for f=1:length(fv)
        fsize=size(fv(f).faces,1);
        vsize=size(fv(f).vertices,1);
        afv.faces(foffset:foffset+fsize-1,:)=...
            fv(f).faces(:,1:3)+voffset-1;
        afv.vertices(voffset:voffset+vsize-1,:)=...
            fv(f).vertices;

        foffset=foffset+fsize;
        voffset=voffset+vsize;

        if isfield(fv(f),'facevertexcdata')
            if size(fv(f).facevertexcdata,1)<size(fv(f).facevertexcdata,2)
                fv(f).facevertexcdata=fv(f).facevertexcdata';
            end

            fcsize=size(fv(f).facevertexcdata,1);
            fcdim=size(fv(f).facevertexcdata,2);
            if ~isfield(afv,'facevertexcdata')
                afv.facevertexcdata=zeros(sum(facecollen),fcdim);
            end

            afv.facevertexcdata(fcoffset:fcoffset+fcsize-1,:)=...
                fv(f).facevertexcdata;

            fcoffset=fcoffset+fcsize;
        end
        ea_dispercent(f/length(fv));
    end
    ea_dispercent(1,'end');
else
    for f=1:length(fv)
        if f==1
            afv=fv(1);
        else
            [afv.vertices,afv.faces]=surfboolean(fv(f).vertices,fv(f).faces,'xor',afv.vertices,afv.faces);
        end
    end
end
