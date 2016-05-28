function afv=ea_concatfv(fv)

afv.faces=[];
afv.vertices=[];

for f=1:length(fv)
   afv.faces=[afv.faces;fv(f).faces+length(afv.vertices)];
   afv.vertices=[afv.vertices;fv(f).vertices];
end


if isfield(fv,'facevertexcdata')
    afv.facevertexcdata=[];
    for f=1:length(fv) 
        afv.facevertexcdata=[afv.facevertexcdata;fv(f).facevertexcdata];
    end
    
end