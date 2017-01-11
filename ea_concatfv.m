function afv=ea_concatfv(fv,usecork)

if ~exist('usecork','var')
    usecork=0;
end


if ~usecork
    afv.faces=[];
    afv.vertices=[];
    
    for f=1:length(fv)
        afv.faces=[afv.faces;fv(f).faces+length(afv.vertices)];
        afv.vertices=[afv.vertices;fv(f).vertices];
    end
    
    
    if isfield(fv,'facevertexcdata')
        afv.facevertexcdata=[];
        for f=1:length(fv)
            
            if size(fv(f).facevertexcdata,1)<size(fv(f).facevertexcdata,2)
                fv(f).facevertexcdata=fv(f).facevertexcdata';
            end
            afv.facevertexcdata=[afv.facevertexcdata;fv(f).facevertexcdata];
        end
        
    end
    
else
    
    
    for f=1:length(fv)
        if f==1
            afv=fv(1);
        else
            [afv.vertices,afv.faces]=surfboolean(fv(f).vertices,fv(f).faces,'xor',afv.vertices,afv.faces);
        end
    end
    
    
    
end