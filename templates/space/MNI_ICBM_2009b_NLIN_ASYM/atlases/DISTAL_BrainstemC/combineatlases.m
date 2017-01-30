a1=load('atlas_index1');

a2=load('atlas_index2');




atlases=a1.atlases;

offset=length(a1.atlases.names);
for natl=1:length(a2.atlases.names);
    atlases.names{offset+natl}=a2.atlases.names{natl};
    
    atlases.types(offset+natl)=a2.atlases.types(natl);
    atlases.colors(offset+natl)=a2.atlases.colors(natl);
    for side=1:2
    atlases.fv{offset+natl,side}=a2.atlases.fv{natl,side};
    atlases.cdat{offset+natl,side}=a2.atlases.cdat{natl,side};
    atlases.XYZ{offset+natl,side}=a2.atlases.XYZ{natl,side};
    atlases.pixdim{offset+natl,side}=a2.atlases.pixdim{natl,side};
    atlases.colorc{offset+natl,side}=a2.atlases.colorc{natl,side};
    atlases.normals{offset+natl,side}=a2.atlases.normals{natl,side};
    end
    
end
save('atlas_index','atlases','-v7.3');
