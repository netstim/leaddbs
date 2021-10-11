function activeidx=ea_getactiveidx(S,side,centroids,mesh,elfv,elspec,meshregions)

emesh=[mesh.tet,meshregions];
nmesh=mesh.pnt;
load([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname])

% init activeidx:

for s=1:4
    for c=1:electrode.numel
        activeidx(s).con(c).ix=[];
        activeidx(s).con(c).perc=0;
        activeidx(s).con(c).pol=0;
    end
end


active=find(S.activecontacts{side});

switch side
    case 1
        sidec='R';
    case 2
        sidec='L';
end

for reg=1:size(centroids,1)
    % in this case only check for contacts
    
    % a - check contacts:
    for con=active;
        convin=ea_intriangulation(elfv(con).vertices,elfv(con).faces,centroids(reg,:));
        
        thiscompsnodes=emesh(emesh(1:end,5)==reg,1:4); % get this components nodes
        dirinodes=nmesh(thiscompsnodes,:);
        dirinodes=ea_nudgedirinodes(dirinodes,centroids(reg,:));
        in=double(ea_intriangulation(elfv(con).vertices,elfv(con).faces,dirinodes));
        
        if convin && mean(in)>0.7
                
                % we captured an active contact. need to assign to correct
                % source and polarity
                
                for source=1:4
                    if S.([sidec,'s',num2str(source)]).amp % then this active contact could be from this source since source is active
                        if S.([sidec,'s',num2str(source)]).(['k',num2str(con+electrode.numel*(side-1)-1)]).perc % current captured contact is from this source
                            activeidx(source).con(con).ix=[activeidx(source).con(con).ix;unique(thiscompsnodes(:))];
                            activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con+electrode.numel*(side-1)-1)]).pol;
                            activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con+electrode.numel*(side-1)-1)]).perc;
                        end
                    end
                end
                
                
            
            disp(['Region ',num2str(reg),' captured by contact material.']);
            break
        end
    end
    
end