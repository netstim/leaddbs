function activeidx=ea_getactiveidx(S,side,centroids,mesh,elfv,elspec,meshregions)
% Function used to find the active contacts when creating the first headmodel
%
% Adapted for Aleva directSTIM since the centroids of directSTIM is outside
% the mesh of the concave contacts.

emesh=[mesh.tet,meshregions];
nmesh=mesh.pnt;
load([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname])

% init activeidx:

for s=1:4
    for c=1:electrode.numContacts
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
    if contains(elspec.elmodel, 'Aleva')
        % first check if whether contact or insulator
        thiscompsnodes=emesh(emesh(1:end,5)==reg,1:4); % get this components nodes
        Ntc=size(thiscompsnodes,1);
        if Ntc>2500 % take a representative sample if whole points too large.
            thiscompsnodes=thiscompsnodes(round(linspace(1,Ntc,2500)),:);
        end
        testnode=mean(cat(3,nmesh(thiscompsnodes(:,1),:),nmesh(thiscompsnodes(:,2),:),nmesh(thiscompsnodes(:,3),:),nmesh(thiscompsnodes(:,4),:)),3);
    end

    % in this case only check for contacts
    % a - check contacts:
    for con=active
        if ~contains(elspec.elmodel, 'Aleva')
            convin=ea_intriangulation(elfv(con).vertices,elfv(con).faces,centroids(reg,:));
            thiscompsnodes=emesh(emesh(1:end,5)==reg,1:4); % get this components nodes
            dirinodes=nmesh(thiscompsnodes,:);
            testnode=ea_nudgedirinodes(dirinodes,centroids(reg,:));
        end

        in = double(ea_intriangulation(elfv(con).vertices,elfv(con).faces,testnode));
        if contains(elspec.elmodel, 'Aleva')
            captured = mean(in)>0.7;
        else
            captured = convin && mean(in)>0.7;
        end
        
        if  captured
            % we captured an active contact. need to assign to correct
            % source and polarity

            for source=1:4
                if S.([sidec,'s',num2str(source)]).amp % then this active contact could be from this source since source is active
                    if S.([sidec,'s',num2str(source)]).(['k',num2str(con)]).perc % current captured contact is from this source
                        activeidx(source).con(con).ix=[activeidx(source).con(con).ix;unique(thiscompsnodes(:))];
                        activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con)]).pol;
                        activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con)]).perc;
                    end
                end
            end

            disp(['Region ',num2str(reg),' captured by contact material.']);
            break
        end
    end
end
