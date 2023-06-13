function xyzv=ea_buildxyzv(M,varnum)

% varnum can be an integer referring to the regressor # or a string
% {'contacts','pairs'} to give back contact or contact pair coordinates.

if ~exist('varnum','var')
    varnum='contacts';
end

maxNumContacts = get_maxNumContacts(M.elstruct);

if ischar(varnum) % simply return coordinate list
    switch varnum
        case 'contacts'
            ixx=1:maxNumContacts;
            meannec=0;
        case 'pairs'
            ixx=1:maxNumContacts-1;
            meannec=1;
    end
    cnt=1;
    for pt=1:length(M.patient.list)
        for side=1:2
            for clen=ixx
                if meannec
                    xyzv(cnt,:)=mean(M.elstruct(pt).coords_mm{side}(clen:clen+1,:),1);
                else
                    xyzv(cnt,:)=M.elstruct(pt).coords_mm{side}(clen,:);
                end
                cnt=cnt+1;
            end
        end
    end
else
    cnt=1;
    switch size(M.clinical.vars{varnum},2)
        case {(maxNumContacts-1)*2, maxNumContacts*2} % contacts / contact pairs
            if size(M.clinical.vars{varnum},2)==(maxNumContacts-1)*2
                ixx=1:maxNumContacts-1;
                meannec=1;
            else
                ixx=1:maxNumContacts;
                meannec=0;
            end
            
            for pt=1:length(M.patient.list)
                for side=1:2
                    for clen=ixx
                        sidec=[ixx]+(side-1)*ixx(end);
                        if meannec
                            xyzv(cnt,:)=[mean(M.elstruct(pt).coords_mm{side}(clen:clen+1,:),1),M.clinical.vars{varnum}(pt,sidec(clen))];
                        else
                            xyzv(cnt,:)=[M.elstruct(pt).coords_mm{side}(clen,:),M.clinical.vars{varnum}(pt,sidec(clen))];
                        end
                        cnt=cnt+1;
                    end
                end
            end
        case 1 % patient
            for pt=1:length(M.patient.list)
                    xyzv(pt,:)=[mean([M.elstruct(pt).coords_mm{1}(find(M.S(pt).activecontacts{1}),:);... % right coordinate
                        ea_flip_lr_nonlinear(M.elstruct(pt).coords_mm{2}(find(M.S(pt).activecontacts{2}),:))],1),... % nonlinearly flipped left coordinate
                        M.clinical.vars{varnum}(pt)]; % value
            end
        case 2 % hemisphere
            cnt=1;
            for pt=1:length(M.patient.list)
                    xyzv(cnt,:)=[M.elstruct(pt).coords_mm{1}(find(M.S(pt).activecontacts{1}),:),M.clinical.vars{varnum}(pt,1)]; % right hem
                    xyzv(cnt+1,:)=[ea_flip_lr_nonlinear(M.elstruct(pt).coords_mm{2}(find(M.S(pt).activecontacts{2}),:)),M.clinical.vars{varnum}(pt,2)]; % nonlinearly flipped left coordinate
                    cnt=cnt+2;
            end
    end
end


function maxNumContacts = get_maxNumContacts(elstruct)
coords = {elstruct.coords_mm};
coords = horzcat(coords{:})';
maxNumContacts = max(cellfun(@(x) size(x,1), coords));
