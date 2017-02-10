function xyzv=ea_buildxyzv(M,varnum)


cnt=1;
switch size(M.clinical.vars{varnum},2)
    
    case {6,8} % contacts / contact pairs
        
        if size(M.clinical.vars{varnum},2)==6;
            ixx=1:3;
        else
            ixx=1:4;
        end
        
        for pt=1:length(M.patient.list)
            for side=1:2
                for clen=ixx
                    sidec=[ixx]+(side-1)*ixx(end);

                    xyzv(cnt,:)=[mean(M.elstruct(pt).coords_mm{side}(clen:clen+1,:),1),M.clinical.vars{varnum}(pt,sidec(clen))];
                    cnt=cnt+1;
                end
            end
            
        end
end

