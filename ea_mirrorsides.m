function elstruct=ea_mirrorsides(elstruct)


% function that will duplicate elstruct with mirrored sides
if isstruct(elstruct) % proper elstruct
dim=length(elstruct);
for el=1:dim
    elstruct(el+dim)=elstruct(el);
    for side=1:2
        switch side
            case 1
                oside=2;
            case 2
                oside=1;
        end
        elstruct(el+dim).coords_mm{side}=elstruct(el).coords_mm{oside}; % copy contralateral
        elstruct(el+dim).coords_mm{side}(:,1)=elstruct(el+dim).coords_mm{side}(:,1).*-1;

        elstruct(el+dim).trajectory{side}=elstruct(el).trajectory{oside}; % copy contralateral
        elstruct(el+dim).trajectory{side}(:,1)=elstruct(el+dim).trajectory{side}(:,1).*-1;
        
        elstruct(el+dim).markers(side)=elstruct(el).markers(oside); % copy contralateral
        elstruct(el+dim).markers(side).head(1)=elstruct(el+dim).markers(side).head(1).*-1;
        elstruct(el+dim).markers(side).tail(1)=elstruct(el+dim).markers(side).tail(1).*-1;
        elstruct(el+dim).markers(side).x(1)=elstruct(el+dim).markers(side).x(1).*-1;
        elstruct(el+dim).markers(side).y(1)=elstruct(el+dim).markers(side).y(1).*-1;
        elstruct(el+dim).name=[elstruct(el+dim).name,'_mirrored'];
        
    end

end

else % isomatrix
    isom=elstruct;
    isom{1}=[elstruct{1};elstruct{2}];
    isom{2}=[elstruct{2};elstruct{1}];
    elstruct=isom;
end