function elstruct=ea_mirrorsides(elstruct)


% function that will duplicate elstruct with mirrored sides
if isstruct(elstruct) % proper elstruct
    dim=length(elstruct);


    % 1st loop: Aggregate points to warp!! Important, this loop needs to be
    % structured exactly as the one below!
    cnt=1;
    for el=1:dim
        elstruct(el+dim)=elstruct(el);


        for side=1:2
            switch side
                case 1
                    oside=2;
                case 2
                    oside=1;
            end
            elstruct(el+dim).markers(side)=elstruct(el).markers(oside); % copy contralateral
            towarp(cnt,:)=elstruct(el+dim).markers(side).head; cnt=cnt+1;
            towarp(cnt,:)=elstruct(el+dim).markers(side).tail; cnt=cnt+1;
            towarp(cnt,:)=elstruct(el+dim).markers(side).x; cnt=cnt+1;
            towarp(cnt,:)=elstruct(el+dim).markers(side).y; cnt=cnt+1;

        end
        elstruct(el+dim).name=[elstruct(el+dim).name,'_mirrored'];
    end

    % do warp:

    disp('Nonlinearly mirroring electrode coordinates to the other hemisphere...');

    towarp=ea_flip_lr_nonlinear(towarp);


    % 2nd loop: feed in warps !! Important, this loop needs to be
    % structured exactly as the one above!
    cnt=1;
    for el=1:dim
        for side=1:2
            elstruct(el+dim).markers(side).head=towarp(cnt,:); cnt=cnt+1;
            elstruct(el+dim).markers(side).tail=towarp(cnt,:); cnt=cnt+1;
            elstruct(el+dim).markers(side).x=towarp(cnt,:); cnt=cnt+1;
            elstruct(el+dim).markers(side).y=towarp(cnt,:); cnt=cnt+1;
        end

        % now re-resolve flipped coords:
        options.elmodel=elstruct(el+dim).elmodel;
        options=ea_resolve_elspec(options);
        [coords,trajectory,markers]=ea_resolvecoords(elstruct(el+dim).markers,options);

        elstruct(el+dim).coords_mm=coords;
        elstruct(el+dim).trajectory=trajectory;


    end

else % isomatrix
    if isempty(elstruct)
        return
    end
    for is=1:length(elstruct)
        isom{is}=elstruct{is}; % elstruct variable is the isomatrix here!
        isom{is}{1}=[elstruct{is}{1};elstruct{is}{2}];
        isom{is}{2}=[elstruct{is}{2};elstruct{is}{1}];
    end
    elstruct=isom; % that's why we need to export isom as fake elstruct variable.
end



