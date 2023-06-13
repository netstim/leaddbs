function [elrender,ellabel]=ea_showcorticalstrip(resultfig,elstruct,pt,options)
% This function renders the electrode as defined by options.elspec and
% coords_mm.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

coords_mm=elstruct.coords_mm;
trajectory=elstruct.trajectory;

if ~isfield(elstruct,'elmodel') % usually, elspec is defined by the GUI. In case of group analyses, for each patient, a different electrode model can be selected for rendering.
    elspec=options.elspec;
else % if elspec is defined for each electrode, overwrite options-struct settings here.
    o=ea_resolve_elspec(elstruct);
    elspec=o.elspec; clear o
end

if ~isfield(elstruct,'activecontacts')
    elstruct.activecontacts{1}=zeros(elspec.numel,1);
    elstruct.activecontacts{2}=zeros(elspec.numel,1);
end

jetlist=jet;

for iside=1:length(options.sides)
    side=options.sides(iside);

    trajvector=mean(diff(trajectory{side}));

    trajvector=trajvector/norm(trajvector);

    try
        startpoint=trajectory{side}(1,:)-(2*(coords_mm{side}(1,:)-trajectory{side}(1,:)));
    end

    if options.d3.elrendering<3
        set(0,'CurrentFigure',resultfig);
        % draw patientname
        lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));
        ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);

        % draw trajectory
        cnt=1;
        load([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname])
        A=[electrode.head_position,1;
            electrode.tail_position,1
            electrode.x_position,1
            electrode.y_position,1]; % points in model
        redomarkers=0;

        if ~isfield(elstruct,'markers') % backward compatibility to old electrode format
            redomarkers=1;
        else
            if isempty(elstruct.markers)
                redomarkers=1;
            end
        end

        if redomarkers
            for iside2=1:length(options.sides)
                side2=options.sides(iside2);

                elstruct.markers(side2).head=elstruct.coords_mm{side2}(1,:);
                elstruct.markers(side2).tail=elstruct.coords_mm{side2}(4,:);

                [xunitv, yunitv] = ea_calcxy(elstruct.markers(side2).head, elstruct.markers(side2).tail);
                elstruct.markers(side2).x = elstruct.coords_mm{side2}(1,:) + xunitv*(options.elspec.lead_diameter/2);
                elstruct.markers(side2).y = elstruct.coords_mm{side2}(1,:) + yunitv*(options.elspec.lead_diameter/2);
            end
        end

        B=[elstruct.markers(side).head,1;
           elstruct.markers(side).tail,1;
           elstruct.markers(side).x,1;
           elstruct.markers(side).y,1];
        setappdata(resultfig,'elstruct',elstruct);
        setappdata(resultfig,'elspec',elspec);

        X = mldivide(A,B); X=X';

        if options.d3.elrendering==2 % show a transparent electrode.
            aData=0.1;
        elseif options.d3.elrendering==1 % show a solid electrode.
            aData=1;
        end

        for ins=1:length(electrode.insulation)
            electrode.insulation(ins).vertices=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
            electrode.insulation(ins).vertices=electrode.insulation(ins).vertices(1:3,:)';
            elrender{side}(cnt)=patch(electrode.insulation(ins));

            if isfield(elstruct,'group')
                usecolor=elstruct.groupcolors(elstruct.group,:);
            else
                usecolor=elspec.lead_color;
            end
            ea_specsurf(elrender{side}(cnt),usecolor,aData);
            cnt=cnt+1;
        end

        for con=1:length(electrode.contacts)
            electrode.contacts(con).vertices=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
            electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';
            elrender{side}(cnt)=patch(electrode.contacts(con));

            if options.d3.hlactivecontacts && ismember(con,find(elstruct.activecontacts{side})) % make active red contact without transparency
                ea_specsurf(elrender{side}(cnt),[0.8,0.2,0.2],1);
            else
                ea_specsurf(elrender{side}(cnt),elspec.contact_color,aData);
            end
            cnt=cnt+1;
        end
    else % simply draw pointcloud
        shifthalfup=0;
        % check if isomatrix needs to be expanded from single vector by using stimparams:
        try
            if size(options.d3.isomatrix{1},2)==elspec.numel-1 % number of contact pairs
                shifthalfup=1;
            elseif size(options.d3.isomatrix{1},2)==elspec.numel % number of contacts
                shifthalfup=0;
            else
                ea_cprintf('CmdWinErrors', 'Be careful! Isomatrix might have wrong size, or numbers of contacts are not consistent across patients.\n');
            end
        end

        if options.d3.prolong_electrode
            startpoint=trajectory{side}(1,:)-(options.d3.prolong_electrode*(coords_mm{side}(1,:)-trajectory{side}(1,:)));
        else
            startpoint=trajectory{side}(1,:);
        end

        set(0,'CurrentFigure',resultfig);

        % draw patientname
        lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));

        %ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);
        ellabel=nan;
        pcnt=1;

        % draw contacts
        try
            minval=abs(min(options.d3.isomatrix{side}(:)));
            maxval=max(options.d3.isomatrix{side}(:));
        end

        for cntct=1:elspec.numel-shifthalfup
            if (options.d3.showactivecontacts && ismember(cntct,find(elstruct.activecontacts{side}))) || (options.d3.showpassivecontacts && ~ismember(cntct,elstruct.activecontacts{side}))
                if options.d3.hlactivecontacts && ismember(cntct,find(elstruct.activecontacts{side})) % make active red contact without transparency
                    useedgecolor=[0.8,0.5,0.5];
                    ms=10;
                elseif options.d3.hlactivecontacts && ~ismember(cntct,find(elstruct.activecontacts{side})) % make inactive grey and smaller contact without transparency
                    useedgecolor='none';
                    ms=5;
                else
                    useedgecolor='none';
                    ms=10;
                end

                % define color
                if options.d3.colorpointcloud
                    % draw contacts as colored cloud defined by isomatrix.
                    if ~isnan(options.d3.isomatrix{side}(pt,cntct))
                        usefacecolor=((options.d3.isomatrix{side}(pt,cntct)+minval)/(maxval+minval))*64;
                        usefacecolor=ind2rgb(round(usefacecolor),jetlist);
                    else
                        usefacecolor=nan; % won't draw the point then.
                    end
                else
                    if isfield(elstruct,'group')
                        usefacecolor=elstruct.groupcolors(elstruct.group,:);
                    else
                        usefacecolor=elspec.contact_color;
                    end
                end

                if ~isnan(usefacecolor)
                    set(0,'CurrentFigure',resultfig);
                    if ~shifthalfup
                        elrender{side}(pcnt)=plot3(coords_mm{side}(cntct,1),coords_mm{side}(cntct,2),coords_mm{side}(cntct,3),'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                        pcnt=pcnt+1;
                    else

                        elrender{side}(pcnt)=plot3(mean([coords_mm{side}(cntct,1),coords_mm{side}(cntct+1,1)]),...
                            mean([coords_mm{side}(cntct,2),coords_mm{side}(cntct+1,2)]),...
                            mean([coords_mm{side}(cntct,3),coords_mm{side}(cntct+1,3)]),...
                            'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                        drawnow;
                        pcnt=pcnt+1;
                    end
                end
                hold on
            end
        end
    end
end

if ~exist('elrender','var')
    elrender=nan;
end
