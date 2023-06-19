function [elrender,ellabel,eltype,eltext]=ea_showelectrode(obj,cmd,options)
% This function renders the electrode as defined by options.elspec and coords_mm.
% _______________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

switch cmd
    case 'dbs'
        resultfig=obj.plotFigureH;
        elstruct=obj.elstruct;
        pt=obj.pt;
    case 'plan'
        resultfig=obj.plotFigureH;
        elstruct=obj.plan2elstruct;
        elstruct.elmodel = obj.plan2elstruct_model;
        pt=1;
end

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

if isfield(options.d3,'pntcmap')
    cmap = options.d3.pntcmap;
elseif isfield(options.d3, 'regressorcolormap')
    cmap = options.d3.regressorcolormap;
else % use default blue to red colormap
    cmap = ea_colorgradient(length(gray), [0,0,1], [1,1,1], [1,0,0]);
end

try
    cmap=evalin('base','custom_cmap');
end

for side=options.sides
    try
        startpoint=trajectory{side}(1,:)-(2*(coords_mm{side}(1,:)-trajectory{side}(1,:)));
    catch
        keyboard
    end

    if options.d3.elrendering<3
        set(0,'CurrentFigure',resultfig);

        % draw patientname
        lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));
        %ellabel=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);

        lp=[trajectory{side}(end,1),trajectory{side}(end,2),trajectory{side}(end,3)];
        ap=[trajectory{side}(1,1),trajectory{side}(1,2),trajectory{side}(1,3)];
        lp=lp+(lp-ap);

        ellabel=text(lp(1),lp(2),lp(3),ea_underscore2space(elstruct.name),'Color',[1,1,1]);

        % draw trajectory
        cnt=1;

        err=1;
        for tries=1:2
            [X,electrode,err]=ea_mapelmodel2reco(options,elspec,elstruct,side,resultfig);
            if isfield(options,'nowrite')
                if options.nowrite % means elstruct was manually manipulated after it had been plotted, we dont want to save that to disk.
                    err = 0;
                end
            end
            if ~err
                break
            elseif ~options.d3.mirrorsides
                try
                    % Recalculate as the tolerance/precision was not
                    % satisfactory for this use case (will be loaded at tries==2)
                    if ~isfield(options,'patient_list') % single subject mode
                        % Recalc then reload depending on options.native flag
                        ea_recalc_reco([],[],[options.root,options.patientname]);
                        [coords_mm,trajectory,markers]=ea_load_reconstruction(options);
                    else
                        % Recalc then reload depending on options.native flag
                        ea_recalc_reco([],[],[options.patient_list{pt},filesep]);
                        [coords_mm,trajectory,markers]=ea_load_reconstruction(options);
                    end
                    elstruct.markers=markers;
                    elstruct.coords_mm=coords_mm;
                    elstruct.trajectory=trajectory;
                catch
                    if ~isfield(options,'patient_list') % single subject mode
                        warning(['There seems to be some inconsistency with the reconstruction of ',options.patientname,' that could not be automatically resolved. Please check data of this patient.']);
                    else
                        warning(['There seems to be some inconsistency with the reconstruction of ',options.patient_list{pt},' that could not be automatically resolved. Please check data of this patient.']);
                    end
                end
            end
        end
        if err
            if ~isfield(options,'patient_list') % single subject mode
                warning(['There seems to be some inconsistency with the reconstruction of ',options.patientname,' that could not be automatically resolved. Please check data of this patient.']);
            else
                warning(['There seems to be some inconsistency with the reconstruction of ',options.patient_list{pt},' that could not be automatically resolved. Please check data of this patient.']);
            end
        end
        if options.d3.elrendering==2 % show a transparent electrode.
            aData=0.1;
        elseif options.d3.elrendering==1 % show a solid electrode.
            aData=1;
        end

        if isfield(elstruct, 'name') && ~isempty(elstruct.name)
            nameprefix = [elstruct.name, '_'];
        else
            nameprefix = '';
        end

        for ins=1:length(electrode.insulation)
            electrode.insulation(ins).vertices=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
            electrode.insulation(ins).vertices=electrode.insulation(ins).vertices(1:3,:)';
            elrender(cnt)=patch(electrode.insulation(ins));

            elrender(cnt).Tag = [nameprefix, 'Insulation', num2str(cnt), '_Side', num2str(side)];

            if isfield(options,'sidecolor')
                switch side
                    case 1
                        usecolor=[ 0.9797    0.9831    0.5185];
                    case 2
                        usecolor=[ 0.4271    0.7082    0.8199];
                    otherwise
                        usecolor=elspec.lead_color;
                end
            else
                switch cmd
                    case 'dbs'
                        if isfield(elstruct,'group')
                            usecolor=elstruct.groupcolors(elstruct.group,:);
                        else
                            usecolor=elspec.lead_color;
                        end
                    case 'plan'
                        usecolor=obj.color;
                end
            end
            ea_specsurf(elrender(cnt),usecolor,aData,'insulation');
            cnt=cnt+1;
        end

        eltext=getappdata(resultfig,'eltext');
        [contactnames,directional]=ea_getelcontactnames(elspec,side);
        for con=1:size(coords_mm{side},1)
            % add text:
            centroid=coords_mm{side}(con,:)+0.01;

            % find intersection point S on line defined by tail and head
            Xpt = centroid;
            Ppt = elstruct.markers(side).head;
            Qpt = elstruct.markers(side).tail;

            % Calculate the unit vector for the line segment PQ
            u = (Qpt - Ppt) / norm(Qpt - Ppt);

            % Calculate the vector from P to Xpt
            v = Xpt - Ppt;

            % Calculate the distance from Xpt to the line segment PQ
            d = norm(v - dot(v,u)*u);

            % Calculate the point X on the line segment PQ that is closest to Y
            Spt = Ppt + dot(v,u)*u;

            normv=norm(centroid-Spt);
            if directional(con)
                pointfortext=centroid+0.9*((centroid-Spt)/normv);
            else
                pointfortext=centroid+1.8*((centroid-Spt)/normv);
            end
            eltext(side,con)=text(pointfortext(1),pointfortext(2),pointfortext(3),contactnames{con},'FontWeight','bold','FontSize',14,'Color',[0,0,0],'HorizontalAlignment','center','VerticalAlignment','middle');
            set(eltext(side,con), 'Visible','off');
        end
        setappdata(resultfig,'eltext',eltext);

        for con=1:length(electrode.contacts)
            electrode.contacts(con).vertices=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
            electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';
            elrender(cnt)=patch(electrode.contacts(con));
            



            elrender(cnt).Tag = [nameprefix, 'Contact', num2str(con), '_Side', num2str(side)];
            eltype(cnt)=1;
            if ~isempty(options.colorMacroContacts)
                ea_specsurf(elrender(cnt),options.colorMacroContacts(con,:),1,'metal');
            else
                if contains(nameprefix, '_mirrored_')
                    % Mirror highlight contact
                    if side == 1
                        mside = 2;
                    elseif side == 2
                        mside = 1;
                    end
                else
                    mside = side;
                end

                if options.d3.hlactivecontacts && ismember(con,find(elstruct.activecontacts{mside})) % make active red contact without transparency
                    ea_specsurf(elrender(cnt),[0.8,0.2,0.2],1,'metal');
                else
                    ea_specsurf(elrender(cnt),elspec.contact_color,aData,'metal');
                end
            end
            cnt=cnt+1;
        end

        % arrows for directional leads
        if isfield(options.prefs.d3,'showdirarrows') && options.prefs.d3.showdirarrows
            switch options.elmodel
                case {'Medtronic B33005'
                      'Medtronic B33015'
                      'Boston Scientific Vercise Directed'
                      'Boston Scientific Vercise Cartesia HX'
                      'Boston Scientific Vercise Cartesia X'
                      'Abbott Directed 6172 (short)'
                      'Abbott Directed 6173 (long)'}
                    % Marker position relative to head position along z axis
                    markerposRel = options.elspec.markerpos-electrode.head_position(3);
                    dothearrows = 1;
                otherwise
                    dothearrows = 0;
            end
            if dothearrows
                % Calc stretch factor since lead could be stretched due to non-linear transformation
                stretchfactor = norm(elstruct.markers(side).tail - elstruct.markers(side).head) / (electrode.tail_position(3)-electrode.head_position(3));
                % Direction of the lead
                unitvector = (elstruct.markers(side).tail - elstruct.markers(side).head) / norm(elstruct.markers(side).tail - elstruct.markers(side).head);
                % Calc stick location
                stxmarker = elstruct.markers(side).head + stretchfactor * markerposRel * unitvector;
                arrowtip = stxmarker + 5 * (elstruct.markers(side).y - elstruct.markers(side).head);
                elrender(cnt) = mArrow3(stxmarker,arrowtip,'color',[.3 .3 .3],'tipWidth',0.2,'tipLength',0,'stemWidth',0.2,'Tag','DirectionMarker');
                ea_specsurf(elrender(cnt),[.3 .3 .3],1,'metal');
                cnt = cnt+1;
            end
        end
    else % simply draw pointcloud
        shifthalfup=0;
        % check if isomatrix needs to be expanded from single vector by using stimparams:
        try % sometimes isomatrix not defined.
            if size(options.d3.isomatrix{1}{1},2)==elspec.numel-1 % number of contact pairs
                shifthalfup=1;
            elseif size(options.d3.isomatrix{1}{1},2)==elspec.numel % number of contacts
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
        % lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));
        % ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);

        ellabel=[];
        eltype=[];
        pcnt=1;

        % draw contacts
        try
            normalisomatrix{side}=options.d3.isomatrix{1}{side};
            normalisomatrix{side}(:)=ea_normal(normalisomatrix{side}(:));
            minval=ea_nanmin(normalisomatrix{side}(:));
            maxval=ea_nanmax(normalisomatrix{side}(:));
            %minval=-1;
            %maxval=1;
        end
        for cntct=1:elspec.numel-shifthalfup
            if (options.d3.showactivecontacts && ismember(cntct,find(elstruct.activecontacts{side}))) || (options.d3.showpassivecontacts && ~ismember(cntct,find(elstruct.activecontacts{side})))
                if options.d3.hlactivecontacts && ismember(cntct,find(elstruct.activecontacts{side})) % make active red contact without transparency
                    useedgecolor=[0.9,0.9,0.7];
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
                    if ~isnan(options.d3.isomatrix{1}{side}(pt,cntct))
                        usefacecolor=((normalisomatrix{side}(pt,cntct)-minval)/(maxval-minval))*(length(cmap)-1);
                        % % Add some contrast (remove these lines for linear mapping)
                        % usefacecolor=usefacecolor-20;
                        % usefacecolor(usefacecolor<1)=1;
                        % usefacecolor=usefacecolor*2;
                        % usefacecolor(usefacecolor>64)=64;

                        usefacecolor = usefacecolor + 1;
                        usefacecolor = ind2rgb(round(usefacecolor),cmap);
                    else
                        usefacecolor=nan; % won't draw the point then.
                    end
                else
                    if options.d3.hlactivecontacts && ismember(cntct,find(elstruct.activecontacts{side})) % make active red contact without transparency
                        usefacecolor=[0.9,0.1,0.1];
                    else
                        if isfield(elstruct,'group')
                            usefacecolor=elstruct.groupcolors(elstruct.group,:);
                        else
                            usefacecolor=[1,1,1];
                        end
                    end
                end

                if ~any(isnan(usefacecolor))
                    set(0,'CurrentFigure',resultfig);
                    if ~shifthalfup
                        if strcmpi(options.prefs.d3.pointcloudstyle, '3d')
                            elrender(pcnt)=ea_plotsphere(coords_mm{side}(cntct,:), ms/10, usefacecolor, useedgecolor);
                            elrender(pcnt).Tag = 'PointCloud';
                        else
                            elrender(pcnt)=plot3(coords_mm{side}(cntct,1),coords_mm{side}(cntct,2),coords_mm{side}(cntct,3),'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                        end
                    else
                        if strcmpi(options.prefs.d3.pointcloudstyle, '3d')
                            elrender(pcnt)=ea_plotsphere(mean([coords_mm{side}(cntct,:);coords_mm{side}(cntct+1,:)],1), ms/10, usefacecolor, useedgecolor);
                            elrender(pcnt).Tag = 'PointCloud';
                        else
                            elrender(pcnt)=plot3(mean([coords_mm{side}(cntct,1),coords_mm{side}(cntct+1,1)]),...
                                mean([coords_mm{side}(cntct,2),coords_mm{side}(cntct+1,2)]),...
                                mean([coords_mm{side}(cntct,3),coords_mm{side}(cntct+1,3)]),...
                                'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                        end
                    end
                    drawnow;
                    pcnt=pcnt+1;
                else

                end
                hold on
            end
        end
    end
end

if ~exist('elrender','var')
    elrender=nan;
end
