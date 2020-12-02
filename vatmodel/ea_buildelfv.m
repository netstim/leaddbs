function [elfv,tissuetype,Y,electrode]=ea_buildelfv(elspec,elstruct,side)


load([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname])

    A=[electrode.head_position,1;
        electrode.tail_position,1
        electrode.x_position,1
        electrode.y_position,1]; % points in model

    B=[elstruct.markers(side).head,1;
        elstruct.markers(side).tail,1;
        elstruct.markers(side).x,1;
        elstruct.markers(side).y,1];
    Y = mldivide(A,B); Y=Y';

    cnt=1;

    % overwrite head and tail of model with actual values for mesh generation lateron:
    electrode.head_position=B(1,1:3);
    electrode.tail_position=B(2,1:3);
    % add contacts to mesh
    for con=1:length(electrode.contacts)
        % do not rotate electrode for now.
        %electrode.meshel.con{con}.vertices=X*[electrode.meshel.con{con}.vertices,ones(size(electrode.meshel.con{con}.vertices,1),1)]';
        %electrode.meshel.con{con}.vertices=electrode.meshel.con{con}.vertices(1:3,:)';
        % only rotate contacts.
        electrode.contacts(con).vertices=Y*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';

        elfv(cnt).faces=electrode.contacts(con).faces;
        elfv(cnt).vertices=electrode.contacts(con).vertices;
        tissuetype(cnt)=3;
        t=surfinterior(elfv(cnt).vertices,elfv(cnt).faces);
        cnt=cnt+1;
    end

    % add insulation to mesh
    for ins=1:length(electrode.insulation)
        % do not rotate insulation for now
        %electrode.meshel.ins{ins}.vertices=X*[electrode.meshel.ins{ins}.vertices,ones(size(electrode.meshel.ins{ins}.vertices,1),1)]';
        %electrode.meshel.ins{ins}.vertices=electrode.meshel.ins{ins}.vertices(1:3,:)';
        % only rotate insulation.
        electrode.insulation(ins).vertices=Y*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
        electrode.insulation(ins).vertices=electrode.insulation(ins).vertices(1:3,:)';


        elfv(cnt).faces=electrode.insulation(ins).faces;
        elfv(cnt).vertices=electrode.insulation(ins).vertices;
        t=surfinterior(elfv(cnt).vertices,elfv(cnt).faces);

        tissuetype(cnt)=4;
        cnt=cnt+1;
    end
