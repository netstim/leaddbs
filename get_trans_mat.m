function [trans_mat,elfv,xg,yg,zg] = get_trans_mat(electrode,electrode_patient,grid_vec,side,electrode_native)
% inputs:
% electrode: is the electrode model
% electrode_patient: is the location of the patient's electrode in mni
% space
% grid_vec: is the standard grid points, on whihc the standard efield is
% calculated
% side: 1 is right side of the brain and 2 is left side.

% output:
% transmat: is the transformation matrix from the electrode model to
% patietn electrode in mni space
% elfv: are the faces and vertices of the patient electrode in mni space
% xg, yg, zg, are the transformed coordinated of the electric field in mni
% space


%% Calculate the transformation matrix from the stamdard space to patient mni space
% 
A = [electrode.head_position,1;
electrode.tail_position,1
electrode.x_position,1
electrode.y_position,1]; % points in model

B = [electrode_patient.markers(side).head,1;
electrode_patient.markers(side).tail,1;
electrode_patient.markers(side).x,1;
electrode_patient.markers(side).y,1];

Y = mldivide(A,B); Y = Y';
trans_mat = Y; 

%% move the grid to patient mni space
[Y1,X1,Z1] = meshgrid(grid_vec{1},grid_vec{2},grid_vec{3});

temp = [X1(:),Y1(:),Z1(:)];
ind_mat4 = trans_mat*[temp,ones(size(temp,1),1)]';
ind_mat4 = ind_mat4(1:3,:)';

sz = size(X1);
xg = reshape(ind_mat4(:,1),sz);
yg = reshape(ind_mat4(:,2),sz);
zg = reshape(ind_mat4(:,3),sz);

%% move the electrode to patient mni space

cnt = 1;
    for con = 1:length(electrode.contacts)

        electrode.contacts(con).vertices = trans_mat*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices = electrode.contacts(con).vertices(1:3,:)';

        elfv(cnt).faces = electrode.contacts(con).faces;
        elfv(cnt).vertices = electrode.contacts(con).vertices;
        tissuetype(cnt) = 3;
        t = surfinterior(elfv(cnt).vertices,elfv(cnt).faces);
        cnt=cnt+1;
    end


    for ins = 1:length(electrode.insulation)

        electrode.insulation(ins).vertices = trans_mat*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
        electrode.insulation(ins).vertices = electrode.insulation(ins).vertices(1:3,:)';


        elfv(cnt).faces = electrode.insulation(ins).faces;
        elfv(cnt).vertices = electrode.insulation(ins).vertices;
        t = surfinterior(elfv(cnt).vertices,elfv(cnt).faces);
        cnt = cnt+1;
    end

end