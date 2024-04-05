function till_creategridforelectrode(reco,side,space,options)
% clearvars -except reco side space
% options = ea_getptopts('/media/netstim/E/Till/Testsession/derivatives/leaddbs/sub-BiDirect3');
%% load electrode model
try
    load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname,'.mat'],'electrode')
catch
    try
        options.elmodel = reco.props(side).elmodel;
        options = ea_resolve_elspec(options);
        load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname,'.mat'],'electrode')
    catch
        keyboard
    end
end
%% calculate transformation from template to space
A = [electrode.head_position,1;
    electrode.tail_position,1
    electrode.x_position,1
    electrode.y_position,1]; % points in model
B = [reco.(space).markers(side).head,1;
    reco.(space).markers(side).tail,1;
    reco.(space).markers(side).x,1;
    reco.(space).markers(side).y,1];
tmat = mldivide(A,B)';
clear A B
%% make dir for connectome
dirout = [options.subj.subjDir,filesep,'connectomes',filesep,'dMRI_MultiTract',filesep,'OSSDBSgrid',filesep];
if ~exist(dirout,'dir')
    mkdir(dirout)
end
%%
griddistance=.5;
gridnumber= 30;
gridorigin = (electrode.head_position+electrode.tail_position)/1;
%% parallelgrid
x=[gridorigin(1)-((gridnumber/1)*griddistance):griddistance:gridorigin(1)+((gridnumber/1)*griddistance)]';
y=[gridorigin(2)-((gridnumber/1)*griddistance):griddistance:gridorigin(2)+((gridnumber/1)*griddistance)]';
z=[gridorigin(3)-((gridnumber/1)*griddistance):griddistance:gridorigin(3)+((gridnumber/1)*griddistance)]';
[X,Y]=meshgrid(x,y);
fibers=[];
for i=1:numel(X)
    fibact=horzcat(repmat(X(i),numel(z),1),repmat(Y(i),numel(z),1),z,repmat(i,numel(z),1));
    fibers=vertcat(fibers,fibact);
end
idx = diff(find([true,diff(fibers(:,4)')~=0,true]))';

fiberstmp= tmat * vertcat(fibers(:,1:3)',ones(1,size(fibers,1)));
fibers(:,1:3)=fiberstmp(1:3,:)';
save([dirout 'Zgrid.mat'],'fibers','idx');
clear fibers idx fibersact X Y

%% Ygid
x=[gridorigin(1)-((gridnumber/1)*griddistance):griddistance:gridorigin(1)+((gridnumber/1)*griddistance)]';
y=[gridorigin(2)-((gridnumber/1)*griddistance):griddistance:gridorigin(2)+((gridnumber/1)*griddistance)]';
z=[gridorigin(3)-((gridnumber/1)*griddistance):griddistance:gridorigin(3)+((gridnumber/1)*griddistance)]';
[X,Y]=meshgrid(x,z);
fibers=[];
for i=1:numel(X)
    fibact=horzcat(repmat(X(i),numel(y),1),y,repmat(Y(i),numel(y),1),repmat(i,numel(y),1));
    fibers=vertcat(fibers,fibact);
end
idx = diff(find([true,diff(fibers(:,4)')~=0,true]))';

fiberstmp= tmat * vertcat(fibers(:,1:3)',ones(1,size(fibers,1)));
fibers(:,1:3)=fiberstmp(1:3,:)';
save([dirout 'Ygrid.mat'],'fibers','idx');
clear fibers idx fibersact X Y

%% Xgid
x=[gridorigin(1)-((gridnumber/1)*griddistance):griddistance:gridorigin(1)+((gridnumber/1)*griddistance)]';
y=[gridorigin(2)-((gridnumber/1)*griddistance):griddistance:gridorigin(2)+((gridnumber/1)*griddistance)]';
z=[gridorigin(3)-((gridnumber/1)*griddistance):griddistance:gridorigin(3)+((gridnumber/1)*griddistance)]';
[X,Y]=meshgrid(y,z);
fibers=[];
for i=1:numel(X)
    fibact=horzcat(x,repmat(X(i),numel(x),1),repmat(Y(i),numel(x),1),repmat(i,numel(x),1));
    fibers=vertcat(fibers,fibact);
end
idx = diff(find([true,diff(fibers(:,4)')~=0,true]))';

fiberstmp= tmat * vertcat(fibers(:,1:3)',ones(1,size(fibers,1)));
fibers(:,1:3)=fiberstmp(1:3,:)';
save([dirout 'Xgrid.mat'],'fibers','idx');
clear fibers idx fibersact X Y
end







% end