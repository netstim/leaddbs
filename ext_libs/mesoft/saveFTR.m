function [ftr lens] = saveFTR(finame,datastruct,I,vox,info,saveit)

if nargin == 5,
    saveit = true;
end;

cc = {};
csc = {};
cD = {};
lens = [];
 
if isempty(vox),
    warning('Empty field vox: guessing');
    vox = [2 2 2];
end;


reparam_step = 0.9*min(vox); %mm

P = datastruct.state;

labels = zeros(length(I),1);

if ~isempty(P)

    trof = inv(diag(vox));

    lens = [];
    csc = cell(length(I),1);
    cc = cell(length(I),1);
    lens = zeros(length(I),1);
    cnt = 1;
    for k = 1:length(I)
        if ~isempty(I{k}),
            i = I{k}(1,:)+1;
            cD{cnt} = single(P(11:15,i));
            label(k) = P(7,i(1));
            fibpoints = P(1:3,i);
            if size(fibpoints,2) == 1,
                ll = P(14,i);
                fibpoints = [P(1:3,i)+ll.*P(4:6,i) P(1:3,i)-ll.*P(4:6,i)];
                cD{cnt} = single([P(11:15,i) P(11:15,i)]);
            end;

            [fibpoints leng] =  reparametrize_arclen(single(cat(1,fibpoints,cD{cnt})),double(reparam_step));
            cD{cnt} = fibpoints(4:end,:);
            fibpoints = fibpoints(1:3,:);
            
            leng = size(fibpoints,2);
            
            fibpoints = single(trof * fibpoints);
            fibpoints = fibpoints + 1; %% convert to 1-based coordinates used in fibertools
                                    
            csc{cnt} = fibpoints';
            
            lens(k) = leng;
            
            cc{cnt} = cnt;
            cnt = cnt + 1;
        end;
    end;
end;

thestate.P = P;
thestate.vf = datastruct.vfmap;
thestate.S2 = datastruct.S2;
thestate.meansignal = datastruct.meansignal;
thestate.b0avg = datastruct.b0avg;

ftr.connectCell = cc;
ftr.curveSegCell = csc;
ftr.curveD = cD;
ftr.posSegCell = {};
ftr.vox = vox;
ftr.patient = info.name;
ftr.user = thestate;
ftr.dtdType = 'unused';
ftr.fiber = {};
ftr.algoName = 'mesoGT_V1';
ftr.trackParam = info;
ftr.trackDate = date;
ftr.logData = [];
ftr.labels = label;
if isstruct(info.edges),
    E = [info.edges.hdr.hist.srow_x ; info.edges.hdr.hist.srow_y ; info.edges.hdr.hist.srow_z]; E = [E ; zeros(1,4)]; E(4,4) = 1;
    ftr.hMatrixNifti = E;
    Q = diag([-1 -1 1 1]);    
    ftr.hMatrix =    Q*E; % differnet scanner coordinates in mrstructs....

    
else
    ftr.hMatrix = info.edges;
end;
ftr.version = 'V1.1';
if saveit,
    save(finame,'-struct','ftr');
end;
