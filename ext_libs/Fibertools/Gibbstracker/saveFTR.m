function [ftr lens] = saveFTR(finame,P,I,offset,info,saveit)

if nargin == 5,
    saveit = true;
end;

cc = {};
csc = {};
lens = [];
 
if isempty(offset),
    warning('Empty field offset: guessing');
    offset = eye(4);
    offset(1,1) = 2;     offset(2,2) = 2;     offset(3,3) = 2; 
end;

vox = [offset(1,1) offset(2,2) offset(3,3)];


reparam_step = 0.9*min(vox); %mm

saveP = P;

if ~isempty(P)

    trof = inv(offset(1:3,1:3));

%     P(1:3,:) = trof*P(1:3,:);
%     P(4:6,:) = trof*P(4:6,:);
% 
%     P(1,:) = P(1,:) + offset(1,4);
%     P(2,:) = P(2,:) + offset(2,4);
%     P(3,:) = P(3,:) + offset(3,4);
 
    lens = [];
    csc = cell(length(I),1);
    cc = cell(length(I),1);
    lens = zeros(length(I),1);
    cnt = 1;
    for k = 1:length(I)
        if ~isempty(I{k}),
            i = I{k}(1,:)+1;
            %lens(k) = length(I{k});
            fibpoints = P(1:3,i);

            if size(fibpoints,2) == 1,
                ll = 2;
                fibpoints = [P(1:3,i)+ll*P(4:6,i) P(1:3,i)-ll*P(4:6,i)];
                
            end;
%             %% special treatment of endpoints
%             if P(9,i(1)) == -1, si = -1; else, si = 1; end;
%             fibpoints(1,:) = fibpoints(1,:) + si*P(8,i(1))*P(4:6,i(1))';
%             if P(9,i(end)) == -1, si = -1; else, si = 1; end;
%             fibpoints(end,:) = fibpoints(end,:) + si*P(8,i(end))*P(4:6,i(end))';
%             %%%
% 

            [fibpoints leng] = reparametrize_arclen(single(fibpoints),double(reparam_step));
            fibpoints = trof * fibpoints;
            fibpoints = fibpoints + 1; %% convert to 1-based coordinates used in fibertools
            
            fibpoints(1,:) = fibpoints(1,:) + offset(1,4);
            fibpoints(2,:) = fibpoints(2,:) + offset(2,4);
            fibpoints(3,:) = fibpoints(3,:) + offset(3,4);
                        
            csc{cnt} = fibpoints';
            
            lens(k) = leng;
            %csc{cnt} = fibpoints;
            
            cc{cnt} = cnt;
            cnt = cnt + 1;
        end;
    end;
end;

ftr.connectCell = cc;
ftr.curveSegCell = csc;
ftr.posSegCell = {};
ftr.vox = vox;
ftr.patient = info.name;
ftr.user = saveP;
ftr.dtdType = 'DTD';
ftr.fiber = {};
ftr.algoName = 'fiberGT_release2';
ftr.trackParam = info;
ftr.trackDate = date;
ftr.logData = [];
ftr.hMatrix = info.edges;
ftr.version = 'V1.1';
if saveit,
    save(finame,'-struct','ftr');
end;
