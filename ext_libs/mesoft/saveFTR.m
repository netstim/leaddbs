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


smoothker = permute(fspecial('gaussian',[5 1],2),[2 3 1]);

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
            
           % cD{cnt} = single(P(11:15,i));
           % cD{cnt} = single(datastruct.segparams(i,:)');
            
            label(k) = P(7,i(1));
            fibpoints = P(1:3,i);
            if size(fibpoints,2) == 1,
                ll = P(14,i)*1;
                fibpoints = [floor(P(1:3,i)/2)*2+ll.*P(4:6,i) floor(P(1:3,i)/2)*2-ll.*P(4:6,i)];
%                cD{cnt} = single([P(11:15,i) P(11:15,i)]);
                cD{cnt} = single([datastruct.segparams(i,:)' datastruct.segparams(i,:)']);
            end;

            
            
            [fibpoints leng] =  reparametrize_arclen(single(fibpoints),double(reparam_step));
%            [fibpoints leng] =  reparametrize_arclen(single(cat(1,fibpoints,cD{cnt})),double(reparam_step));
            
%             
%             cD{cnt} = reshape(fibpoints(4:end,:),[size(datastruct.segparams,2) size(datastruct.segparams,3) size(fibpoints,2)]);
%                         
%             cD{cnt} = imfilter(cD{cnt},smoothker,'replicate');
%             
%             M = permute(cD{cnt},[3 1 2]);           
%             M = M(:,:,2:end) .* repmat(1./(eps+M(:,:,1)),[1 1 size(M,3)-1]);
% 
%             modelstruc = evalin('base','modelstruc');
% 
%             est = (modelstruc.apply(M))';
%             cD{cnt}= [est(1,:) ;est(2,:)+est(3,:); est(3,:); est(4,:)];
% 
            
            
% 
% 
% lmax = 4;      
%       
% 
% M = permute(cD{cnt},[3 1 2]);           
% M = M(:,1:(lmax/2+1),:);
% 
% 
% M = M(:,:,2:end) .* repmat(1./(eps+M(:,:,1)),[1 1 size(M,3)-1]);
%             
%              
% 
% tmp = squeeze(M(:,2:end,1:end-1)./repmat(sum(M(:,2:end,:),3),[1 1 size(M,3)-1]));
% M = [squeeze(M(:,1,:)) tmp(:,:) ];
% 
% M = M(:,:);
% maxN = size(M,2);
% 
% 
% Q = [M(:,1)*0+1 ];
% M = (M);
% for kk = 1:maxN,
%     Q = cat(2,Q,M(:,kk));
%     for j = kk:maxN,
%         Q = cat(2,Q,M(:,kk).*M(:,j));
%         for r = j:maxN,
%             Q = cat(2,Q,M(:,kk).*M(:,j).*M(:,r));
%         end
%     end;        
% end;
% 
% alpha = evalin('base','alpha');
% %%
% est =  (Q*alpha)';



            
%%            
            
            
            
            
            
            
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
thestate.stderr = datastruct.stderr;

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
ftr.trackParam.bvals = datastruct.bvals;
ftr.trackDate = date;
ftr.logData = [];
%ftr.labels = label;
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
    save(finame,'-struct','ftr','-v7.3');
end;
