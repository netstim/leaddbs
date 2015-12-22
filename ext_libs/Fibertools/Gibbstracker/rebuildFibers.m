function rebuildFibers(ftrname,minlen,maxlen)

ftr = ftrstruct_read(ftrname);

display('building fibers')
ftr = recomputeFibs(ftr,minlen,maxlen);


[path name ext] = fileparts(ftrname);
newfile = fullfile(path,[name(1:end-4) '_' num2str(minlen) '_' num2str(maxlen) '_FTR.mat']);

ftrstruct_write(ftr, newfile);

function l = len(x)

l = sum(sqrt(sum((x(2:end,:)-x(1:end-1,:)).^2,2)));




function [ftr lens] = recomputeFibs(ftr,minlen,maxlen)

P = ftr.user;
I = BuildFibres(P,minlen);
   


cc = {};
csc = {};
lens = [];

vox = ftr.vox;

reparam_step = 0.9*min(vox); %mm

saveP = P;

if ~isempty(P)

   
    %lens = [];
    %csc = cell(length(I),1);
    %cc = cell(length(I),1);
    %lens = zeros(length(I),1);
    cnt = 1;
    
    trof = diag(1./vox);
    
    for k = 1:length(I)
        if ~isempty(I{k}),
            i = I{k}(1,:)+1;
            if length(I{k}) <= maxlen,
                fibpoints = P(1:3,i);


                [fibpoints leng] = reparametrize_arclen(single(fibpoints),double(reparam_step));
                fibpoints = trof * fibpoints;
                fibpoints = fibpoints + 1/2; %% convert to 1-based coordinates used in fibertools

                fibpoints(1,:) = fibpoints(1,:);
                fibpoints(2,:) = fibpoints(2,:);
                fibpoints(3,:) = fibpoints(3,:);

                csc{cnt} = fibpoints';

                lens(k) = leng;
                %csc{cnt} = fibpoints;

                cc{cnt} = cnt;
                cnt = cnt + 1;
            end;
        end;
    end;
end;

ftr.connectCell = cc;
ftr.curveSegCell = csc;











