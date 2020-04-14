function [fibcell,fibsin,XYZmm,niivx,valsmm]=ea_discfibers_getfibcell(obj,cfile)
if exist(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo','baseinfo.mat'),'file')
    D=load(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo','baseinfo.mat'));
    if isfield(obj.results,ea_conn2connid(obj.connectome)) && isfield(obj.results.(ea_conn2connid(obj.connectome)),'fibcell')
        fibcell = obj.results.(ea_conn2connid(obj.connectome)).fibcell;
        fibsin=fibcell2fibsin(fibcell);
        XYZmm=D.XYZmm;
        niivx=D.niivx;
        valsmm=D.valsmm;
        return
    end
else
    D=[];
end

if obj.M.ui.detached
    pthprefix=[fileparts(obj.leadgroup),filesep];
else
    pthprefix='';
end
allroilist=cell(length(obj.allpatients)*2,2);
cnt=1;
options.native = 0;

for sub=1:length(obj.allpatients) % all patients - for connected fibers selection ? and always flip
    allroilist(cnt,:)={[pthprefix,obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat','_efield','_right.nii'],[pthprefix,obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat','_efield','_left.nii']};
    cnt=cnt+1;
end

for sub=1:length(obj.allpatients) % all patients - for connected fibers selection ? and always flip
    ea_genflippedjointnii([pthprefix,obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat','_efield','_right.nii'],[pthprefix,obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat','_efield','_left.nii']);
    allroilist(cnt,:)={[pthprefix,obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'fl_','vat','_efield','_left.nii'],[pthprefix,obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'fl_','vat','_efield','_right.nii']};
    cnt=cnt+1;
end
allroilist={allroilist};

fibers=load(cfile);
fn=fieldnames(fibers);
try
    fibers=fibers.fibers;
catch
    fibers=fibers.(fn{1});
end


[fibsin,XYZmm,niivx,valsmm]=ea_discfibers_genroilist_connfibers(fibers, allroilist, D);
[fibcell,fibsin]=fibsin2fibcell(fibsin);

obj.results.(ea_conn2connid(obj.connectome)).fibcell=fibcell;
ea_mkdir(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo'));
save(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo','baseinfo.mat'),...
    'XYZmm','niivx','valsmm','-v7.3');


function [fibcell,fibsin]=fibsin2fibcell(fibsin)
% now color fibsin based on predictive value of improvement
ea_dispt('');
for side=1:2
    % Reformat to cell:
    [~,fibiaxfirst]=unique(fibsin{side}(:,4),'first');
    [~,fibiaxlast]=unique(fibsin{side}(:,4),'last');
    fiblen = fibiaxlast - fibiaxfirst + 1;
    fibcell{side} = mat2cell(fibsin{side}(:,1:3),fiblen);
    % repair fibsin to be incrementing from 1 to x:
    for f=1:length(fibcell{side})
        fibsin{side}(fibiaxfirst(f):fibiaxlast(f),4)=f;
    end
end


function fibsin=fibcell2fibsin(fibcell)
for side=1:2
    fibsin{side}=cell2mat(fibcell{side});
    fibsin{side}=[fibsin{side},zeros(size(fibsin{side},1),1)];
    cnt=1;
    for fib=1:length(fibcell{side})
        thisfiblen=length(fibcell{side}{fib});
        fibsin{side}(cnt:cnt+thisfiblen-1,4)=fib;
        cnt=cnt+thisfiblen;
    end
end
