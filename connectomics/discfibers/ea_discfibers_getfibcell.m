function [fibcell,fibsin,XYZmm,niivx,valsmm]=ea_discfibers_getfibcell(obj,cfile)
if exist(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo','baseinfo.mat']),'file')
    load(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo','baseinfo.mat']));
    if isfield(obj.results.(ea_conn2connid(obj.connectome)),'fibcell')
        fibcell = obj.results.(ea_conn2connid(obj.connectome)).fibcell;
        fibsin=fibcell2fibsin(fibcell);
        return
    end
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

if isfield(obj.results.(ea_conn2connid(obj.connectome)),'fibcell')
    fibcell = obj.results.(ea_conn2connid(obj.connectome)).fibcell;
    fibsin=fibcell2fibsin(fibcell);
else
    [fibsin,XYZmm,nii,valsmm]=ea_discfibers_genroilist_connfibers(fibers, allroilist);
    for vta=1:numel(nii)
        niivx(vta,:)=nii{vta}.voxsize;
    end
    [fibcell,fibsin]=fibsin2fibcell(fibsin);
end
obj.results.(ea_conn2connid(obj.connectome)).fibcell=fibcell;
niivx=mean(niivx);
ea_mkdir(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo'));
save(fullfile(fileparts(obj.leadgroup),'disctracts','baseinfo','baseinfo.mat'),...
    'XYZmm','niivx','valsmm','-v7.3');




function [fibcell,fibsin]=fibsin2fibcell(fibsin)
% now color fibsin based on predictive value of improvement
ea_dispt('');


% Reformat to cell:
[~,fibiaxfirst]=unique(fibsin(:,4),'first');
[~,fibiaxlast]=unique(fibsin(:,4),'last');
fiblen = fibiaxlast - fibiaxfirst + 1;
fibcell = mat2cell(fibsin(:,1:3),fiblen);
% repair fibsin to be incrementing from 1 to x:
for f=1:length(fibcell)
    fibsin(fibiaxfirst(f):fibiaxlast(f),4)=f;
end

function fibsin=fibcell2fibsin(fibcell)
fibsin=cell2mat(fibcell);
fibsin=[fibsin,zeros(size(fibsin,1),1)];
cnt=1;
for fib=1:length(fibcell)
    thisfiblen=length(fibcell{fib});
    fibsin(cnt:cnt+thisfiblen-1,4)=fib;
    cnt=cnt+thisfiblen;
end
