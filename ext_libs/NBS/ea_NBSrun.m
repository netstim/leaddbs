function ea_NBSrun(varargin)
%NBSrun Reads user inputs, performs basic error checking and then performs 
%operations specified by user inputs.
%
%   NBSrun(UI) operates on the user inputs in the structure UI.
%
%   NBSrun(UI,S) attempts to write progress to the graphical user interface
%   with handles S, otherwise writes progress to the screen. S is created
%   with NBS. 
%
%   All fields in a UI strucutre must be a string. 
%   A UI structure contains the following fields:
%   UI.method.ui:       'Run NBS' | 'Run FDR'
%                       Perform NBS or FDR?
%
%   UI.test.ui:         'One Sample' | 't-test' | 'F-test'
%                       Statistical test to perform. 
%                       See also NBSglm
%
%   UI.size.ui:         'Extent' | 'Intensity'
%                       Use intensity or extent to measure size of a 
%                       network component?
%                       [optional if UI.method.ui=='Run FDR']
%                       See also NBSstats
%
%   UI.thresh.ui:       Scalar 
%                       Primary test statistic threshold. 
%                       [optional if UI.method.ui=='Run FDR']
%                       See also NBSstats
%
%   UI.perms.ui:        Scalar integer
%                       Number of permutations. 
%                       See also NBSglm
%
%   UI.alpha.ui:        Scalar 
%                       Significance (alpha threshold). 
%                       See also NBSstats
%   
%   UI.contrast.ui      1 x p numeric array specifying contrast, where p 
%                       is the number of independent variables in the GLM.
%                       Must be specified as a valid Matlab expression
%                       for a 1 x p array
%                       See also NBSglm
%
%   UI.design.ui        n x p numeric array specifying a design matrix, 
%                       including a column of ones if necessary. p is the 
%                       number of independent variables in the GLM, n is 
%                       the number of observations. 
%                       Can be specified either as a:
%                       1. Valid Matlab expression for an n x p array
%                       2. Text file containing numeric data arranged into
%                          n rows and p columns
%                       3. A binary Matlab file (.mat) storing an n x p
%                          numeric array
%                       See also NBSglm
%
%   UI.matrices.ui      N x N numeric array specifying a symmetric 
%                       connectivity matrix for each of M observations 
%                       (e.g. subjects), where N is the number of nodes. 
%                       Can be specified either as a:
%                       1. Valid Matlab expression for an N x N x M array
%                       2. A total of M seperate text files stored in a 
%                          common directory, where each text file contains 
%                          numeric data arranged into N rows and N columns. 
%                          Specify only one such text file and the others 
%                          within the same directory will be identified 
%                          automatically. 
%                       3. A binary Matlab file (.mat) storing an N x N x M
%                          numeric array
%                        
%   UI.exchange.ui:     n x 1 numeric array specifying exchange blocks
%                       to constrain permutation for a repeated measures 
%                       design, where n is the number of observations in 
%                       the GLM 
%                       [optional]
%                       Can be specified either as a:
%                       1. Valid Matlab expression for an n x 1 array
%                       2. Text file containing numeric data arranged into
%                          n rows 
%                       3. A binary Matlab file (.mat) storing an n x 1 
%                          numeric array
%                       See also NBSglm
%
%   UI.node_coor.ui:    N x 3 numeric array specifying node coordinates 
%                       in MNI space, where N is the number of nodes
%                       [optional]
%                       Can be specified either as a:
%                       1. Valid Matlab expression for an N x 3 array
%                       2. Text file containing numeric data arranged into
%                          a N rows and 3 columns
%                       3. A binary Matlab file (.mat) storing an N x 3 
%                          numeric array
%                       See also NBSview
%
%   UI.node_label.ui:   N x 1 cell array of strings providing node labels, 
%                       where N is the number of nodes 
%                       [optional]
%                       Can be specified either as a:
%                       1. Valid Matlab expression for an N x 1 cell array 
%                          of strings
%                       2. Text file containing data arranged into N rows 
%                       3. A binary Matlab file (.mat) storing an N x 1 
%                          cell array of strings
%                       See also NBSview    
%
%   UI structure corresponding to the example data provided:
%         UI.method.ui='Run NBS'; 
%         UI.test.ui='t-test';
%         UI.size.ui='Extent';
%         UI.thresh.ui='3.1';
%         UI.perms.ui='5000';
%         UI.alpha.ui='0.05';
%         UI.contrast.ui='[-1,1]'; 
%         UI.design.ui='SchizophreniaExample\designMatrix.txt';
%         UI.exchange.ui=''; 
%         UI.matrices.ui='SchizophreniaExample\matrices\subject01.txt';
%         UI.node_coor.ui='SchizophreniaExample\COG.txt';                         
%         UI.node_label.ui='SchizophreniaExample\nodeLabels.txt';
%
%   Remarks:
%       This function can be used as a command line version of NBS: 
%           1. Specify inputs in structure UI 
%           2. Run NBSrun(UI)
%           3. At completion, results stored in structure nbs
%              Type 'global nbs' before attempting to access the structure.
%
%   See also NBS
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   azalesky@unimelb.edu.au


%Declare the nbs structure global to avoid passing between NBS and NBSrun
global nbs

%Don't precompute randomizations if the number of test statistics populates 
%a matrix with more elements than Limit. Slows down computation but saves 
%memory.   
Limit=10^8/3;

%Waitbar position in figure
WaitbarPos=[0.69 0.021 0.05 0.21];

%User inputs
UI=varargin{1}; 

%Handles to GUI objects to enable progress updates to be written to GUI
if nargin==2
    S=varargin{2};
    try set(S.OUT.ls,'string',[]); catch; end
end

%Assume UI is valid to begin with
%Can be set to zero after reading UI or performing error checking
UI.method.ok=1;
UI.design.ok=1;
UI.contrast.ok=1;
UI.thresh.ok=1;
UI.test.ok=1;
UI.matrices.ok=1;
UI.node_coor.ok=1;
UI.node_label.ok=1;
UI.perms.ok=1;
UI.alpha.ok=1;
UI.size.ok=1;
UI.exchange.ok=1;
    
%Read UI and assign to appropriate structure
%Connectivity matrices
if exist(fileparts(UI.matrices.ui),'dir')
    [nbs.GLM.y,UI.matrices.ok,DIMS]=read_matrices([fileparts(UI.matrices.ui),filesep]);
    if ~UI.matrices.ok
        [nbs.GLM.y,UI.matrices.ok,DIMS]=read_matrices(UI.matrices.ui);
    end
else
    [nbs.GLM.y,UI.matrices.ok,DIMS]=read_matrices(UI.matrices.ui);
end
%Design matrix
[nbs.GLM.X,UI.design.ok,DIMS]=read_design(UI.design.ui,DIMS); 
%Contrast

[nbs.GLM.contrast,UI.contrast.ok]=read_contrast(UI.contrast.ui,DIMS);
%Node coordinates [optional, but mandatory for NBSview]
[tmp,UI.node_coor.ok]=read_node_coor(UI.node_coor.ui,DIMS); 
if UI.node_coor.ok
   nbs.NBS.node_coor=tmp; 
elseif isfield(nbs,'NBS')
    if isfield(nbs.NBS,'node_coor')
        %Delete field if was not able to read
        %Covers the situation when was able to read in previous run, but not in
        %the current run
        nbs.NBS=rmfield(nbs.NBS,'node_coor');
    end
end
%Node labels [optional]
[tmp,UI.node_label.ok]=read_node_label(UI.node_label.ui,DIMS); 
if UI.node_label.ok
   nbs.NBS.node_label=tmp;
elseif isfield(nbs,'NBS')
    if isfield(nbs.NBS,'node_label')
        nbs.NBS=rmfield(nbs.NBS,'node_label');
    end
end
%Exchange blocks for permutation [optional]
[tmp,UI.exchange.ok]=read_exchange(UI.exchange.ui,DIMS);
if UI.exchange.ok
    nbs.GLM.exchange=tmp; 
elseif isfield(nbs,'GLM')
    if isfield(nbs.GLM,'exchange')
        nbs.GLM=rmfield(nbs.GLM,'exchange');
    end
end
%Test statistic
try nbs.GLM.test=UI.test.ui; 
    if strcmp(nbs.GLM.test,'One Sample')
        nbs.GLM.test='onesample';
    elseif strcmp(nbs.GLM.test,'t-test')
        nbs.GLM.test='ttest';
    elseif strcmp(nbs.GLM.test,'F-test')
        nbs.GLM.test='ftest';
    end
catch; UI.test.ok=0; end
%Number of permutations
try nbs.GLM.perms=str2num(UI.perms.ui); catch; UI.perms.ok=0; end 
try if ~isnumeric(nbs.GLM.perms) || ~(nbs.GLM.perms>0)
    UI.perms.ok=0; end
catch; UI.perms.ok=0; end
%Test statistic threshold
try nbs.STATS.thresh=str2num(UI.thresh.ui); catch; UI.thresh.ok=0; end
try if ~isnumeric(nbs.STATS.thresh) || ~(nbs.STATS.thresh>0)
    UI.thresh.ok=0; end
catch; UI.thresh.ok=0; end
%Corrected p-value threshold
try nbs.STATS.alpha=str2num(UI.alpha.ui); catch; UI.alpha.ok=0; end 
try if ~isnumeric(nbs.STATS.alpha) || ~(nbs.STATS.alpha>0)
    UI.alpha.ok=0; end
catch; UI.alpha.ok=0; end
%Component size 
try nbs.STATS.size=UI.size.ui; catch; UI.ok=0; end 
%Number of nodes
nbs.STATS.N=DIMS.nodes;

%Do error checking on user inputs
[msg,stop]=errorcheck(UI,DIMS,S);
%Attempt to print result of error checking to listbox. If this fails, print
% to screen
try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[msg;tmp]); drawnow;
catch
    for i=1:length(msg)
        fprintf('%s\n',msg{i}); 
    end
end
%Do not proceed with computation if mandatory user inputs are missing or
%cannot be read
if stop 
    return
end

%Test statistic
%Only compute the test statistic if it has not been previously computed or
%the relevant inputs have changed 
repeat=1;
if isfield(nbs,'UI') 
    %There has been a previous run
    %Note that after a successful run, UI is copied to nbs.UI
    %Check what inputs have changed since previous run
    if ~repeat_glm(UI,nbs.UI)
        repeat=0;
    end
end
if repeat
    %Check whether the number of elements in the test_stat matrix exceeds
    %Limit
    if (DIMS.nodes*(DIMS.nodes-1)/2)*(nbs.GLM.perms)<Limit
        %Precompute if the number of elements in test_stat is less than
        %Limit
        str='Randomizing data...';
        try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[{str};tmp]); drawnow;
        catch;  fprintf([str,'\n']); end 
        %Present a waitbar on the GUI showing progress of the randomization process
        %Parent of the waitbar is the figure          
        try S.OUT.waitbar=uiwaitbar(WaitbarPos,S.fh); drawnow;  
        catch; S.OUT.waitbar=[]; end
        nbs.STATS.test_stat=zeros(nbs.GLM.perms+1,DIMS.nodes*(DIMS.nodes-1)/2); 
        nbs.STATS.test_stat=ea_NBSglm(nbs.GLM,S.OUT.waitbar);  
        delete(S.OUT.waitbar); 
    else
        %Too big to precompute 
        str='Too many randomizations to precompute...';
        try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[{str};tmp]); drawnow;
        catch;  fprintf([str,'\n']); end 
        nbs.STATS.test_stat=[]; 
    end
end

%Do NBS
if strcmp(UI.method.ui,'Run NBS')
    str='Computing network components...';
    try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[{str};tmp]); drawnow;
    catch;  fprintf([str,'\n']); end 
    try [nbs.NBS.n,nbs.NBS.con_mat,nbs.NBS.pval]=ea_NBSstats(nbs.STATS,S.OUT.ls,nbs.GLM);
    catch; [nbs.NBS.n,nbs.NBS.con_mat,nbs.NBS.pval]=ea_NBSstats(nbs.STATS,-1,nbs.GLM); end
%Do FDR
elseif strcmp(UI.method.ui,'Run FDR')
    str='False Discovery Rate...';
    try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[{str};tmp]); drawnow;
    catch;  fprintf([str,'\n']); end 
    %Show waitbar if test statistics have not been precomputed
    if isempty(nbs.STATS.test_stat)
        try S.OUT.waitbar=uiwaitbar(WaitbarPos,S.fh); drawnow; 
        catch; S.OUT.waitbar=[]; end
        [nbs.NBS.n,nbs.NBS.con_mat,nbs.NBS.pval]=ea_NBSfdr(nbs.STATS,S.OUT.waitbar,nbs.GLM);
        delete(S.OUT.waitbar); 
    else
        [nbs.NBS.n,nbs.NBS.con_mat,nbs.NBS.pval]=ea_NBSfdr(nbs.STATS);
    end
end
    
%Update the UI in the nbs structure to the UI that has just been used for
%the current run
nbs.UI=UI; 

%Copy test statistics to NBS strucutre so that they can be displayed with
%each link
test_stat=zeros(2,DIMS.nodes*(DIMS.nodes-1)/2);
if isempty(nbs.STATS.test_stat)
    K=nbs.GLM.perms;
    %Temporarily set to 1 to save computation
    nbs.GLM.perms=1;
    test_stat=ea_NBSglm(nbs.GLM);
    %Set back to original value
    nbs.GLM.perms=K;
else
    test_stat=nbs.STATS.test_stat(1,:);
end
ind_upper=find(triu(ones(DIMS.nodes,DIMS.nodes),1));
nbs.NBS.test_stat=zeros(nbs.STATS.N,nbs.STATS.N);
nbs.NBS.test_stat(ind_upper)=test_stat(1,:); 
nbs.NBS.test_stat=nbs.NBS.test_stat+nbs.NBS.test_stat';

%Display significant results with NBSview only if node coordinates provided
if nbs.NBS.n>0 && UI.node_coor.ok
    ea_NBSview(nbs.NBS); str=[];
elseif nbs.NBS.n>0 && ~UI.node_coor.ok
    str='Significant result - specify Node Coordinates to view';
    
else
    str='No significant result'; 
end
if ~isempty(str)
    try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[{str};tmp]); drawnow;
    catch;  fprintf([str,'\n']); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read connectivity matrices and vectorize the upper triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,ok,DIMS]=read_matrices(Name)
    ok=1;
    data=ea_readUI(Name);
    if ~isempty(data)
        [nr,nc,ns]=size(data);
        if nr==nc && ns>0 && ~iscell(data) && isnumeric(data)
            ind_upper=find(triu(ones(nr,nr),1));
            y=zeros(ns,length(ind_upper));
            %Collapse matrices
            for i=1:ns
                tmp=data(:,:,i);
                y(i,:)=tmp(ind_upper);
            end
        elseif iscell(data)
            [nr,nc]=size(data{1});
            ns=length(data);
            if nr==nc && ns>0
                ind_upper=find(triu(ones(nr,nr),1));
                y=zeros(ns,length(ind_upper));
                %Collapse matrices
                for i=1:ns
                    [nrr,ncc]=size(data{i});
                    if nrr==nr && ncc==nc && isnumeric(data{i})
                        y(i,:)=data{i}(ind_upper);
                    else
                        ok=0; y=[]; 
                        break
                    end
                end
            else
                ok=0; y=[];
            end
        end
    else
        ok=0; y=[];
    end
    if ok==1
        %Number of nodes
        DIMS.nodes=nr;
        %Number of matrices
        DIMS.observations=ns;
    else
        %Number of nodes
        DIMS.nodes=0;
        %Number of matrices
        DIMS.observations=0;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read design matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,ok,DIMS]=read_design(Name,DIMS)
ok=1;
data=ea_readUI(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.observations && nc>0 && ns==1 && isnumeric(data) 
        X=data; 
    else
        ok=0; X=[];
    end
else
    ok=0; X=[];
end
clear data
if ok==1
    %Number of predictors
    DIMS.predictors=nc;
else
    DIMS.predictors=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read contrast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [contrast,ok]=read_contrast(Name,DIMS)
ok=1; 
data=ea_readUI(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data); 
    if nr==1 && nc==DIMS.predictors && ns==1 && isnumeric(data) 
        contrast=data; 
    else
        ok=0; contrast=[];
    end
else
    ok=0; contrast=[];
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read node coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_coor,ok]=read_node_coor(Name,DIMS)
ok=1;
data=ea_readUI(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.nodes && nc==3 && ns==1 && isnumeric(data)
        node_coor=data; 
    else
        ok=0; node_coor=[];
    end
else
    ok=0; node_coor=[];
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read node labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_label,ok]=read_node_label(Name,DIMS)
ok=1;
data=ea_readUI(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.nodes && nc==1 && ns==1
        node_label=data; 
    else
        ok=0; node_label=[];
    end
else
    ok=0; node_label=[]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read permutation exchange blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [exchange,ok]=read_exchange(Name,DIMS)
ok=1;
data=ea_readUI(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.observations && nc==1 && ns==1
        exchange=data'; 
    else
        if nc==DIMS.observations && nr==1 && ns==1
            exchange=data;
        else
        ok=0; exchange=[];
        end
    end
else
    ok=0; exchange=[];
end
    %Set up exchange blocks
    blks=unique(exchange); 
    %Number of blocks
    n_blks=length(blks);
    %Number of observations per block
    sz_blk=length(exchange)/n_blks;
    if rem(sz_blk,1)>0
        ok=0;
        exchange=[];
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%There is not need to solve the GLM again if all of the GLM inputs have
%remained the same
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function repeat=repeat_glm(UI,UIprev)
    repeat=~all([strcmp(UI.matrices.ui,UIprev.matrices.ui),... 
                 strcmp(UI.design.ui,UIprev.design.ui),... 
                 strcmp(UI.contrast.ui,UIprev.contrast.ui),...
                 strcmp(UI.test.ui,UIprev.test.ui),... 
                 strcmp(UI.perms.ui,UIprev.perms.ui),...
                 strcmp(UI.exchange.ui,UIprev.exchange.ui)]);

function [msg,stop]=errorcheck(UI,DIMS,S)
    stop=1;
    %Mandatroy UI
    %UI.method.ok %no need to check
    if ~UI.matrices.ok
        msg={'Stop: Connectivity Matrices not found or inconsistent'};
        try set(S.DATA.matrices.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.design.ok
        keyboard
        msg={'Stop: Design Matrix not found or inconsistent'};
        try set(S.STATS.design.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.contrast.ok 
        msg={'Stop: Contrast not found or inconsistent'};
        try set(S.STATS.contrast.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.thresh.ok
        msg={'Stop: Threshold not found or inconsistent'};
        try set(S.STATS.thresh.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.test.ok 
        msg={'Stop: Statistical Test not found or inconsistent'};
        try set(S.STATS.test.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.perms.ok
        msg={'Stop: Permutations not found or inconsistent'};
        try set(S.ADV.perms.text,'ForegroundColor','red');
        catch; end
        return;
    end        
    if ~UI.alpha.ok
        msg={'Stop: Significance not found or inconsistent'};
        try set(S.ADV.alpha.text,'ForegroundColor','red');
        catch; end
        return;
    end
    if ~UI.size.ok
        msg={'Stop: Component Size not found or inconsistent'};
        try set(S.ADV.size.text,'ForegroundColor','red');
        catch; end
        return;
    end
    stop=0;
    
    msg=[{sprintf('Nodes: %d',DIMS.nodes)};...
             {sprintf('Observations: %d',DIMS.observations)};...
             {sprintf('Predictors: %d',DIMS.predictors)}]; 
    
    %Optional, but mandatory for NBSview
    if ~UI.node_coor.ok
        msg=[msg;{'Node Coordinates: No'}];
        try set(S.DATA.node_coor.text,'ForegroundColor','red');
        catch; end
    end
    
    %Optional
    if ~UI.exchange.ok
        msg=[msg;{'Exchange Blocks: No'}];
        try set(S.ADV.exchange.text,'ForegroundColor','red');
        catch; end
    else
        msg=[msg;{'Exchange Blocks: Yes'}];
    end
    if ~UI.node_label.ok
        msg=[msg;{'Node Labels: No'}];
        try set(S.DATA.node_label.text,'ForegroundColor','red');
        catch; end
    end