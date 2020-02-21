function data=ea_readUI(UI)
%readUI Attempt to read a user input. 
%
%   DATA = readUI(UI) reads a user input described in UI and returns the
%   relevant data. An empty array is returned if UI could not be read.
%
%   Remark:
%       Attempts to read UI are made in the following order:
%           1. Binary matlab file containing one variable
%           2. Text file containing numeric data
%           3. Text file containing a column of strings
%           4. Valid Matlab expression  
%           5. Multiple text files in a common directory. The numeric data
%              in the i'th text file is placed in data{i}. 
%
%   See also NBSrun
%
%   azalesky@unimelb.edu.au
data=[]; 
[path,name,ext]=fileparts(UI);
if strcmp(ext,'.mat')
    %Assume data stored as a Matlab file
    try tmp=load(UI);
        f=fieldnames(tmp);
        if length(f)>1
            fprintf('%s contains multiple variables',[name,ext]);
        end
        try data=eval(['tmp.',f{1}]); clear tmp f; catch end
        catch end
elseif ~isempty(name)
    %Assume data stored in an ASCII file
    %First try numeric data, then try a column of strings
    try data=dlmread(UI); catch 
        try data=textread(UI,'%s'); catch 
           %Try to evaluate as a Matlab expression
            try data=eval(UI); catch end
        end
    end 
elseif isdir(path)
    %Assume data stored as ASCII files in the directory specified by path
    try tmp=dir(UI);
        cnt=0;
        for i=1:length(tmp)
            try mat=dlmread([path filesep tmp(i).name]);
               cnt=cnt+1;
               data{cnt}=mat; clear mat; 
            catch end
        end
    catch end
end