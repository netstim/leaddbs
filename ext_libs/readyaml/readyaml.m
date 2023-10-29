function results = readyaml(filePath)
% readyaml() reads out YAML files and writes the data to a struct
%
% Syntax:  readyaml('myyamlfile.yaml')
%
% Inputs:
%    filepath        - string or char that gives the path to the yaml file
%
% Outputs:
%    results         - resulting struct from the yaml file
% 
% Author: Maarten J. Jongeneel, Ph.D. researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: info@maartenjongeneel.nl  
% October 2023; Last revision: 08-October-2023
%--------------------------------------------------------------------------

%We open the file and read out the yaml file to a data variable.
fid = fopen(filePath, 'r');
try
    data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
catch
    error(append("Cannot read the YAML file in path: ",filePath,", please check if the file exist and the path is correct."));
end
fclose(fid);

% remove empyy lines
data = deblank(data{1});
for ii = 1:length(data)
    %Comments begin with the number sign (#), can start anywhere on a line
    %and continue until the end of the line we remove them from the data
    if startsWith(strtrim(data{ii}),'#') 
        data{ii} = [];
    end
end

data(cellfun('isempty', data)) = [];
key = "data";
data{end+1}=''; %Add empty string to be able to write last row

%Here, we build the structlvl array which tells us the levels of a struct
for jj = 1:length(data); test{jj,1}=strrep(data{jj},'-',' ');end
for jj=1:length(data); structlvl(jj) = (length(test{jj})-length(strtrim(test{jj})))/2; end

dit = [0 diff(structlvl)];
ii =1;
while ~isempty(find(dit(ii:end)>1,1,'first')+ii-1) %as long as I can find an indentation that is bigger than 1
    strindex = find(dit(ii:end)>1,1,'first')+ii-1;
    endindex = strindex+find(structlvl(strindex:end)<=structlvl(strindex-1),1,"first")-2;
    index = (strindex:endindex);
    structlvl(index) = structlvl(index)-1;
    dit = [0 diff(structlvl)];
    if sum(dit(index)>1)==0
        ii = strindex; %This is just to update the end of struct
    end
end

%Now we have our stripped data and we know the structlevels, we go get the data
results = createStruct(data,1,key,structlvl);
end




function [results,ii,structlvl] = createStruct(data,ii,key,structlvl)
%If we just entered this function, we want to reset some values
lstnr = 0;
results = struct();
lst.bool= false;
lst.structlvl=[];
lst.data=[];

%This is the current object level. Important to track if we go down or back
%up one level in the struct
objlevel=structlvl(ii);

%Loop through all the lines
while ii<length(data)

    %Extract this line
    thisLine = data{ii};

    %Go back out of the function if we need to go back up one level
    if ii>1
        if structlvl(ii) < objlevel
            return;
        end
    end

    Ijustgotbackup = 0;
    %Check if we go down one level
    if ii>1 && structlvl(ii)>objlevel
        %Here we go down one level into the struct
        [value,ii,structlvl] = createStruct(data,ii,key,structlvl);

        %Here we just went back up, and we check if the value we got is
        %numeric (in case our list was a matrix for example)
        try
            if isnumeric(cell2mat(value))
            value = cell2mat(value);
            end
        catch
            value = value;
        end
        Ijustgotbackup=1;
    else
        %If we do not go down one level, we want to get some info from this line
        [key,value,flag,lst,lstnr,structlvl] = KeyValueFlag(thisLine,data,ii,lst,lstnr,structlvl,key);
        ii=ii+1;
    end

    %Write the data to the struct
    if Ijustgotbackup %in this case we just got back up one level
        if lst.bool %If I was building a list on this level, I should assign whatever I just got as value, to my list
            lst.data{lstnr}.(key{1}) = value; %Key is the key I had when I went down, so that is where I assign it
            results = lst.data; %And then I can update the result with my updated list
        elseif flag %In this case, I got just back up, and I was not writing a list on this level, so I can direcltly assign my value to the key I had
            keycell = cellstr(key);
            results = setfield(results,keycell{:},value);
        end
    elseif flag %In this case, we did not just got up, but we're gonna write the value to the given key
        keycell = cellstr(key);
        results = setfield(results,keycell{:},value);
    elseif lst.bool %in this case we did not go up, but we're bulding a list, so I'm saving the whole list in each step to the output
        results= lst.data;
    end
end
end






function [key,value,flag,lst,lstnr,structlvl] = KeyValueFlag(thisLine,data,ii,lst,lstnr,structlvl,key)
if startsWith(strtrim(thisLine),"- ") || lst.bool
    %In this case we are in a list
    if ~lst.bool
        lst.structlvl = structlvl(ii); 
        lst.bool = true; %boolean for list
    end       

    if  startsWith(strtrim(thisLine),"- ")
        lstnr = lstnr+1; %only in this case we go to a new item on the list   
        
        % Trim the line to only the data
        sepIndex = find(data{ii}=='-', 1, 'first');
        array = strsplit(data{ii}(sepIndex+2:end), '#');
        array = strtrim(array{1});

        if startsWith(array,'[')
            %Arrays like - [john, Martin], but also - [name, age]: [Rae Smith, 4] 

            if contains(array,']: [')
                %In this case we have a situation like - [name, age]: [Rae Smith, 4]
                %Get the seperation index
                sepIndex = find(array==':', 1, 'first');
                key = GetArr(strtrim(array(1:sepIndex-1)));
                value = GetArr(strtrim(array(sepIndex+1:end)));

                if length(key)~=length(value); error(append('Number of keys and values not equal in line: ',array)); end
                for ll = 1:length(key); lst.data{lstnr,1}.(key{ll})=value{ll}; end
            else
                %In this case we have an ordinary array like - [john, Martin]
                valarr = GetArr(array);
                value = [];
                lst.data{lstnr,1}=valarr;
            end
        elseif startsWith(array,'{')
            %Inline blocks link - {name: John Smith, age: 33}
            sepIndex1 = find(array=='{', 1, 'first');
            sepIndex2 = find(array=='}', 1, 'last');
            valarr = GetArr(array);
            for ll = 1:length(valarr)
                [key,value] = LineToKeyValue(valarr{ll});
                lst.data{lstnr,1}.(key{1})=value;
            end

        elseif contains(array,':')
            %array of objects like - name: john
            [key,value] = LineToKeyValue(array);
            lst.data{lstnr,1}.(key{1})=value;
        else
            %In this case we have an ordinary list - John
            value = [];
            lst.data{lstnr,1} = array;
        end
        flag = false;
        
    elseif contains(thisLine,':') %in this we have a key:value situation
        [key,value] = LineToKeyValue(thisLine);
        flag = false; 
        lst.data{lstnr}.(key{1})=value;
    end    
else
    [key,value] = LineToKeyValue(thisLine);
    
    flag = true;
end
end


function [valarr] = GetArr(array)
%This function loops through all the data of an array and writes it to its
%associated data file (string, numeric, etc.)
aridx = 1;
string = false;
stridx(1) = 2;
for ia = 1:length(array) %Let's loop through the array
    if array(ia) =='"' && ~string
        string = true; %Beginning of string detected
    elseif array(ia) =='"' && string
        string = false; %End of string detected
    end
    if (array(ia) =="," || array(ia) =="]" || ia==length(array)) && ~string %End of array element
        aridx=aridx+1;
        stridx(aridx) = ia; %seperation index
        if aridx==2
            valarr{1,aridx-1}=array(stridx(aridx-1):stridx(aridx)-1);
        else
            valarr{1,aridx-1}=array(stridx(aridx-1)+2:stridx(aridx)-1);
        end
        if ~isnan(str2double(valarr{aridx-1}))
            valarr{1,aridx-1} = str2double(valarr{aridx-1});
        else
            valarr{1,aridx-1} = erase(valarr{1,aridx-1},["'",'"']);
        end
    end
end
% attempt to convert value to numeric type
try
    if isnumeric(cell2mat(valarr))
        valarr = cell2mat(valarr);
    end
catch
end
end

function [key,value] = LineToKeyValue(thisLine)
%This function is called when we know the line has a key, value pair and we
%want to seperate the key from the value and store them in their associated
%arrays
% find the seperator between key and value
sepIndex = find(thisLine==':', 1, 'first');

% get the key name (remove whitespace)
% key{structlvl(ii)} = strtrim(thisLine(1:sepIndex-1));
key{1} = strtrim(thisLine(1:sepIndex-1));

% get the value, ignoring any comments (remove whitespace)
value = strsplit(thisLine(sepIndex+1:end), '#');
value = strtrim(value{1});

if isempty(value)
    % In this case, the value is empty, we will enter one way deeper
    value = value;
elseif  startsWith(value,'[')
    % In this case, the value is an array
    value = GetArr(value);
elseif startsWith(value,"|")
    %Paragraphs of text literal style
    %Check where the text ends
    %         idxend = find(structlvl(ii+1:end)==structlvl(ii),1,'first')+ii-1;
    %         for aa = 1:idxend-ii+1
    %             text
elseif startsWith(value,'>')
    %Paragraphs of text folded styel
elseif startsWith(value,'"') && endsWith(value,'"')
    %In this case, the value is a string, which we store in a cell to be
    %easily compatible with HDF5 format.
    value = value(2:end-1);
elseif startsWith(value,"'") && endsWith(value,"'")
     %In this case, the value is a string, which we store in a cell to be
    %easily compatible with HDF5 format.
    value = value(2:end-1);
end

try
    [convertedValue, success] = str2num(value);
    if success
        value = convertedValue;
    end
    if ~success && ~isempty(value)
        %In this case value is text, and we save it as a string
        value = (value);
    end
catch
    value= value;
end

end