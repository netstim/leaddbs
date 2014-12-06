function z = ea_normal(data, varargin)
% Sacha van Albada, Dec 11 2006
%----------------------------------------------------------------------
% This function takes a set of observations and transforms these to the
% normal distribution.
%
% If the data contain ties, the default is to map these to equal values,
% so that the output is not completely normal.
% However, there is an option to split ties.
%----------------------------------------------------------------------
% If you use this code, please reference:
%
% Van Albada, S.J., Robinson P.A. (2006) Transformation of arbitrary
% distributions to the normal distribution with application to EEG
% test-retest reliability. J Neurosci Meth (in press)
%----------------------------------------------------------------------
% Syntax:
% normal(data,n_datasets,n_headerlines,col_sep,means,stdevs,split_ties)
%----------------------------------------------------------------------
% Inputs:
%
% data:             An array or datafile containing nonnormal data. If
%                   there is more than one data set, the data should appear
%                   in columns.
%
% n_datasets:       The number of data sets, from 1 to 20. Defaults to 1.
%
% n_headerlines:    If the data are in a file, the number of header lines
%                   can be specified. Defaults to 0, and ignored if the
%                   data are in array form.
%
% col_sep:          A column delineator. Defaults to a space, and ignored
%                   if the data are in array form.
%
% means:            The desired mean(s) of the output data.
%                   Either a single number or an array of values.
%                   Defaults to 0 for each data set.
%
% stdevs:           The desired standard deviation(s) of the output data.
%                   Either a single number or an array of values.
%                   Defaults to 1 for each data set.
%
% split_ties:       If set, this maps ties arbitrarily to different values,
%                   so that the output data are as close as possible to the
%                   normal distribution.
%----------------------------------------------------------------------
% Output:
%
% The array of normalized values
%----------------------------------------------------------------------
% Examples:
%
% a = [6.1;4.2;5.4;2.9;4.5;5.4;2.3;4.7;6.5;7.1;2.0;7.5;3.2;1.9;5.4;3.6];
%
% normal(a)
%
% ans =
%
%    0.7764
%   -0.2372
%    0.2372
%   -0.7764
%   -0.0784
%    0.2372
%   -1.0100
%    0.0784
%    1.0100
%    1.3180
%   -1.3180
%    1.8627
%   -0.5791
%   -1.8627
%    0.2372
%   -0.4023
%
% This creates a distribution that is very close to the standard normal
% distribution. In order to split ties, and obtain mean 4.3 and std 2,
% type the following:
%
% normal(a, 1, 4.3, 2, 'TRUE')
%
% Here, the second argument indicates that the data consist of a single
% column.
%
% ans =
%
%    5.8528
%    3.8256
%    4.7744
%    2.7472
%    4.1432
%    5.1045
%    2.2800
%    4.4568
%    6.3200
%    6.9360
%    1.6640
%    8.0255
%    3.1417
%    0.5745
%    5.4583
%    3.4955
%
% Similarly, with multiple columns:
%
% b=[[5;2;3;4;5;2;4;3;6],[3;4;1;5;2;6;2;3;7]];
%
% normal(b, 2, [5 2], [3 4], 'FALSE')
%
% produces
%
% ans =
%
%    7.2021    1.1858
%    0.9313    3.3862
%    3.8151   -3.9251
%    5.5086    4.5839
%    7.2021   -1.4854
%    0.9313    6.0574
%    5.5086   -1.4854
%    3.8151    1.1858
%   10.0859    8.4970
%----------------------------------------------------------------------
% varargin can contain: n_datasets, n_headerlines, col_sep, means, stdevs,
% split_ties

error(nargchk(1,7, nargin))   % Prints an error message and returns if the
% nr of input arguments supplied is less than
% 1 or more than 7
%----------------------------------------------------------------------
switch nargin
    case 1 % Only 'data' given
        % one data set, no header lines, mean 0, stdev 1, ties not split
        
        if(ischar(data))
            fid  = fopen(data);
            if(fid == -1)
                sprintf('No valid data file selected')
                return
            end
            distrib = fscanf(fid, '%f', inf);  % Read data to end of file
            fclose(fid);
            if(size(distrib,2)==0)
                sprintf('File does not start with numeric data')
                return
            end
            
            if(size(distrib,2)>1)
                distrib = distrib(:,1);
                str1 = 'Warning: array dimensions do not match specified';
                str2 = ' number of data sets.\n';
                str3 = 'The first column of data read in.';
                
                sprintf(strcat(str1,str2,str3))
            end
        elseif(isnumeric(data))
            if(size(data,2) == 1)    % Check for correct nr of columns
                distrib = data;
            elseif(size(data,2) > 1)
                distrib = data(:,1);
                str1 = 'Warning: array dimensions do not match specified';
                str2 = ' number of data sets.\n';
                str3 = 'The first column of data read in.';
                
                sprintf(strcat(str1,str2,str3))
            else
                sprintf('Array does not contain data.')
                return
            end
            
        else
            sprintf('Input data should be given as a filename or array.')
            return
        end
        
        z = no_split(distrib);
        
        %----------------------------------------------------------------------
    case 2              % Two arguments: data and n_datasets
        % no header lines, mean(s) 0, stdev(s) 1, ties not split
        n_datasets = varargin{1};
        n_headerlines = 0;
        col_sep = ' ';
        
        distrib = read(data, n_datasets,n_headerlines, col_sep);
        z = no_split(distrib);
        
        %----------------------------------------------------------------------
    case 3   % Three arguments: either data, n_datasets, n_headerlines
        % or data, n_datasets, split_ties.
        % mean(s) 0, stdev(s) 1
        n_datasets = varargin{1};
        col_sep = ' ';
        
        if(isnumeric(varargin{2}))     % The third argument is the number of
            % header lines; ties are not split
            n_headerlines = varargin{2};
            
            distrib = read(data, n_datasets, n_headerlines, col_sep);
            z = no_split(distrib);
            
        elseif(strcmp(varargin{2},'TRUE') | strcmp(varargin{2},'FALSE'))
            split_ties = varargin{2};
            n_headerlines = 0;
            distrib = read(data, n_datasets, n_headerlines, col_sep);
            
            if(strcmp(split_ties,'TRUE'))
                z = split(distrib);
            else                                 % If split_ties == 'FALSE'
                z = no_split(distrib);
            end
        else
            str1 = 'Please supply a valid third argument:\n';
            str2 = 'either the number of header lines of your data file\n';
            str3 = 'or \''TRUE\'' or \''FALSE\'' indicating';
            str4 = ' whether ties should be split.';
            sprintf(strcat(str1,str2,str3,str4))
            return
        end
        
        %----------------------------------------------------------------------
    case 4   % Four arguments: either data, n_datasets, n_headerlines,
        % col_sep, or data, n_datasets, means, stdevs.
        % Ties are not split
        n_datasets = varargin{1};
        if(isnumeric(varargin{3}))    % The fourth argument gives stdev
            if(length(varargin{2})==length(varargin{3}) &...
                    length(varargin{2})==n_datasets)
                
                means = varargin{2};
                stdevs = varargin{3};
                n_headerlines = 0;
                col_sep = ' ';
                
                distrib = read(data, n_datasets,n_headerlines,col_sep);
                ns = no_split(distrib);
                
                z = zeros(size(ns,1), size(ns,2));
                for(i = 1:size(ns,2))
                    z(:,i) = means(i)-(stdevs(i)/std(ns(:,i)))*mean(ns(:,i))+...
                        (stdevs(i)/std(ns(:,i)))*ns(:,i);
                    % Since ns will not have mean=0 and sd=1 in case of ties
                end
                
            else sprintf('Dimensions specified by arguments do not agree.')
                return
            end
        elseif(ischar(varargin{3}))   % The fourth argument gives col_sep
            col_sep = varargin{3};
            n_headerlines = varargin{2};
            
            distrib = read(data, n_datasets, n_headerlines, col_sep);
            z = no_split(distrib);
        else   sprintf('Fourth argument not valid.')
            return
        end
        %----------------------------------------------------------------------
    case 5    % Five arguments: either data, n_datasets,n_headerlines,
        % col_sep, split_ties
        % or data, n_datasets, means, stdevs, split_ties
        n_datasets = varargin{1};
        if (strcmp(varargin{4},'TRUE')|strcmp(varargin{4},'FALSE'))
            split_ties = varargin{4};
        else
            sprintf('Fifth argument not valid')
            return
        end
        
        if(isnumeric(varargin{3}))    % The fourth argument gives stdev
            if(length(varargin{2})==length(varargin{3}) & ...
                    length(varargin{2})==n_datasets)
                
                means = varargin{2};
                stdevs = varargin{3};
                n_headerlines = 0;
                col_sep = ' ';
                
                distrib = read(data, n_datasets, n_headerlines, col_sep);
                
                if(strcmp(varargin{4},'TRUE'))
                    s = split(distrib);
                    z = zeros(size(s,1), size(s,2));
                    for(i = 1:size(s,2))
                        z(:,i)=means(i)+stdevs(i)*s(:,i);
                    end
                    
                else                         % If split_ties == 'FALSE'
                    ns = no_split(distrib);
                    z = zeros(size(ns,1), size(ns,2));
                    for(i = 1:size(ns,2))
                        z(:,i) = means(i)-(stdevs(i)/std(ns(:,i)))*mean(ns(:,i))+...
                            (stdevs(i)/std(ns(:,i)))*ns(:,i);
                    end
                end
                
            else
                sprintf('Dimensions specified by arguments do not agree.')
                return
            end
        elseif(ischar(varargin{3}))   % The fourth argument gives col_sep
            
            col_sep = varargin{3};
            n_headerlines = varargin{2};
            
            distrib = read(data, n_datasets, n_headerlines, col_sep);
            if(strcmp(varargin{4},'TRUE'))
                z = split(distrib);
            else
                z = no_split(distrib);
            end
            
        else   sprintf('Fourth argument not valid.')
            return
        end
        %----------------------------------------------------------------------
    case 6        % Six arguments: data, n_datasets, n_headerlines,
        % col_sep, means, stdevs
        n_datasets = varargin{1};
        n_headerlines = varargin{2};
        col_sep = varargin{3};
        means = varargin{4};
        stdevs = varargin{5};
        
        if(length(varargin{2})==length(varargin{3}) & ...
                length(varargin{2})==n_datasets)
            distrib = read(data, n_datasets, n_headerlines, col_sep);
            ns = no_split(distrib);
            z = zeros(size(ns,1), size(ns,2));
            for(i = 1:size(ns,2))
                z(:,i) = means(i)-(stdevs(i)/std(ns(:,i)))*mean(ns(:,i))+...
                    (stdevs(i)/std(ns(:,i)))*ns(:,i);
            end
        else sprintf('Dimensions specified by arguments do not agree.')
            return
        end
        %----------------------------------------------------------------------
    case 7        % Seven arguments: data, n_datasets, n_headerlines,
        % col_sep, means, stdevs, and split_ties
        n_datasets = varargin{1};
        n_headerlines = varargin{2};
        col_sep = varargin{3};
        means = varargin{4};
        stdevs = varargin{5};
        
        if(length(varargin{4})==length(varargin{5}) & ...
                length(varargin{5})==n_datasets)
            distrib = read(data, n_datasets, n_headerlines, col_sep);
            
            if (strcmp(varargin{6},'TRUE') | strcmp(varargin{6},'FALSE'))
                
                if (strcmp(varargin{6},'TRUE'))
                    s = split(distrib);
                    z = zeros(size(s,1), size(s,2));
                    for(i = 1:size(s,2))
                        z(:,i)=means(i)+stdevs(i)*s(:,i);
                    end
                    
                else                                 % If split_ties == 'FALSE'
                    ns = no_split(distrib);
                    z = zeros(size(ns,1), size(ns,2));
                    for(i = 1:size(ns,2))
                        z(:,i) = means(i)-(stdevs(i)/std(ns(:,i)))*mean(ns(:,i))+...
                            (stdevs(i)/std(ns(:,i)))*ns(:,i);
                    end
                end
            else
                sprintf('Last argument not valid')
                return
            end
            
        else
            sprintf('Dimensions specified by arguments do not agree.')
            return
        end
end

%----------------------------------------------------------------------
% Functions

function distrib = read(data, n_datasets, n_headerlines, col_sep)
if(ischar(data))              % If the data are in a file
    fid  = fopen(data);
    
    if(n_datasets>=1 & n_datasets <= 20)
        i=1;
        while(i <= n_headerlines)  % Bring position to end of header
            fgetl(fid);
            i=i+1;
        end
        
        j=1;
        scan_str = '';
        while (j <= n_datasets-1)
            scan_str = [scan_str,'%f',col_sep];
            j=j+1;
        end
        scan_str = [scan_str,'%f\n'];
        
        distrib = fscanf(fid, scan_str,[n_datasets inf]);
        % The columns are converted into rows
        fclose(fid);
        
    else
        sprintf('Please specify a number of data sets between 1 and 20.')
        return
    end
    distrib = distrib';       % Convert back into columns
elseif(isnumeric(data))            % If the data are in array form
    
    if(n_datasets == 1)
        if(size(data,2) == 1)   % Check for correct nr of columns
            if(n_headerlines ~= 0)
                sprintf('Argument n_headerlines ignored.')
            end
            if(~strcmp(col_sep, ' ') & ~strcmp(col_sep, ''))
                sprintf('Argument col_sep ignored.')
            end
            distrib = data;
        elseif(size(data,2) > 1)
            distrib = data(:,1);
            str1 = 'Warning: array dimensions do not match specified';
            str2 = ' number of data sets.\n'
            str3 = 'The first column of data read in.'
            
            sprintf(strcat(str1,str2,str3))
        else
            sprintf('Array does not contain data.')
            return
        end
        
    elseif(n_datasets <= 20)
        if(size(data,2) == n_datasets) % Check for correct nr of columns
            distrib = data;
        elseif(size(data,2) > n_datasets)
            distrib = data(:,1:n_datasets);
            str1 = 'Warning: array dimensions do not match specified';
            str2 = ' number of data sets.\n';
            str3 = 'The first %d columns of data read in.';
            
            sprintf(strcat(str1,str2,str3), n_datasets)
        else
            sprintf('Array does not contain sufficient data.')
            return
        end
    else
        sprintf('Please specify a number of data sets between 1 and 20.')
        return
    end
    
else      sprintf('Input data should be given as a filename or array.')
    return
end

%-----------------------------------------------------------------------
function ns = no_split(distrib)
ns = zeros(size(distrib,1), size(distrib,2));

for(i = 1:size(distrib,2))
    
    [b,m,n] = unique(distrib(:,i)); % returns unique values of distrib,
    % where b = distrib(m) and distrib = b(n)
    % n are rank numbers for the unique observations
    [sorted, index]=sort(n);
    new_sorted = sorted;
    
    ties=0;
    % Make the maximum rank the total nr of observations,
    % not the nr of unique observations
    for(j = 1:(size(sorted)-1))
        if(sorted(j)==sorted(j+1))
            ties=ties+1;
        else
            new_sorted(j+1:size(sorted))=sorted(j+1:size(sorted))+ties;
        end
    end
    
    [x, index2] = sort(index);
    rank = new_sorted(index2);
    % The cumulative distribution function
    
    CDF = rank./length(distrib) - 1/(2*length(distrib));
    
    % Apply the inverse Rosenblatt transformation
    % z will not be completely normal if there are ties
    ns(:,i) = sqrt(2)*erfinv(2*CDF - 1);
    
end

%-----------------------------------------------------------------------
function s = split(distrib)
s = zeros(size(distrib,1), size(distrib,2));

for(i = 1:size(distrib,2))
    % Sort the data in ascending order and retain permutation indices
    [sorted,index]  = sort(distrib(:,i));
    % Make 'rank' the rank number of each observation
    [x, rank]       = sort(index);
    % The cumulative distribution function
    CDF = rank./length(distrib) - 1/(2*length(distrib));
    % Apply the inverse Rosenblatt transformation
    % z is now normally distributed
    s(:,i) = sqrt(2)*erfinv(2*CDF - 1);
    
end