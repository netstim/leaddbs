% mrCache.m
%
%        $Id$
%      usage: mrCache(command,cache,id,data)
%         by: justin gardner
%       date: 05/29/07
%    purpose: implements a simple cache for use with
%             refreshMLRDisplay in mrLoadRet
%       e.g.:
%             c = mrCache('init',10); % init a cache of size 10
%             data = 'whatever';
%             c = mrCache('add',c,'dataid',data); % add data
%             [cachedData c] = mrCache('find',c,'dataid');
%             
%             to clear the cache of all id's with string 'str'
%             c = mrCache('clear',c,'str');
function [retval1 retval2] = mrCache(command,cache,id,data)

% check arguments
if ~any(nargin == [1 2 3 4])
  help mrCache
  return
end

switch lower(command(1))
  % init cache
  case {'i'}
    %see if user passed in a value
    % for how large the cache shohuld be
    if exist('cache','var')
      n = cache;
      clear cache;
      cache.n = n;
    else
      cache.n = 1;
    end
    % init other fields
    cache.data = cell(1,cache.n);
    cache.id = cell(1,cache.n);
    cache.accessed = zeros(1,cache.n); 
    % and return
    retval1 = cache;
   % add an item to the cache
   case {'a'}
    % see if it is already in cache
    cacheNum = find(strcmp(cache.id,id));
    if isempty(cacheNum)
      % if not, add it in to the place that has
      % been accessed the longest time ago
      cacheNum = first(find(max(cache.accessed)==cache.accessed));
    end
    % set the cache data and id
    if ~isempty(cacheNum)
      cache.id{cacheNum} = id;
      cache.data{cacheNum} = data;
      % set the accessed field so that everyone
      % else gets 1 added 
      cache.accessed = cache.accessed+1;
      cache.accessed(cacheNum) = min(cache.accessed)-1;
    end
    % return the cache
    retval1 = cache;
   % find an item in the cache
   case {'f'}
    cacheNum = find(strcmp(cache.id,id));
    if isempty(cacheNum)
      retval1 = [];
      retval2 = cache;
    else
      % set the accessed field so that this one is youngest
      cache.accessed = cache.accessed+1;
      cache.accessed(cacheNum) = min(cache.accessed)-1;
      % return the data and cache
      retval1 = cache.data{cacheNum};
      retval2 = cache;
    end
  % clear cache
  case {'c'}
   % look for any id that contains the string
   for i = 1:length(cache.id)
      % and clear it
      if ~isempty(findstr(cache.id{i},id))
	cache.id{i} = '';
	cache.data{i} = [];
	cache.accessed(i) = max(cache.accessed)+1;
      end
    end
    retval1 = cache;
end

