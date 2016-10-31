% disppercent.m
%
%         by: justin gardner
%       date: 10/05/04
%      usage: disppercent(percentdone,message)
%    purpose: display percent done
%             Start by calling with a negative value:
%             disppercent(-inf,'Message to display');
%
%             Update by calling with percent done:
%             disppercent(0.5);
% 
%             Finish by calling with inf (elapsedTime is in seconds):
%             elapsedTime = disppercent(inf);
%
%             If you want to change the message before calling with inf:
%             disppercent(0.5,'New message to display');
% 
%             Also, if you have an inner loop within an outer loop, you
%             can call like the following:
%             n1 = 15;n2 = 10;
%             disppercent(-1/n1); % init with how much the outer loop increments
%             for i = 1:n1
%               for j = 1:n2
%                  pause(0.1);
%                  disppercent((i-1)/n1,j/n2);
%               end
%               disppercent(i/n1,sprintf('Made it through %i/%i iterations of outer loop',i,n1));
%             end
%             disppercent(inf);
%
%       e.g.:
%
%disppercent(-inf,'Doing stuff');for i =  1:30;pause(0.1);disppercent(i/30);end;elapsedTime = disppercent(inf);
function retval = disppercent(percentdone,mesg)

retval = nan;
% check command line arguments
if ((nargin ~= 1) && (nargin ~= 2))
  help disppercent;
  return
end

% global for disppercent
global gDisppercent;

% if this is an init then remember time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (percentdone < 0)
  % global to turn off printing
  global gVerbose;
  gDisppercent.verbose = gVerbose;
  if ~gVerbose,return,end
  % set starting time
  gDisppercent.t0 = clock;
  % default to no message
  if (nargin < 2)
    mrDisp(sprintf('00.00%% (00 hrs 00 min 00 sec)'));
    gDisppercent.mesg = '';
  else
    mrDisp(sprintf('%s 00.00%% (00 hrs 00 min 00 sec)',mesg));
    gDisppercent.mesg = mesg;
  end    
  if isinf(percentdone)
    gDisppercent.increment = 0;
  else
    gDisppercent.increment = abs(percentdone);
  end
    
% display total time at end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (percentdone == inf)
  if ~gDisppercent.verbose,return,end
  % reprint message if necessary
  if (nargin == 2) && ischar(mesg)
    reprintMessage(mesg);
  end
  % get elapsed time
  elapsedTime = etime(clock,gDisppercent.t0);
  % separate seconds and milliseconds
  numms = round(1000*(elapsedTime-floor(elapsedTime)));
  numsecs = floor(elapsedTime);
  % if over a minute then display minutes separately
  if numsecs>60
    nummin = floor(numsecs/60);
    numsecs = numsecs-nummin*60;
    % check if over an hour
    if nummin > 60
      numhours = floor(nummin/60);
      nummin = nummin-numhours*60;
      timestr = sprintf('%i hours %i min %i secs %i ms',numhours,nummin,numsecs,numms);
    else
      timestr = sprintf('%i min %i secs %i ms',nummin,numsecs,numms);
    end
  else
    timestr = sprintf('%i secs %i ms',numsecs,numms);
  end
  % display time string
  mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\btook %s\n',timestr));
  retval = elapsedTime;
% otherwise show update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  if ~gDisppercent.verbose,return,end
  
  % keep track of what increments percent done is being called in
  if (gDisppercent.increment == 0) & (percentdone > 0)
    gDisppercent.increment = percentdone;
  end
  % a negative value on percent done, simply means
  % how much each increment of percent done is.
  if percentdone < 0
    gDisppercent.increment = -percentdone;
    percentdone = 0;
  end
  % see if we should reprint message
  newmesg = '';
  if nargin == 2
    if ischar(mesg)
      newmesg = reprintMessage(mesg);
    % otherwise if the second argument is a number,
    % it means we have a secondary value for percent done
    % i.e. the first number is the large increments and
    % the second number is what percentage of the large
    % increments has been done (useful for when you are
    % doing loops within loops).
    elseif isscalar(mesg)
      % make percent done into value computed by summing
      % percentdone with the increment passed in mesg.
      percentdone = percentdone + gDisppercent.increment*mesg;
    end
  end
      
  % avoid things that will end up dividing by 0
  if (percentdone >= 1)
    percentdone = .9999;
  elseif (percentdone <= 0)
    percentdone = 0.0001;
  end
  elapsedTime = etime(clock,gDisppercent.t0);

  % display percent done and estimated time to end
  if ~isempty(newmesg)
    % always display if there is a new message
    mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s%05.2f%% (%s)',newmesg,floor(10000*percentdone)/100,disptime(elapsedTime*(1/percentdone - 1))));
  % display only if we have update by a least a percent or if at least 1 second has elapsed since last display
  elseif (gDisppercent.percentdone ~= floor(100*percentdone)) || floor(elapsedTime)~=gDisppercent.elapsedTime
    mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05.2f%% (%s)',floor(10000*percentdone)/100,disptime(elapsedTime*(1/percentdone - 1))));
  end
end
% remember current percent done
gDisppercent.percentdone = floor(100*percentdone);
gDisppercent.elapsedTime = floor(etime(clock,gDisppercent.t0));


%%%%%%%%%%%%%%%%%%
%%   disptime   %%
%%%%%%%%%%%%%%%%%%
function retval = disptime(t)

hours = floor(t/(60*60));
minutes = floor((t-hours*60*60)/60);
seconds = floor(t-hours*60*60-minutes*60);

retval = sprintf('%02i hrs %02i min %02i sec',hours,minutes,seconds);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   reprintMessage   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function newmesg = reprintMessage(mesg)

global gDisppercent;

newmesg = '';

if ~strcmp(mesg,gDisppercent.mesg)
  % first clear old message
  mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s%s                             ',repmat(sprintf('\b'),1,length(gDisppercent.mesg)+1),repmat(sprintf(' '),1,length(gDisppercent.mesg)+1)));
  % print <or return> new message
  if nargout == 1
    newmesg = sprintf('%s%s ',repmat(sprintf('\b'),1,length(gDisppercent.mesg)+1),mesg);
  else
    mrDisp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s%s                              ',repmat(sprintf('\b'),1,length(gDisppercent.mesg)+1),mesg));
  end
  gDisppercent.mesg = mesg;
end
