% myhighpass.m
%
%      usage: eventRelatedHighpass(d,cutoff)
%         by: justin gardner
%       date: 06/27/05
%
function d = eventRelatedHighpass(d,cutoff,runtype)

if ~any(nargin == [2 3])
  help eventRelatedHighpass;
  return
end
if exist('runtype')~=1,runtype = 'both';,end

% check runtype
if ~any(strcmp(runtype,{'both','init','filter'}))
  disp(sprintf('UHOH: Runtype must be: both, init or filter'));
  return
end
dispfig = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the filter if we are run with
% runtype set to both or to init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(runtype,'both') || strcmp(runtype,'init')
  % first get length of filter
  n = d.dim(4);

  % get nyquist frequencies
  freqmax = 1/(2*d.tr);
  freqmin = 1/(n*d.tr);
  freqdelta = 1/(n*d.tr);

  % get the frequencies
  freqs = 0:freqdelta:(freqdelta*(n-1));

  % get the times
  times = 0:d.tr:(d.tr*(n-1));

  % now make a high-pass cutoff
  hipassfilter = ones(1,n);

  % knock out frequncies below cutoff
  hipassfilter(freqs<cutoff) = 0;

  % and smoothly go to 1 as a gaussian for edges
  smoothedge = 1-eventRelatedGauss([1 cutoff cutoff/2 0],freqs);

  % and add that smooth edge to the square filter
  hipassfilter(find(freqs>cutoff)) = smoothedge(freqs>cutoff);

  % now add the mirror reverse to take care of
  % the negative frequencies, don't need to copy the
  % first element since that is for the dc component
  if (isodd(n))
    hipassfilter(n:-1:(round(n/2)+1)) = hipassfilter(2:round(n/2));
  else
    hipassfilter(n:-1:(n/2+2)) = hipassfilter(2:n/2);
  end  

  % notch out highest frequency that we get with sense processing
  if d.notchFilterForTSense
    if iseven(n)
      if d.notchFilterForTSense == 2
	hipassfilter((n/2)+1) = 0;
	disp(sprintf('(eventRelatedHighpass) Applying notch filter for acceleration of x%i',d.notchFilterForTSense));
      elseif d.notchFilterForTSense == 4
	hipassfilter((n/2)+1) = 0;
	hipassfilter((n/4)+1) = 0;
	hipassfilter((3*n/4)+1) = 0;
	disp(sprintf('(eventRelatedHighpass) Applying notch filter for acceleration of x%i',d.notchFilterForTSense));
      else
	mrWarnDlg(sprintf('(eventRelatedHighpass) Notch for tSense not implemented yet for acceleration of %i',d.notchFilterForTSense));
      end
    else
      mrWarnDlg(sprintf('(eventRelatedHighpass) Notch for tSense expects that the time series is even (because it was accelerated) - still running, but filtering out the highest frequency component for an odd length time series may not completely remove the artifact'));
      hipassfilter((n+1)/2:((n+1)/2)+1)= 0;
    end
  end
  % give back the DC
  hipassfilter(1) = 1;
  %hipassfilter(1) = 0;

  %remember filter
  d.hipassfilter = hipassfilter;
  d.hipasscutoff = cutoff;

  % display the filter if called for
  if (dispfig)
    mlrSmartfig('myhighpass');
    subplot(1,2,1);
    plot(freqs,hipassfilter);
    hold on
    vline(cutoff,'r-');
    vline(freqs(first(find(hipassfilter > 0.999))),'g-');
    xlabel('frequency (Hz)');
    ylabel('magnitude');
  
    % display the convolution kernel of the filter
    subplot(1,2,2);
    plot(times,abs(fft(hipassfilter)));
    xlabel('time (sec)');
    ylabel('filter magnitude');
    xaxis(0,60);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter if
% runtype set to both or filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(runtype,'both') || strcmp(runtype,'filter')
  if ~isfield(d,'hipassfilter')
    disp(sprintf('(eventRelatedHighpass): Must run myhighpass(d,cutoff,''init'') first'));
    return
  end
  hipassfilter = d.hipassfilter;
  if isfield(d,'hipasscutoff')
    % go through the data, detrend and apply filter in fourier domain
    disppercent(-inf,sprintf('(eventRelatedHighpass) Applying temporal hipass filter (cutoff=%0.03f Hz)',d.hipasscutoff));
  else
    disppercent(-inf,sprintf('(eventRelatedHighpass) Applying temporal hipass filter'))
  end    
  if (0)
    for i = 1:d.dim(1)
      disppercent(i/d.dim(1));
      for j = 1:d.dim(2)
	for k = 1:d.dim(3)
	  timecourse = squeeze(d.data(i,j,k,:));
	  timecourse = eventRelatedDetrend(timecourse);
	  timecourse = ifft(fft(timecourse) .* hipassfilter');
	  d.data(i,j,k,:) = real(timecourse);%abs(timecourse).*cos(angle(timecourse));
	end
      end
    end
  elseif (~isempty(d.data))
    % faster way, do across first dimension
    % replicate the high pass filter
    for i = 1:d.dim(1)
      newhipassfilter(i,:) = hipassfilter;
    end
    % detrend and high pass filter
    for k = 1:d.dim(3)
      disppercent(k/d.dim(3));
      for j = 1:d.dim(2)
	timecourses = squeeze(d.data(:,j,k,:));
	timecourses = eventRelatedDetrend(timecourses')';
	timecourses = ifft(fft(timecourses,[],2) .* newhipassfilter,[],2);
	d.data(:,j,k,:) = real(timecourses);
      end
    end
  elseif (isfield(d,'roidata'))
    for i = 1:size(d.roidata,1)
      disppercent(i/size(d.roidata,1))
      timecourse = squeeze(d.roidata(i,:));
      timecourse = eventRelatedDetrend(timecourse);
      timecourse = ifft(fft(timecourse) .* hipassfilter');
      d.roidata(i,:) = real(timecourse);
    end
  end
  disppercent(inf);
end


 
