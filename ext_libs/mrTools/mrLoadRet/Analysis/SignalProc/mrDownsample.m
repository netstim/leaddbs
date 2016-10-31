function ds = mrDownsample(s, factor, sampleNumber)
  % downsample a signal by factor. the original and downsampled signals
  % have equal integrals. sampleNumber, if it exists, is the old sample at which the new value is estimated
  % (1<=sampleNumber<=factor; default: sampleNumber=factor)
  %        $Id$
  % to downsample conserving the amplitude use: 
  %     ds = downsample(s, factor) / factor;
  
  if size(s,1)==1 %if only one row, we've been passed a row vector
     transpose_vector = 1;
     s = s';
  else
     transpose_vector = 0;
  end
  if ieNotDefined('sampleNumber') || sampleNumber < 1 || sampleNumber > factor
    sampleNumber = floor(factor/2)+1;
  end
  %because we're sampling from the integral of the signal, we need to pad with zeros
  % so that sampleNumber samples the correct sample (the middle of the values
  %used in the difference)
  a = [zeros(floor(factor/2)+1,size(s,2));cumsum(s, 1)];
  ds = diff(a(sampleNumber:factor:end,:));

  if transpose_vector
     ds = ds';
  end