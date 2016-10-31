function hdrSingleValue(thisView,overlayNum,scan,x,y,s,roi)
% hdrSingleValue.m
%
%       $Id$
%      usage: hdrSingleValue(thisView,overlayNum,scan,x,y,s,roi)
%         by:  julien besle
%       date: 20/05/2010
%    purpose: extract overlays from GLM/event-related analysis
%

% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data
% get the analysis structure
analysis = viewGet(thisView,'analysis');
if ~isfield(analysis,'d') || isempty(analysis.d)
  disp(sprintf('(eventRelatedPlot) No Event-Related analysis'));
  return
end


set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;

nScans=viewGet(thisView,'nscans');
preSelect = zeros(1,nScans);
for iScan = 1:nScans
  if ~isempty(analysis.d{iScan})
    preSelect(iScan)=1;
  end
end
if sum(preSelect)>1
  scanlist = selectInList(thisView,'scan','Select scans to include');
else
  scanlist = find(preSelect);
end
if isempty(scanlist)
    return;
end

d = analysis.d{scanlist(1)};
params.sampleRange = ['1:' num2str(size(d.ehdr,5))];
valuteTypeMenu = {'er2beta','average','dotProduct','max','delay','diffGammaDelay'};

params = {...
            {'valueType',valuteTypeMenu,'type=popupmenu'},...
            {'sampleRange',params.sampleRange,'Type of the combination function to apply.  s, or a custom function'},...
         };

params = mrParamsDialog(params, 'HDR to single value parameters');
% Abort if params empty
if ieNotDefined('params'),return,end

params.sampleRange = eval(params.sampleRange);

stimNames = d.stimNames;
maxNumberConditions = 0;

singleValue = cell(1,nScans);
for iScan = scanlist
   d = analysis.d{iScan};
   if ~isempty(d)
      dimensions = size(d.ehdr);
      maxNumberConditions = max (maxNumberConditions,dimensions(4)); %this is in case there are different number of conditions between scans

      switch(params.valueType)
         case 'average'  %(averaged on sampleRange)

            singleValue{iScan} = mean(d.ehdr(:,:,:,:,params.sampleRange),5);
            
         case 'max'  %max value (across sampleRange)

            singleValue{iScan} = max(d.ehdr(:,:,:,:,params.sampleRange),[],5);

         case {'dotProduct', 'er2beta'}

            if strcmp(params.valueType,'er2beta')
               number_of_samples = size(d.scm,1);
               design_matrix = zeros(number_of_samples,dimensions(4));
               for i_event = 1:dimensions(4)
                  design_matrix(d.stimvol{i_event},i_event) = 1;
               end
            end
            singleValue{iScan} = NaN(dimensions(1:4));
            h_wait = waitbar(0,['Computing ' params.valueType ' for scan ' num2str(iScan)]);
            for z = 1:dimensions(3)
               waitbar(z/dimensions(3),h_wait);

               if strcmp(params.valueType,'er2beta')
                   tseries = squeeze(loadTSeries(thisView,iScan,z));
                   tseries = tseries - repmat(mean(tseries,3),[1 1 size(tseries,3)]);
                   deconvolution_matrix = nan(number_of_samples,dimensions(4));
               end

               for x = 1:dimensions(1)
                  for y = 1:dimensions(2)
                     if any(any(isnan(d.ehdr(x,y,z,:,:))))
                     else
                        data = squeeze(d.ehdr(x,y,z,:,params.sampleRange));
                        switch (params.valueType)
                           case 'dotProduct'
                              mean_data = mean(data); %compute mean response across conditions for this voxel
                                                     %compute the inner_product of each condition with the mean response,
                              singleValue{iScan}(x,y,z,:) = data*mean_data'/(norm(mean_data)^2);    %normalize by the square norm and chose the highest value among conditions 
                           
                           case 'er2beta'
                              for i_event = 1:dimensions(4)
                                 this_column = conv(design_matrix(:,i_event)',mean(data))';
                                 deconvolution_matrix(:,i_event) = this_column(1:number_of_samples);
                              end
                              singleValue{iScan}(x,y,z,:) = pinv(deconvolution_matrix)*squeeze(tseries(x,y,:));
                        end
                              
                     end
                  end
               end
            end

            close(h_wait);
            
         otherwise
            mrWarnDlg(['(hdrSingleValue) Value type ''' params.valueType ''' is not yet implemented']);
            return

      end
   else
       mrWarnDlg(sprintf('(getOverlayFromGlmAnalysis) skipping scan %i because analysis is empty',iScan));
   end

end

maxValue = -inf;
minValue = inf;
for iScan = 1:size(singleValue,1)
   maxValue = max(maxValue,max(singleValue{iScan}(:)));
   minValue = min(minValue,min(singleValue{iScan}(:)));
end

outputOverlay.groupName = viewGet(thisView,'groupName');
outputOverlay.function = 'hdrSingleValue';
outputOverlay.type = 'combination';
outputOverlay.params = params;
for iOverlay = 1:dimensions(4)
   outputOverlay.name = [params.valueType '(' stimNames{iOverlay} ')'];
   outputOverlay.data = cell(1,nScans);

   for iScan = 1:nScans
     if ~isempty(singleValue{iScan})
      outputOverlay.data{iScan} = singleValue{iScan}(:,:,:,iOverlay);
     end
   end

   outputOverlay.type = 'beta';
   outputOverlay.date = datestr(now);
   outputOverlay.range = [minValue maxValue];
   outputOverlay.clip = [minValue maxValue];
   thisView = viewSet(thisView,'newoverlay',outputOverlay);  
   
end

set(viewGet(thisView,'figNum'),'Pointer','arrow');drawnow;
