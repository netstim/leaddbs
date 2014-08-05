function ea_show_normalization(options)

for export=1:3-(options.modality-1) % if CT, only do 1:2, if MR, do 1:3.
    try
    switch export
        case 1
            checkf=[options.root,options.prefs.patientdir,filesep,options.prefs.prenii,',1'];
            
            outf=['check_',options.prefs.prenii];
            
            suff='_pre_tra';
        case 2
            checkf=        [options.root,options.prefs.patientdir,filesep,options.prefs.tranii,',1'];
            
            outf=['check_',options.prefs.tranii];
            
            suff='_tra';
        case 3
            checkf=        [options.root,options.prefs.patientdir,filesep,options.prefs.cornii,',1'];
            outf=['check_',options.prefs.cornii];
            suff='_cor';
    end
    matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'templates',filesep,'mni_wires.nii,1']
        checkf};
    matlabbatch{1}.spm.util.imcalc.output = outf;
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = ['5*i1+i2'];
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;
   
    
    nii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,outf]);
 
    h1=figure('name',['Normalization results for ',options.prefs.patientdir,'_',outf],'NumberTitle','off');
    set(gcf,'color','w')
    imagesc(flipud(squeeze(nii.img(:,:,round(end/2)))'));
    axis('equal')
    axis('off')
    colormap gray
    
    tightfig;
    
    h2=figure('name',['Normalization results for ',options.prefs.patientdir,'_',outf],'NumberTitle','off');

    
    subplot(2,1,1);
    imat=scale_image(flipud(squeeze(nii.img(:,round(end/2),:))'),[2,1]);
    imagesc(imat);
    axis('equal')
    axis('off')
    subplot(2,1,2);
    imat=scale_image(flipud(squeeze(nii.img(round(end/2),:,:))'),[2,1]);
    imagesc(imat);
    axis('equal')
    axis('off')
    colormap gray

    tightfig;
    saveas(h1,[options.root,options.prefs.patientdir,filesep,'eauto_normalization_check',suff,'_axial.png']);
    saveas(h2,[options.root,options.prefs.patientdir,filesep,'eauto_normalization_check',suff,'_corsag.png']);
    end
end





function hfig = tightfig(hfig)

% Copyright (c) 2011, Richard Crozier
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
%
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
%
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

if nargin == 0
    hfig = gcf;
end

% There can be an issue with tightfig when the user has been modifying
% the contnts manually, the code below is an attempt to resolve this,
% but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
set(hfig, 'WindowStyle', 'normal');

% 1 point is 0.3528 mm for future use

% get all the axes handles note this will also fetch legends and
% colorbars as well
hax = findall(hfig, 'type', 'axes');

% get the original axes units, so we can change and reset these again
% later
origaxunits = get(hax, 'Units');

% change the axes units to cm
set(hax, 'Units', 'centimeters');

% get various position parameters of the axes
if numel(hax) > 1
    %         fsize = cell2mat(get(hax, 'FontSize'));
    ti = cell2mat(get(hax,'TightInset'));
    pos = cell2mat(get(hax, 'Position'));
else
    %         fsize = get(hax, 'FontSize');
    ti = get(hax,'TightInset');
    pos = get(hax, 'Position');
end

% ensure very tiny border so outer box always appears
ti(ti < 0.1) = 0.15;

% we will check if any 3d axes are zoomed, to do this we will check if
% they are not being viewed in any of the 2d directions
views2d = [0,90; 0,0; 90,0];

for i = 1:numel(hax)
    
    set(hax(i), 'LooseInset', ti(i,:));
    %         set(hax(i), 'LooseInset', [0,0,0,0]);
    
    % get the current viewing angle of the axes
    [az,el] = view(hax(i));
    
    % determine if the axes are zoomed
    iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
    
    % test if we are viewing in 2d mode or a 3d view
    is2d = all(bsxfun(@eq, [az,el], views2d), 2);
    
    if iszoomed && ~any(is2d)
        error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.')
    end
    
end

% we will move all the axes down and to the left by the amount
% necessary to just show the bottom and leftmost axes and labels etc.
moveleft = min(pos(:,1) - ti(:,1));

movedown = min(pos(:,2) - ti(:,2));

% we will also alter the height and width of the figure to just
% encompass the topmost and rightmost axes and lables
figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);

figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);

% move all the axes
for i = 1:numel(hax)
    
    set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
    
end

origfigunits = get(hfig, 'Units');

set(hfig, 'Units', 'centimeters');

% change the size of the figure
figpos = get(hfig, 'Position');

set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);

% change the size of the paper
set(hfig, 'PaperUnits','centimeters');
set(hfig, 'PaperSize', [figwidth, figheight]);
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition',[0 0 figwidth figheight]);

% reset to original units for axes and figure
if ~iscell(origaxunits)
    origaxunits = {origaxunits};
end

for i = 1:numel(hax)
    set(hax(i), 'Units', origaxunits{i});
end

set(hfig, 'Units', origfigunits);



function scaled = scale_image(imat,scale_zoom)

  oldSize = size(imat);                               % Old image size
  newSize = max(floor(scale_zoom(1:2).*oldSize(1:2)),1);  % New image size
  newX = ((1:newSize(2))-0.5)./scale_zoom(2)+0.5;  % New image pixel X coordinates
  newY = ((1:newSize(1))-0.5)./scale_zoom(1)+0.5;  % New image pixel Y coordinates
  oldClass = class(imat);  % Original image type
  imat = double(imat);      % Convert image to double precision for interpolation
  scaled = interp2(imat,newX,newY(:),'cubic');
  scaled = cast(scaled,oldClass);  % Convert back to original image type



