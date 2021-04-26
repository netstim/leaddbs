function ea_show_normalization(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn
if ~isfield(options, 'leadprod')
    callingfunction='normalization dbs';
else
    callingfunction=['normalization ', options.leadprod];
end

if options.modality==1
    expdo=1;
    subdir=[options.root,options.patientname,filesep];

else
    expdo=1;
    subdir=[options.root,options.patientname,filesep];
end

disp('Preparing images to show Normalization...');

for export=expdo % if CT, only do 1, if MR, do 1:3.

    try
        switch export
            case 1
                checkf=[options.root,options.prefs.patientdir,filesep,options.prefs.gprenii,',1'];
                checkfn=options.prefs.gprenii;
                outf=['check_',options.prefs.prenii];
                [~, fname] = fileparts(options.prefs.gprenii);
                addstr=['MNI ',upper(options.primarytemplate),' (wireframes) & Preoperative MRI (',fname,')'];
                suff='_pre_tra';
            case 2
                if options.modality==1
                    checkf=[options.root,options.prefs.patientdir,filesep,options.prefs.gtranii,',1'];
                    checkfn=options.prefs.gtranii;
                    [~, fname] = fileparts(options.prefs.gtranii);
                    addstr=['MNI ',upper(options.primarytemplate),' (wireframes) & Postoperative MRI (',fname,')'];
                    suff='_tra';

                elseif options.modality==2
                    checkf=[options.root,options.prefs.patientdir,filesep,'tp_',options.prefs.gctnii,',1'];
                    checkfn=['tp_',options.prefs.gctnii];
                    [~, fname] = fileparts(options.prefs.gctnii);
                    addstr=['MNI ',upper(options.primarytemplate),' (wireframes) & Postoperative (tonemapped) CT (',fname,')'];
                    suff='_ct';
                end
                outf=['check_',options.prefs.tranii];
            case 3
                checkf=[options.root,options.prefs.patientdir,filesep,options.prefs.gcornii,',1'];
                checkfn=options.prefs.gcornii;
                outf=['check_',options.prefs.cornii];
                [~, fname] = fileparts(options.prefs.cornii);
                addstr=['MNI ',upper(options.primarytemplate),' (wireframes) & Postoperative coronal MRI (',fname,')'];
                suff='_cor';
            case 4
                checkf=[options.root,options.prefs.patientdir,filesep,'gl',options.prefs.fa2anat,',1'];
                checkfn=['gl',options.prefs.fa2anat];
                outf=['check_',options.prefs.fa2anat];
                [~, fname] = fileparts(options.prefs.fa2anat);
                addstr=['MNI ',upper(options.primarytemplate),' (wireframes) & FA (',fname,')'];
                suff='_fa';
        end

        w=load([ea_space(options),'wires.mat']);
        pt=ea_load_nii(checkf);

        if ~isequal(size(w.wires),size(pt.img))
            matlabbatch{1}.spm.util.imcalc.input = {[ea_space(options),options.primarytemplate,'.nii'];
                                                     checkf};
            matlabbatch{1}.spm.util.imcalc.output = checkfn;
            matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
            matlabbatch{1}.spm.util.imcalc.expression = 'i2';
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            jobs{1}=matlabbatch;
            spm_jobman('run',jobs);
            clear matlabbatch jobs;
            pt=ea_load_nii(checkf);
        end
        %mni.img(:)=zscore(mni.img(:));

        w.wires=single(w.wires);
        w.wires=w.wires/255;
        w.wires=w.wires.*0.2;
        w.wires=w.wires+0.8;
        switch suff
            case {'_fa','_ct'} % do no windowing for now.
                mni_img=ea_load_nii([ea_space(options),'fa.nii']);
                mni_img.img=(mni_img.img-min(mni_img.img(:)))/(max(mni_img.img(:))-min(mni_img.img(:)));
            otherwise
                pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
                pt.img(pt.img>0.5) = 0.5;
                pt.img=(pt.img-min(pt.img(:)))/(max(pt.img(:)));
                if ~exist('mni_img','var')
                    mni_img=ea_load_nii([ea_space(options),options.primarytemplate,'.nii']);
                    mni_img.img(:)=zscore(mni_img.img(:));
                    mni_img.img=(mni_img.img-min(mni_img.img(:)))/(max(mni_img.img(:))-min(mni_img.img(:)));
                end
        end
        %joint_im=0.5*wires.img+pt.img;
        joint_im=pt.img;
        if strcmp(options.prefs.dev.profile,'se') % do siobhan specific stuff (here don't show wires to enable any template to be used as normalization template)
                        joint_im=joint_im.*w.wires; %shows white wires, if commented out, normalizations are shown without the wires, useful if other templates than the MNI are used to normalize images to
        else
        joint_im=joint_im.*w.wires; %shows white wires, if commented out, normalizations are shown without the wires, useful if other templates than the MNI are used to normalize images to
        end
        %joint_im(w.wires>0.9)=1;
        stand=std(joint_im(:));
        joint_im=joint_im-min(joint_im(:));
            joint_im=joint_im./(max(joint_im(:))-stand);
            joint_im(joint_im>1)=1;
        %joint_im=repmat(joint_im,1,1,1,3);
        %jim=cat(4,mni.img,pt.img,mean(cat(4,mni.img,pt.img),4));
        %ea_imshowpair(jim,options,addstr);

        % ----------------------------------------------------------
        % edited by TH 2016-02-17 to add windowed normalization check
        % ----------------------------------------------------------
        pt.img=single(pt.img);
        mni_img.img=single(mni_img.img);
        joint_im=single(joint_im);

        gridf = [options.root,options.patientname,filesep,'glgrid.mat'];

        if ~exist(gridf, 'file') % 'glgrid.mat' doesn't exist, try to generat it here and then delete the glgrid.nii file
            if exist([options.root,options.patientname,filesep,'glgrid.nii'],'file')
                gridnii=ea_load_nii([options.root,options.patientname,filesep,'glgrid.nii,1']);
                if ~isequal(size(w.wires),size(gridnii.img))
                    matlabbatch{1}.spm.util.imcalc.input = {[ea_space(options),options.primarytemplate,'.nii'];
                                                            [options.root,options.patientname,filesep,'glgrid.nii,1']};
                    matlabbatch{1}.spm.util.imcalc.output = 'glgrid.nii';
                    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
                    matlabbatch{1}.spm.util.imcalc.expression = 'i2';
                    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                    jobs{1}=matlabbatch;
                    spm_jobman('run',jobs);
                    clear matlabbatch jobs;
                    gridnii=ea_load_nii([options.root,options.patientname,filesep,'glgrid.nii,1']);
                end
                gridnii.img=gridnii.img/max(gridnii.img(:));
                gridnii.img=gridnii.img.*255;
                grid=uint8(gridnii.img);
                save(gridf, 'grid');
                delete([options.root,options.patientname,filesep,'glgrid.nii']);
            else
                % Commented this - since by default no grid files are
                % being generated.
                %fprintf('No glgrid.nii file found!\n');
            end
        end

        if exist(gridf, 'file') % 'glgrid.mat' does exist, append the grid image in wim
            g=load(gridf);
            grid=single(g.grid);
            grid=grid-min(grid(:));
            grid=grid./max(grid(:));
            grid=grid.*0.2;
            grid=grid+1;
            grid_im=pt.img.*grid;
            %grid_im(grid>0.7)=1;
            grid_im=grid_im-min(grid_im(:));
            grid_im=grid_im./max(grid_im(:));
            wim = cat(4,pt.img,mni_img.img,joint_im,grid_im);
        else
            wim = cat(4,pt.img,mni_img.img,joint_im);
        end

        clear joint_im pt grid_im
        ea_imshowpair(wim,options,addstr,callingfunction);

%       if strcmp(options.prefs.dev.profile,'se')
%         ;
%       else
    catch ME
        warning(ME.message);
        fprintf('Skip showing normalization of %s\n', checkf);
    end

end

disp('Done.');



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
