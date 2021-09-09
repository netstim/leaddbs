function [roll_out] = ea_diode_main(options)
%% Determine Orientation for BSCI directed leads from postoperative CT
% has an unsupervised and a supervised version
folder = [options.root options.patientname filesep];

if  options.modality == 1 % check for electrode type and postoperative imaging
    msg = sprintf(['Automatic rotation detection works only for postoperative CT images.']);
    choice = questdlg(msg,'No postoperative CT!','Abort','Abort');
    roll_out = [];
    return
end
%% import CTs and choose which CT to use
if exist([folder options.prefs.ctnii_coregistered],'file') == 2
    ct_reg = ea_load_nii([folder options.prefs.ctnii_coregistered]);
    tmat_reg = ct_reg.mat;
else
    error(['No coregistered CT (',options.prefs.ctnii_coregistered,') found in folder: ' folder])
end
tol=0.0001; % tolerance of (rounding) difference in qform and sform matrices.
if exist([folder options.prefs.rawctnii_unnormalized],'file') == 2
    ct_org = ea_load_nii([folder 'postop_ct.nii']);
    if sum(sum(abs(ct_org.mat-ct_org.private.mat0)))<tol
        tmat_org = ct_org.mat;
    else
        msg = sprintf('Different sForm and qForm matrices found in:\n%s\nPlease select the matrix you want to use.', ct_org.fname);
        choice = questdlg(msg, 'Warning!','sForm','qForm','sForm');
        switch choice
            case 'sForm'
                tmat_org = ct_org.mat;
            case 'qForm'
                tmat_org = ct_org.private.mat0;
        end
    end
    ct = ct_org;
    
else
    msg = sprintf(['No postop_ct.nii found in folder: ' folder '\nScript will run using coregistered rpostop_ct.nii which may lead to inaccurate results.']);
    choice = questdlg(msg,'Warning!','Continue','Abort','Abort');
    switch choice
        case 'Continue'
            disp(['Using rpostop_ct.nii as reference image.'])
            ct = ct_reg;
        case 'Abort'
            error('Aborted by user')
    end
end

%% import transformation matrices for CT coregistration
tmat_reg2org=eye(4); % default.
try
    if strcmp(options.prefs.reco.mancoruse,'postop')
        load([folder 'ea_coregctmethod_applied.mat']);
        switch coregct_method_applied{end}
            case 'ea_coregpostopct_fsl'
                %             tmat_reg2org = dlmread([folder 'anat_t12postop_ct_flirt1.mat']));
                disp(['Warning: Temporary fix to use DiODe algorithm with FLIRT. rpostop_ct is used so results may be slightly less accurate.'])
                ct = ct_reg;
            otherwise
                [tmat_reg2org,ctfname] = ea_getrawct2preniimat(options,1);
                ct=ea_load_nii(ctfname);
        end
    else
        ct = ct_reg;
    end
catch
    reg2org = load([folder 'Postop_CT_2_T1.mat']);
    tmat_reg2org =ea_antsmat2mat(reg2org.AffineTransform_double_3_3,reg2org.fixed);
    tmat_reg2org = inv(tmat_reg2org);
end



sides = {'right','left','3','4','5','6','7','8'};
for side = options.elside
    disp(['Reconstructing rotation of ' sides{side} ' lead!'])
    % import lead information
    load([folder 'ea_reconstruction.mat']); % included in for-loop to make independent ea_save_reconstruction for both sides
    
    %% transform head/tail coordinates from native to image coordinates
    head_native = [reco.native.markers(side).head 1]';
    tail_native = [reco.native.markers(side).tail 1]';
    CTname = find(ct.fname==filesep);
    CTname = ct.fname(CTname(end):end);
    
    if strcmp(CTname,[filesep,options.prefs.rawctnii_unnormalized]) || strcmp(CTname,[filesep,'postop_ct_resliced.nii'])
        % transform rpostop_ct -> postop_ct
        head_mm = (tmat_reg2org) * head_native;
        tail_mm = (tmat_reg2org) * tail_native;
        % transform postop_ct mm -> voxel
        head_vx = inv(tmat_org) * head_mm;
        tail_vx = inv(tmat_org) * tail_mm;
        tmat_vx2mm = tmat_org;
    elseif strcmp(CTname,[filesep,options.prefs.ctnii_coregistered])
        head_mm = head_native;
        tail_mm = tail_native;
        % transfrom rpostop_ct mm -> voxel
        head_vx = inv(tmat_reg) * head_mm;
        tail_vx = inv(tmat_reg) * tail_mm;
        tmat_vx2mm = tmat_reg;
    end
    
    unitvector_mm = (tail_mm - head_mm)/norm(tail_mm - head_mm); % vector along the lead axis with 1mm length
    
    %% launch DiODe for different leads
    if ismember(options.elmodel, {'Boston Scientific Vercise Directed', ...
            'St. Jude Directed 6172 (short)', ...
            'St. Jude Directed 6173 (long)'})
        [roll_y,y,~] = ea_diode_auto(side,ct,head_mm,unitvector_mm,tmat_vx2mm,options.elspec);
    elseif ismember(options.elmodel, {'Medtronic B33005', 'Medtronic B33015'})
        [roll_y,y,~] = ea_diode_medtronic(side,ct,head_mm,unitvector_mm,tmat_vx2mm,options.elspec);
    else  % check for electrode type and postoperative imaging
        msg = sprintf(['No Valid Directional Lead Selected!']);
        choice = questdlg(msg,'No Directional Lead!','Abort','Abort');
        roll_out = [];
        return
    end
    
    if ~isempty(y)
        %% transform y to native space and back
        if strcmp(CTname,[filesep 'postop_ct.nii']) || strcmp(CTname,[filesep 'postop_ct_resliced.nii'])
            % transform postop_ct_mm -> rpostop_ct_mm
            y = tmat_reg2org \ y;
        elseif strcmp(CTname,[filesep 'rpostop_ct.nii'])
            y = y;
        end
        
        head = head_native(1:3)';
        tail = tail_native(1:3)';
        y = y(1:3)' - head;
        
        %% Calculate direction of x and y markers
        [xunitv, yunitv] = ea_calcxy(head, tail, y);
        
        y = head + yunitv * (options.elspec.lead_diameter / 2);
        x = head + xunitv * (options.elspec.lead_diameter / 2);
        
        reco.native.markers(side).y = y;
        reco.native.markers(side).x = x;
        
        %% for direct saving into manual reconstruction
        [coords,trajectory,markers]=ea_resolvecoords(reco.native.markers,options);
        ea_save_reconstruction(coords,trajectory,markers,options.elmodel,1,options)
        
        % %% for transfering to ea_manualreconstruction
        yunitv(3) = 0;
        roll_out = rad2deg(atan2(norm(cross([0 1 0],yunitv)),dot([0 1 0],yunitv)));
        if markers(side).y(1) > markers(side).head(1) % negative 90 points to right, positive 90 points to left
            roll_out = - roll_out;
        end
        disp(['Corrected roll angle roll = ' num2str(rad2deg(roll_y)) ' deg, has been converted to orientation angle = ' num2str(roll_out) ' for compatibility with ea_mancorupdatescene.'])
    else
        roll_out = [];
    end
    %% methods dump:
    ea_methods(options,...
        ['Orientation of directional DBS leads was determined using the algorithm published by Dembek et al. 2019 as implemented in Lead-DBS software.'],...
        {'T.A. Dembek, M. Hoevels, A. Hellerbach, A. Horn, J.N. Petry-Schmelzer, J. Borggrefe, J. Wirths, H.S. Dafsari, M.T. Barbe, V. Visser-Vandewalle & H. Treuer (2019). Directional DBS leads show large deviations from their intended implantation orientation. Parkinsonism Relat Disord. 2019 Oct;67:117-121. doi: 10.1016/j.parkreldis.2019.08.017.'});
    
end
