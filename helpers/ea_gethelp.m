function ea_gethelp(seltype,hObject)
% function mediating help calls from GUI objects
hhg=0; % help has been given?
switch seltype
    case 'normal'
        return
    case 'alt'
        switch hObject.Tag
            case 'normalize_checkbox'
                hhg=ea_web('https://leaddbs.gitbooks.io/leaddbs-manual/content/Normalization%20of%20Images.html');
            case 'patdir_choosebox'
                hhg=ea_web('https://leaddbs.gitbooks.io/leaddbs-manual/content/Lead-DBS%20Main%20Window%20and%20Loading%20Images.html');
            case 'recent'
                hhg=ea_web('https://leaddbs.gitbooks.io/leaddbs-manual/content/Lead-DBS%20Main%20Window%20and%20Loading%20Images.html');
            case 'electrode_model_popup'
                
            case 'MRCT'
            case 'coregct_checkbox'
                
            case 'coregctmethod'
            case 'coregctcheck'
            case 'coregmrmethod'
               hhg=ea_web('https://leaddbs.gitbooks.io/leaddbs-manual/content/Normalization%20of%20Images.html');

            case 'normmethod'
                
            case 'normcheck'
                
            case 'doreconstruction'
                
            case 'targetpopup'
                
            case 'maskwindow_txt'
                
            case 'refinelocalization'
                
            case 'include_lead_connectome_subroutine'
                
            case 'openleadconnectome'
                
            case 'atlassetpopup'
                
            case 'vizspacepopup'
                
            case 'writeout2d_checkbox'
                
            case 'render_checkbox'
                
            case 'exportservercheck'
                
            case 'run_button'
                
            case 'exportcode'
                
            case 'openpatientdir'
                
            case 'viewmanual'
                
                hhg=ea_web('https://leaddbs.gitbooks.io/leaddbs-manual/content');
            case 'updatebutn'
                
            case 'dlgroupc'
        
        end
end


if ~hhg
    msgbox('Unfortunately, no help information is available for this item.','Oops!','Help');
end


function hg=ea_web(string)
hg=1;
web(string);
