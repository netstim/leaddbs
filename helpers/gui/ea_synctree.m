function ea_synctree(handles)
ea_busyaction('on',handles.atlasselect,'atlcontrol');

import com.mathworks.mwswing.checkboxtree.*

h=getappdata(handles.atlasselect,'h');
jtree=getappdata(handles.atlasselect,'jtree');
% will sync tree and surfaces based on toggle buttons
sels=ea_storeupdatemodel(jtree,h);
for branch=1:length(sels.branches)
    for leaf=1:length(sels.leaves{branch})
        for side=1:length(sels.sides{branch}{leaf})
            sidec=getsidec(length(sels.sides{branch}{leaf}),side);
            [ixs,ixt]=ea_getsubindex(h.sgsub{branch}{leaf},sidec,h.atlassurfs,h.togglebuttons);
            switch h.togglebuttons(ixt).State
                
                case 'on'
                    
                    set(h.sgsubside{branch}{leaf}{side},'SelectionState',SelectionState.SELECTED)
                    % check if siblings are also selected
                    mixed=0;
                    for chil=1:length(h.sgsubside{branch}{leaf})
                        if ~strcmp(h.sgsubside{branch}{leaf}{chil}.getSelectionState,'selected')
                            mixed=1;
                            break
                        end
                    end
                    % now set parent node accordingly
                    if ~mixed
                        set(h.sgsub{branch}{leaf},'SelectionState',SelectionState.SELECTED)
                    else
                        set(h.sgsub{branch}{leaf},'SelectionState',SelectionState.MIXED)
                    end
                    if ~strcmp(h.atlassurfs(ixs).Visible,'on'); % also make sure surface is right
                        h.atlassurfs(ixs).Visible='on';
                    end
                    
                case 'off'
                    set(h.sgsubside{branch}{leaf}{side},'SelectionState',SelectionState.NOT_SELECTED)
                    % check if siblings are also selected
                    mixed=0;
                    for chil=1:length(h.sgsubside{branch}{leaf})
                        if ~strcmp(h.sgsubside{branch}{leaf}{chil}.getSelectionState,'not selected')
                            mixed=1;
                            break
                        end
                    end
                    % now set parent node accordingly
                    if ~mixed
                        set(h.sgsub{branch}{leaf},'SelectionState',SelectionState.NOT_SELECTED)
                    else
                        set(h.sgsub{branch}{leaf},'SelectionState',SelectionState.MIXED)
                    end
                    if ~strcmp(h.atlassurfs(ixs).Visible,'on'); % also make sure surface is right
                        h.atlassurfs(ixs).Visible='on';
                    end                    
                    if ~strcmp(h.atlassurfs(ixs).Visible,'off'); % also make sure surface is right
                        h.atlassurfs(ixs).Visible='off';
                    end
            end
            
        end
        
    end
end
jtree.updateUI
ea_busyaction('off',handles.atlasselect,'atlcontrol');

function sidec=getsidec(sel,side)

if sel==2
    switch side
        case 1
            sidec='_right';
        case 2
            sidec='_left';
    end
elseif sel==1
    sidec='_midline';
end
