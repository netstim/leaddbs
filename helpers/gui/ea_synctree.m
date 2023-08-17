function ea_synctree(handles)
ea_busyaction('on',handles.atlasselect,'atlcontrol');

import com.mathworks.mwswing.checkboxtree.*

h=getappdata(handles.atlasselect,'h');
jtree=getappdata(handles.atlasselect,'jtree');
atlases=getappdata(handles.atlasselect,'atlases');
% will sync tree and surfaces based on toggle buttons
sels=ea_storeupdatemodel(jtree,h);
for branch=1:length(sels.branches)
branchsel=[];
    for leaf=1:length(sels.leaves{branch})
        for side=1:length(sels.sides{branch}{leaf})
            sidec=getsidec(length(sels.sides{branch}{leaf}),side,atlases.types(leaf));
            [ixs,ixt]=ea_getsubindex(h.sgsubfi{branch}{leaf},sidec,h.atlassurfs,h.togglebuttons,0);
            switch h.togglebuttons(ixt).State
                case 'on'
                    set(h.sgsubside{branch}{leaf}{side},'SelectionState',SelectionState.SELECTED)
                                            branchsel(end+1)=1;
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
                    if ~strcmp(h.atlassurfs{ixs}.Visible,'on') % also make sure surface is right
                        h.atlassurfs{ixs}.Visible='on';
                    end

                case 'off'
                    set(h.sgsubside{branch}{leaf}{side},'SelectionState',SelectionState.NOT_SELECTED)
                                                                branchsel(end+1)=0;
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
                    if ~strcmp(h.atlassurfs{ixs}.Visible,'on'); % also make sure surface is right
                        h.atlassurfs{ixs}.Visible='on';
                    end
                    if ~strcmp(h.atlassurfs{ixs}.Visible,'off'); % also make sure surface is right
                        h.atlassurfs{ixs}.Visible='off';
                    end

            end

        end


    end

    % set master set checkbox:
        if (any(branchsel<0)) || (any(branchsel>0) && ~(all(branchsel>0))) % mixed
            set(h.sg{branch},'SelectionState',SelectionState.MIXED)
        elseif all(branchsel>0) % all on
            set(h.sg{branch},'SelectionState',SelectionState.SELECTED)

        elseif all(branchsel==0) % all off
            set(h.sg{branch},'SelectionState',SelectionState.NOT_SELECTED)
        end
end
jtree.updateUI
ea_busyaction('off',handles.atlasselect,'atlcontrol');


function sidec=getsidec(sel,side,type)

if sel==2
    switch side
        case 1
            sidec='_right';
        case 2
            sidec='_left';
    end
elseif sel==1
    switch type
        case 1
            sidec='_right';
        case 2
            sidec='_left';
        case 5
            sidec='_midline';
    end
end
