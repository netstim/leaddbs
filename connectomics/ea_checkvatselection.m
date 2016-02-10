function [usevat,dimensionality,currentseed,sides]=ea_checkvatselection(handles)
% small helper function that will check whether left, right or both vat
% checkboxes are true.

if (get(handles.rvatcheck,'Value') && strcmp(get(handles.rvatcheck,'Enable'),'on')) && ...
        (get(handles.lvatcheck,'Value') && strcmp(get(handles.lvatcheck,'Enable'),'on'))
    %preparecombinedvat(directory,stim);
    usevat={'right','left'};
    dimensionality=2; % how many ROI.
    currentseed=[1,2];
    sides=[1,2];
    
elseif (get(handles.rvatcheck,'Value') && strcmp(get(handles.rvatcheck,'Enable'),'on')) && ...
        ~(get(handles.lvatcheck,'Value') && strcmp(get(handles.lvatcheck,'Enable'),'on'))
    usevat={'right'};
    dimensionality=1; % how many ROI.
    currentseed=1;
    sides=[1];
elseif ~(get(handles.rvatcheck,'Value') && strcmp(get(handles.rvatcheck,'Enable'),'on')) && ...
        (get(handles.lvatcheck,'Value') && strcmp(get(handles.lvatcheck,'Enable'),'on'))
    usevat={'left'};
    dimensionality=1; % how many ROI.
    currentseed=[1];
    sides=2;
else
    usevat={};
    dimensionality=0;
    currentseed=0;
    sides=[];
end