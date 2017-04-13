function [choice,index] = ea_popupquest(prompt,varargin)

% SYNTAX:
%   ea_popupquest(prompt,varargin)
%       
% INPUTS: 
%       1. prompt is a character string and
%       2. varargin is any number of menu items to be selected from a
%             drop-down menu. Add them as comma separated character strings
%
% EXAMPLE: 
%   ea_popupquest('In what space were your electrode locations determined?',...
%         'Preop','Postop')

% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh (UPMC), Brain Modulation Lab
% Ari Kappel

    d = dialog('Position',[300 300 250 150],'Name','Select One');
    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 90 210 40],...
           'String',prompt);
       
    popup = uicontrol('Parent',d,...
           'Style','popup',...
           'Position',[75 70 100 25],...
           'String',varargin,...
           'Callback',@popup_callback);
       
    btn = uicontrol('Parent',d,...
           'Position',[89 20 70 25],...
           'String','Choose',...
           'Callback','delete(gcf)');
       
    % choice = 'Preop';
       
    % Wait for d to close before running to completion
    uiwait(d);
   
       function popup_callback(popup,callbackdata)
          index = popup.Value;
          popup_items = popup.String;
          % This code uses dot notation to get properties.
          % Dot notation runs in R2014b and later.
          % For R2014a and earlier:
          % idx = get(popup,'Value');
          % popup_items = get(popup,'String');
          choice = char(popup_items(index,:));
       end
end