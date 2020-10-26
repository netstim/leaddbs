function ea_elvisible(hobj,ev,atls,pt,side,onoff,options)

% if(getappdata(hobj.Parent.Parent,'altpressed'))
%
%     eltog=getappdata(hobj.Parent.Parent,'eltog');
%     set(eltog,'State',onoff);
%     for el=1:length(atls)
%         for iside=1:length(options.sides)
%            side=options.sides(iside);
%            try
%                set(atls(el).el_render{side}, 'Visible', onoff);
%            end
%         end
%     end
% else
try
    set(atls, 'Visible', onoff);
catch
    keyboard
end
%end
