function ea_dispt(msg)

eatimer=getappdata(0,'eatimer');
if ~isempty(eatimer) % timer has been set before
    elapsed=toc(eatimer);
    fprintf(1,'\b'); % go back one line
    fprintf(1,'%s\n',[' Done. (',num2str(elapsed),' s)']);
end

if ~isempty(msg)
    disp(msg);
    % re set timer for next item.
    tv=tic;
    setappdata(0,'eatimer',tv);
else
    setappdata(0,'eatimer',[]);
end
