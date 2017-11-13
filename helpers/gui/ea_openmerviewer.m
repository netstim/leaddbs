function mercontrolfig=ea_openmerviewer(hobj,ev,varargin)

if nargin>3
    mercontrolfig = ea_mercontrol(varargin{1}, varargin{2}); % supplying resultfig and options
    setappdata(varargin{1}, 'mercontrolfig', mercontrolfig);
else
    mercontrolfig = ea_mercontrol(varargin{1}); % supplying ea_trajectory object
        setappdata(varargin{1}.plotFigureH, 'mercontrolfig', mercontrolfig);
end
%try WinOnTop(mercontrolfig, true); end