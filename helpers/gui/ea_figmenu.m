function out = figmenu(arg,cmd)
%figmenu: Add a Figures menu to non-GUI figures.
%  From the Figures menu you can bring any figure or the Command Window to
%  the front, and maximize, restore or close all non-GUI figures.  It is
%  especially helpful when figures are covered by other figures.
%
%  If there are more than 20 figures, they are split off into submenus in
%  groups of 10.
%
%  Syntax:
%    figmenu on                     Set the default CreateFcn so that new
%                                      figures will have a Figures menu.
%    figmenu off                    Remove the default CreateFcn.
%    figmenu                        Toggle the state between 'on' and
%                                      'off'.
%    figmenu add                    Add a Figures menu to existing non-GUI
%                                      figures.
%    figmenu remove                 Remove all Figures menus.
%    figmenu state                  Return the current state, 'on' or
%                                      'off'.
%    figmenu(fig_handles,'add')     Add a Figures menu to the specified
%                                      figure(s).
%    figmenu(fig_handles,'remove')  Remove the Figures menu from the
%                                      specified figure(s).
%
%  Tip: Put 'figmenu on' in your startup file.

%  Note: This function uses the undocumented underlying JavaFrame of a
%  figure for maximizing and restoring.  That functionality is encapsulated
%  in a try-catch block in case it goes away in a future version of MATLAB.
%  Thanks to Yair Altman for insight and techniques.

%  Note: The MATLAB documentation for uimenu warns about dynamically
%  changing menu items from the uimenu callback, however, it seems to work
%  fine on Windows 7 and Mac OS X.

% Version: 1.1, 29 May 2013
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})

if nargin == 0
	% Toggle state of figmenu.
	if isequal(get(0,'DefaultFigureCreateFcn'),@figmenu)
		set(0,'DefaultFigureCreateFcn',get(0,'FactoryFigureCreateFcn'))
		disp('figmenu off')
	else
		set(0,'DefaultFigureCreateFcn',@figmenu)
		disp('figmenu on')
	end

elseif nargin == 1
	% Handle simple string arguments: on, off, add, remove, state.
	switch arg
		case 'on'
			set(0,'DefaultFigureCreateFcn',@figmenu)
		case 'off'
			set(0,'DefaultFigureCreateFcn',get(0,'FactoryFigureCreateFcn'))
		case 'add'
			% Add figmenu to non-GUI figures that don't already have one.
			all_figs = findobj('Type','figure','Menubar','figure');
			all_figmenus = findobj('Type','uimenu','Tag','figmenu');
			bad_figs = get(all_figmenus,{'Parent'});
			ok_figs = setdiff(all_figs,[bad_figs{:}]);
			for ii = 1:length(ok_figs)
				h = uimenu(ok_figs(ii),'Label','Figures','Callback',...
					{@build_fig_menu,ok_figs(ii)},'Tag','figmenu');
				% Add menu items with accelerators.
				uimenu(h,'Label','Next Figure','Accelerator','K',...
					'Separator','on','Callback',{@next,ok_figs(ii)})
				uimenu(h,'Label','Previous Figure','Accelerator','J',...
					'Callback',{@previous,ok_figs(ii)})
				uimenu(h,'Label','Recall Last Figure','Accelerator','L',...
					'Callback',@recall_last)
				uimenu(h,'Label','Command Window','Separator','on',...
					'Accelerator','0','Callback',@(~,~)commandwindow)
			end
		case 'remove'
			% Remove all existing figmenus.
			all_figmenus = findobj('Type','uimenu','Tag','figmenu');
			delete(all_figmenus)
		case 'state'
			if isequal(get(0,'DefaultFigureCreateFcn'),@figmenu)
				out = 'on';
			else
				out = 'off';
			end
		otherwise
			error('Unknown command option.')
	end

elseif nargin == 2
	% figmenu called as CreateFcn or explicitly.
	if isempty(cmd) || strcmpi(cmd,'add')
		% Add figmenu to non-GUI figures that don't already have one.
		all_figs = findobj(arg,'Type','figure','Menubar','figure');
		all_figmenus = findobj('Type','uimenu','Tag','figmenu');
		bad_figs = get(all_figmenus,{'Parent'});
		ok_figs = setdiff(all_figs,[bad_figs{:}]);
		for ii = 1:length(ok_figs)
			h = uimenu(ok_figs(ii),'Label','Figures','Callback',...
				{@build_fig_menu,ok_figs(ii)},'Tag','figmenu');
			% Add menu items with accelerators.
			uimenu(h,'Label','Next Figure','Accelerator','K',...
				'Separator','on','Callback',{@next,ok_figs(ii)})
			uimenu(h,'Label','Previous Figure','Accelerator','J',...
				'Callback',{@previous,ok_figs(ii)})
			uimenu(h,'Label','Recall Last Figure','Accelerator','L',...
				'Callback',@recall_last)
			uimenu(h,'Label','Command Window','Separator','on',...
				'Accelerator','0','Callback',@(~,~)commandwindow)
		end
	elseif strcmpi(cmd,'remove')
		delete(findobj(arg,'Type','uimenu','Tag','figmenu'))
	else
		error('Unknown command option.')
	end
end


%--------------------- Callback functions ---------------------

	function build_fig_menu(h,~,this_fig)
		% The menu is built each time the Figures menu is selected so it is
		% always up-to-date.
		if nargin < 2
			this_fig = arg;
		end
		delete(get(h,'Children'))
		figs = struct( ...
			'handle',num2cell(sort(findall(0,'Type','figure'))),...
			'label','',...
			'isgui',false);
		for i = 1:length(figs)
			% Construct label.
			name = get(figs(i).handle,'Name');
			hasName = ~isempty(name);
			hasNumTitle = strcmp(get(figs(i).handle,'NumberTitle'),'on');
			if hasNumTitle
				if hasName
					figs(i).label = sprintf('Figure %g: %s',...
						figs(i).handle,name);
				else
					figs(i).label = sprintf('Figure %g',figs(i).handle);
				end
			else
				if hasName
					figs(i).label = sprintf('%s',name);
				else
					figs(i).label = sprintf('[Figure %g]',figs(i).handle);
				end
			end

			% Determine if figure is a GUI.
			figs(i).isgui = ~strcmp(get(figs(i).handle,'Menubar'),...
				'figure');
		end

		% Build menu items for non-GUI figures.
		group_threshold = 20;
		group_size = 10;
		idx = find(~[figs.isgui]);
		num_figs = length(idx);
		if num_figs <= group_threshold
			for i = idx
				if figs(i).handle == this_fig
					checked = 'on';
				else
					checked = 'off';
				end
				uimenu(h,'Label',figs(i).label,...
					'Callback',{@raise,figs(i).handle},'Checked',checked)
			end
		else
			num_groups = ceil(num_figs/group_size);
			grp_item = zeros(1,num_groups);
			for grp = 1:num_groups
				grp_item(grp) = uimenu(h,'Label',sprintf('Group %d',grp));
			end
			for i = 1:num_figs
				if figs(idx(i)).handle == this_fig
					checked = 'on';
				else
					checked = 'off';
				end
				uimenu(grp_item(ceil(i/group_size)),...
					'Label',figs(idx(i)).label,...
					'Callback',{@raise,figs(idx(i)).handle},...
					'Checked',checked)
			end
		end

		% Build menu items for GUI figures.
		if any([figs.isgui])
			uimenu(h,'Label','GUI Figures','Enable','off',...
				'Separator','on')
			for i = find([figs.isgui])
				if figs(i).handle == this_fig
					checked = 'on';
				else
					checked = 'off';
				end
				uimenu(h,'Label',figs(i).label,...
					'Callback',{@raise,figs(i).handle},'Checked',checked)
			end
		end

		% Add quick-change menu items.
		uimenu(h,'Label','Next Figure','Accelerator','K',...
			'Separator','on','Callback',{@next,this_fig})
		uimenu(h,'Label','Previous Figure','Accelerator','J',...
			'Callback',{@previous,this_fig})
		uimenu(h,'Label','Recall Last Figure','Accelerator','L',...
			'Callback',@recall_last)

		% Add remaining menu items.
		uimenu(h,'Label','Maximize Non-GUI Figures','Separator','on',...
			'Callback',@maximize)
		uimenu(h,'Label','Restore Non-GUI Figures',...
			'Callback',@restore)
		uimenu(h,'Label','Close Non-GUI Figures',...
			'Callback',@close_figs)
		uimenu(h,'Label','Command Window','Separator','on',...
			'Accelerator','0','Callback',@(~,~)commandwindow)
	end


%---------------------

	function raise(~,~,thisfig)
		% Bring the selected figure to the front.
		figure(thisfig)
	end


%---------------------

	function next(~,~,this_fig)
		% Bring the next figure to the front.
		figsc = get(findobj('Tag','figmenu'),{'Parent'});
		figs = sort([figsc{:}]);
		idx = find(figs == this_fig);
		new_idx = mod(idx,length(figs)) + 1;
		figure(figs(new_idx))
	end


%---------------------

	function previous(~,~,this_fig)
		% Bring the previous figure to the front.
		figsc = get(findobj('Tag','figmenu'),{'Parent'});
		figs = sort([figsc{:}]);
		idx = find(figs == this_fig);
		new_idx = mod(idx - 2,length(figs)) + 1;
		figure(figs(new_idx))
	end


%---------------------

	function recall_last(~,~)
		% Recall the last figure to the front.
		figsc = get(findobj('Tag','figmenu'),{'Parent'});
		figs = [figsc{:}];
		if length(figs) > 1
			figure(figs(2))
		end
	end


%---------------------

	function close_figs(~,~)
		% Close all non-GUI figures.
		fig_handles = findobj('Type','figure','Menubar','figure');
		close(fig_handles)
	end


%---------------------

	function maximize(~,~)
		% Maximizes each non-GUI figure.
		these_figs = findobj('Type','figure','Menubar','figure');
		for i = length(these_figs):-1:1
			try
				java_frame = ea_getJavaFrame(handle(these_figs(i)));
				java_frame.setMaximized(true)
			catch
			end
		end
	end


%---------------------

	function restore(~,~)
		% Un-maximizes each non-GUI figure.
		these_figs = findobj('Type','figure','Menubar','figure');
		for i = length(these_figs):-1:1
			try
				java_frame = ea_getJavaFrame(handle(these_figs(i)));
				java_frame.setMaximized(false)
			catch
			end
		end
	end

end
