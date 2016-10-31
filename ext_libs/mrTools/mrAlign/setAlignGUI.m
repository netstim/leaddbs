function setAlignGUI(handles,field,value)
% function value = setAlignGUI(handles,field,value)
%
% djh 7/8/2004

switch field
	case 'nSlices'
		% Reset slice slider
		% This needs to be done when the slice orientation is changed or when
		% the volume is (re-)loaded:
		% - sagittalRadioButton_Callback
		% - coronalRadioButton_Callback
		% - axialRadioButton_Callback
		% - loadVolMenuItem_Callback
		% ALIGN.volSize(ALIGN.sliceOrientation)
		nSlices = value;
		set(handles.sliceSlider,'Min',1);
		set(handles.sliceSlider,'Max',nSlices);
		set(handles.sliceSlider,'sliderStep',[1/(nSlices-1) 10/(nSlices-1)]);
		
	case 'trans'
		% Reset the trans and rot editable text fields
		% This needs to be done when loading or computing an alignment:
		% - initializeMenuItem_Callback
		% - loadAlignMenuItem_Callback
		% - coarseAndFineMenuItem_Callback
		% - coarseMenuItem_Callback
		% - fineMenuItem_Callback
		set(handles.transX,'String',num2str(value(1)));
		set(handles.transY,'String',num2str(value(2)));
		set(handles.transZ,'String',num2str(value(3)));
		
	case 'rot'
		set(handles.rotX,'String',num2str(value(1)));
		set(handles.rotY,'String',num2str(value(2)));
		set(handles.rotZ,'String',num2str(value(3)));

	otherwise
		warndlg('Invalid field for setAlignGUI');
end
