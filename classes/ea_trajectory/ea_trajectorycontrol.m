function varargout = ea_trajectorycontrol(varargin)
% EA_TRAJECTORYCONTROL MATLAB code for ea_trajectorycontrol.fig
%      EA_TRAJECTORYCONTROL, by itself, creates a new EA_TRAJECTORYCONTROL or raises the existing
%      singleton*.
%
%      H = EA_TRAJECTORYCONTROL returns the handle to a new EA_TRAJECTORYCONTROL or the handle to
%      the existing singleton*.
%
%      EA_TRAJECTORYCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_TRAJECTORYCONTROL.M with the given input arguments.
%
%      EA_TRAJECTORYCONTROL('Property','Value',...) creates a new EA_TRAJECTORYCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_trajectorycontrol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_trajectorycontrol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_trajectorycontrol

% Last Modified by GUIDE v2.5 26-Jan-2020 12:09:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_trajectorycontrol_OpeningFcn, ...
    'gui_OutputFcn',  @ea_trajectorycontrol_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ea_trajectorycontrol is made visible.
function ea_trajectorycontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_trajectorycontrol (see VARARGIN)

% Choose default command line output for ea_trajectorycontrol
handles.output = hObject;

set(0,'CurrentFigure',handles.trajectorycontrol);
im=imread([ea_getearoot,'icons',filesep,'logo_lead_dbs.png']);

image(im);
axis off;
axis equal;


% Update handles structure
guidata(hObject, handles);



% UIWAIT makes ea_trajectorycontrol wait for user response (see UIRESUME)
% uiwait(handles.trajectorycontrol);
movegui(hObject,'northwest');



setappdata(handles.trajectorycontrol,'chandles',handles);
obj=varargin{1};
% add menu:
f = uimenu('Label','Tools');
uimenu(f,'Label','Export Plan as Reconstruction...','Callback',{@ea_plan2reconstruction,obj},'Accelerator','E');


h=uimenu(f,'Label','Common DBS targets');

% Horn et al. ACPC based:
ho=uimenu(h,'Label','Horn et al. 2017 NeuroImage, ACPC Coordinates');
ho1=uimenu(ho,'Label','STN, Parkinson''s Disease, Active Contacts (Caire 2013)');
uimenu(ho1,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Caire2013_lh',handles,obj});
uimenu(ho1,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Caire2013_rh',handles,obj});

ho2=uimenu(ho,'Label','GPi, Dystonia, Active Contacts All Leads (Starr 2016)');
uimenu(ho2,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Starr2016_lh',handles,obj});
uimenu(ho2,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Starr2016_rh',handles,obj});

ho3=uimenu(ho,'Label','GPi, Dystonia, Active Contacts Top Responders (Starr 2016)');
uimenu(ho3,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Starr2016tr_lh',handles,obj});
uimenu(ho3,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Starr2016tr_rh',handles,obj});

ho4=uimenu(ho,'Label','VIM, Essential Tremor, Active Contacts (Papavassiliou 2004)');
uimenu(ho4,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Papavassiliou2004_lh',handles,obj});
uimenu(ho4,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Papavassiliou2004_rh',handles,obj});

ho5=uimenu(ho,'Label','SCC, Depression, Standard Contacts (Hamani 2009)');
uimenu(ho5,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Hamani2009_lh',handles,obj});
uimenu(ho5,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Hamani2009_rh',handles,obj});

ho6=uimenu(ho,'Label','SCC, Depression, Active Contacts (Hamani 2009)');
uimenu(ho6,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Hamani2009ac_lh',handles,obj});
uimenu(ho6,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Hamani2009ac_rh',handles,obj});

ho7=uimenu(ho,'Label','ALIC, OCD, Tip (Nuttin 2003), Target (Anderson 2003)');
uimenu(ho7,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Nuttin2003_lh',handles,obj});
uimenu(ho7,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Nuttin2003_rh',handles,obj});

ho8=uimenu(ho,'Label','NAc, OCD, Target Coordinates (Franzini 2010)');
uimenu(ho8,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Franzini2010_lh',handles,obj});
uimenu(ho8,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Franzini2010_rh',handles,obj});

ho9=uimenu(ho,'Label','NAc, Addiction, Target Coordinates (Müller 2009)','Callback',{@ea_getlittarget_horn,'Müller2009',handles,obj});
uimenu(ho9,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Müller2009_lh',handles,obj});
uimenu(ho9,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Müller2009_rh',handles,obj});

ho10=uimenu(ho,'Label','CM/Pv/VOI, Tourette''s Syndrome, Target Coordinates (Ackermans 2011)');
uimenu(ho10,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Ackermans2011_lh',handles,obj});
uimenu(ho10,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Ackermans2011_rh',handles,obj});

ho11=uimenu(ho,'Label','Fornix, Alzheimer''s Disease, Target Coordinates (Ponce 2015)');
uimenu(ho11,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Ponce2015_lh',handles,obj});
uimenu(ho11,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Ponce2015_rh',handles,obj});

ho12=uimenu(ho,'Label','Fornix, Alzheimer''s Disease, Active Contacts (Ponce 2015)');
uimenu(ho12,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn,'Ponce2015ac_lh',handles,obj});
uimenu(ho12,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn,'Ponce2015ac_rh',handles,obj});

% Horn et al. MNI based:
ho_mni=uimenu(h,'Label','Horn et al. 2017 NeuroImage, MNI Coordinates');
ho1_mni=uimenu(ho_mni,'Label','STN, Parkinson''s Disease, Active Contacts (Caire 2013)');
uimenu(ho1_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Caire2013_lh',handles,obj});
uimenu(ho1_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Caire2013_rh',handles,obj});

ho2_mni=uimenu(ho_mni,'Label','GPi, Dystonia, Active Contacts All Leads (Starr 2016)');
uimenu(ho2_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Starr2016_lh',handles,obj});
uimenu(ho2_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Starr2016_rh',handles,obj});

ho3_mni=uimenu(ho_mni,'Label','GPi, Dystonia, Active Contacts Top Responders (Starr 2016)');
uimenu(ho3_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Starr2016tr_lh',handles,obj});
uimenu(ho3_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Starr2016tr_rh',handles,obj});

ho4_mni=uimenu(ho_mni,'Label','VIM, Essential Tremor, Active Contacts (Papavassiliou 2004)');
uimenu(ho4_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Papavassiliou2004_lh',handles,obj});
uimenu(ho4_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Papavassiliou2004_rh',handles,obj});

ho5_mni=uimenu(ho_mni,'Label','SCC, Depression, Standard Contacts (Hamani 2009)');
uimenu(ho5_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Hamani2009_lh',handles,obj});
uimenu(ho5_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Hamani2009_rh',handles,obj});

ho6_mni=uimenu(ho_mni,'Label','SCC, Depression, Active Contacts (Hamani 2009)');
uimenu(ho6_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Hamani2009ac_lh',handles,obj});
uimenu(ho6_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Hamani2009ac_rh',handles,obj});

ho7_mni=uimenu(ho_mni,'Label','ALIC, OCD, Tip (Nuttin 2003), Target (Anderson 2003)');
uimenu(ho7_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Nuttin2003_lh',handles,obj});
uimenu(ho7_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Nuttin2003_rh',handles,obj});

ho8_mni=uimenu(ho_mni,'Label','NAc, OCD, Target Coordinates (Franzini 2010)');
uimenu(ho8_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Franzini2010_lh',handles,obj});
uimenu(ho8_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Franzini2010_rh',handles,obj});

ho9_mni=uimenu(ho_mni,'Label','NAc, Addiction, Target Coordinates (Müller 2009)','Callback',{@ea_getlittarget_horn_mni,'Müller2009',handles,obj});
uimenu(ho9_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Müller2009_lh',handles,obj});
uimenu(ho9_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Müller2009_rh',handles,obj});

ho10_mni=uimenu(ho_mni,'Label','CM/Pv/VOI, Tourette''s Syndrome, Target Coordinates (Ackermans 2011)');
uimenu(ho10_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Ackermans2011_lh',handles,obj});
uimenu(ho10_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Ackermans2011_rh',handles,obj});

ho11_mni=uimenu(ho_mni,'Label','Fornix, Alzheimer''s Disease, Target Coordinates (Ponce 2015)');
uimenu(ho11_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Ponce2015_lh',handles,obj});
uimenu(ho11_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Ponce2015_rh',handles,obj});

ho12_mni=uimenu(ho_mni,'Label','Fornix, Alzheimer''s Disease, Active Contacts (Ponce 2015)');
uimenu(ho12_mni,'Label','Left Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Ponce2015ac_lh',handles,obj});
uimenu(ho12_mni,'Label','Right Hemisphere','Callback',{@ea_getlittarget_horn_mni,'Ponce2015ac_rh',handles,obj});

% Hoeflich et al. MNI based
hoe=uimenu(h,'Label','Höflich et al. 2010 NeuroImage');
hoe1=uimenu(hoe,'Label','ALIC, OCD, (Anderson 2003, Chang 2010, Nuttin 2003)');
uimenu(hoe1,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'ALIC_OCD_lh',handles,obj});
uimenu(hoe1,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'ALIC_OCD_rh',handles,obj});

hoe2=uimenu(hoe,'Label','NC/N.Acc/VS/VC, OCD, (Aouizerate 2004/9, Baker 2007, Franzini 2010, Greenberg 2006/10, Haq 2010, Okun 2004/7, Sturm 2003)');
uimenu(hoe2,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'VCVS_OCD_lh',handles,obj});
uimenu(hoe2,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'VCVS_OCD_rh',handles,obj});

hoe3=uimenu(hoe,'Label','STN, OCD, (Mallet 2008)');
uimenu(hoe3,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'STN_OCD_lh',handles,obj});
uimenu(hoe3,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'STN_OCD_rh',handles,obj});

hoe4=uimenu(hoe,'Label','sgACC, TRD, (Guinjoan 2010, Hamani 2009, Holtzheimer 2012, Lozano 2012, Puigdemont 2012)');
uimenu(hoe4,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'sgACC_TRD_lh',handles,obj});
uimenu(hoe4,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'sgACC_TRD_rh',handles,obj});

hoe5=uimenu(hoe,'Label','NAcc/VS/VC, TRD, (Bewernick 2010, Lujan 2012, Malone 2009, Schlaepfer 2008)');
uimenu(hoe5,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'NAcc_TRD_lh',handles,obj});
uimenu(hoe5,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'NAcc_TRD_rh',handles,obj});

hoe6=uimenu(hoe,'Label','Thalamus, GTS, (Ackermans 2011, Bajwa 2007, Houeto 2005, Kaido 2011, Maciunas 2007, Marceglia 2010, Savica 2012, Servello 2008, Vandewalle 1999/2003, Vernaleken 2009)');
uimenu(hoe6,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'Thalamus_GTS_lh',handles,obj});
uimenu(hoe6,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'Thalamus_GTS_rh',handles,obj});

hoe7=uimenu(hoe,'Label','GPi, GTS, (Dehning 2008, Diederich 2005, Dueck 2009, Martinez-Fernandez 2011)');
uimenu(hoe7,'Label','Left Hemisphere','Callback',{@ea_getlittarget_hoeflich,'GPi_GTS_lh',handles,obj});
uimenu(hoe7,'Label','Right Hemisphere','Callback',{@ea_getlittarget_hoeflich,'GPi_GTS_rh',handles,obj});



set(handles.electrode_model_plan,'String',ea_resolve_elspec);

setappdata(obj.plotFigureH,'trajcontrolfig',handles.trajectorycontrol);
setappdata(handles.trajectorycontrol,'obj',obj);
set(handles.trajectorycontrol,'name','Edit Trajectory');
ea_synctrajectoryhandles(handles,obj)


function ea_getlittarget_horn(~,~,code,handles,obj)

switch code
    case 'Caire2013_lh'
        acpctarget=[-12.02,-1.53,1.91];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Caire2013_rh'
        acpctarget=[12.02,-1.53,1.91];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Starr2016_lh'
        acpctarget=[-20.0,5.8,0.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Starr2016_rh'
        acpctarget=[20.0,5.8,0.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Starr2016tr_lh'
        acpctarget=[-19.8,5.6,-0.6];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Starr2016tr_rh'
        acpctarget=[19.8,5.6,-0.6];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Papavassiliou2004_lh'
        acpctarget=[-12.8,-5.7,-0.8];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Papavassiliou2004_rh'
        acpctarget=[12.8,-5.7,-0.8];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Hamani2009_lh'
        acpctarget=[-5.6,34.2,3.0];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Hamani2009_rh'
        acpctarget=[5.6,34.2,3.0];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Hamani2009ac_lh'
        acpctarget=[-6.3,34.0,2.6];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Hamani2009ac_rh'
        acpctarget=[6.3,34.0,2.6];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Nuttin2003_lh'
        acpctarget=[-14.0,6.0,-6.0];
        set(handles.MCP,'value',0); set(handles.AC,'value',1); set(handles.PC,'value',0);
        obj.planRelative=[1,1,1,1,1];
    case 'Nuttin2003_rh'
        acpctarget=[14.0,6.0,-6.0];
        set(handles.MCP,'value',0); set(handles.AC,'value',1); set(handles.PC,'value',0);
        obj.planRelative=[1,1,1,1,1];
    case 'Franzini2010_lh'
        acpctarget=[-3.0,16.0,2.0];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Franzini2010_rh'
        acpctarget=[3.0,16.0,2.0];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Müller2009_lh'
        acpctarget=[-6.5,2.7,4.5];
        set(handles.MCP,'value',0); set(handles.AC,'value',1); set(handles.PC,'value',0);
        obj.planRelative=[1,1,1,1,1];
    case 'Müller2009_rh'
        acpctarget=[6.5,2.7,4.5];
        set(handles.MCP,'value',0); set(handles.AC,'value',1); set(handles.PC,'value',0);
        obj.planRelative=[1,1,1,1,1];
    case 'Ackermans2011_lh'
        acpctarget=[-5.0,4.0,0.0];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Ackermans2011_rh'
        acpctarget=[5.0,4.0,0.0];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Ponce2015_lh'
        acpctarget=[-4.4,9.8,7.2];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Ponce2015_rh'
        acpctarget=[4.4,9.8,7.2];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Ponce2015ac_lh'
        acpctarget=[-5.6,12.0,1.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
    case 'Ponce2015ac_rh'
        acpctarget=[5.6,12.0,1.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,1];
end
set(handles.space,'value',1);
set(handles.right,'value',1); set(handles.left,'value',0);
set(handles.anterior,'value',1); set(handles.posterior,'value',0);
set(handles.ventral,'value',1); set(handles.dorsal,'value',0);


t.target=acpctarget;
t.entry=acpctarget+([15,15,-55].*([-1,1,1].^(1+double(acpctarget(1)>0))));

set(handles.targetX,'String',num2str(t.target(1)));
set(handles.targetY,'String',num2str(t.target(2)));
set(handles.targetZ,'String',num2str(t.target(3)));

set(handles.entryX,'String',num2str(t.entry(1)));
set(handles.entryY,'String',num2str(t.entry(2)));
set(handles.entryZ,'String',num2str(t.entry(3)));

obj.target=t;

function ea_getlittarget_horn_mni(~,~,code,handles,obj)

switch code
    case 'Caire2013_lh'
        mnitarget=[-12.58,-13.41,-5.87];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Caire2013_rh'
        mnitarget=[12.58,-13.41,-5.87];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Starr2016_lh'
        mnitarget=[-22.37,-5.57,-4.97];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Starr2016_rh'
        mnitarget=[22.37,-5.57,-4.97];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Starr2016tr_lh'
        mnitarget=[-22.71,-5.74,-3.46];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Starr2016tr_rh'
        mnitarget=[22.71,-5.74,-3.46];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Papavassiliou2004_lh'
        mnitarget=[-13.05,-18.38,-2.01];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Papavassiliou2004_rh'
        mnitarget=[13.05,-18.38,-2.01];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Hamani2009_lh'
        mnitarget=[-6.98,23.60,-11.74];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Hamani2009_rh'
        mnitarget=[6.98,23.60,-11.74];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Hamani2009ac_lh'
        mnitarget=[-7.73,23.44,-11.20];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Hamani2009ac_rh'
        mnitarget=[7.73,23.44,-11.20];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Nuttin2003_lh'
        mnitarget=[-15.29,8.08,1.57];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Nuttin2003_rh'
        mnitarget=[15.29,8.08,1.57];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Franzini2010_lh'
        mnitarget=[-3.78,5.08,-7.79];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Franzini2010_rh'
        mnitarget=[3.78,5.08,-7.79];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Müller2009_lh'
        mnitarget=[-7.66,3.61,-10.35];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Müller2009_rh'
        mnitarget=[7.66,3.61,-10.35];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Ackermans2011_lh'
        mnitarget=[-5.54,-15.81,-3.25];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Ackermans2011_rh'
        mnitarget=[5.54,-15.81,-3.25];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Ponce2015_lh'
        mnitarget=[-4.94,-1.52,-13.98];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Ponce2015_rh'
        mnitarget=[4.94,-1.52,-13.98];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Ponce2015ac_lh'
        mnitarget=[-7.02,0.81,-6.43];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Ponce2015ac_rh'
        mnitarget=[7.02,0.81,-6.43];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
end
set(handles.space,'value',3);
set(handles.right,'value',1); set(handles.left,'value',0);
set(handles.anterior,'value',1); set(handles.posterior,'value',0);
set(handles.ventral,'value',1); set(handles.dorsal,'value',0);


t.target=mnitarget;
t.entry=mnitarget+([15,15,55].*([-1,1,1].^(1+double(mnitarget(1)>0))));

set(handles.targetX,'String',num2str(t.target(1)));
set(handles.targetY,'String',num2str(t.target(2)));
set(handles.targetZ,'String',num2str(t.target(3)));

set(handles.entryX,'String',num2str(t.entry(1)));
set(handles.entryY,'String',num2str(t.entry(2)));
set(handles.entryZ,'String',num2str(t.entry(3)));

obj.target=t;

function ea_getlittarget_hoeflich(~,~,code,handles,obj)
set(handles.space,'value',3);

switch code
    case 'ALIC_OCD_lh'
        mnitarget=[-12,10.6,-2.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'ALIC_OCD_rh'
        mnitarget=[11.2,10.6,-2.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'VCVS_OCD_lh'
        mnitarget=[-7.7,6.9,-5.3];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'VCVS_OCD_rh'
        mnitarget=[7.7,6.9,-5.3];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'STN_OCD_lh'
        mnitarget=[-10.3,-16.7,-1];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'STN_OCD_rh'
        mnitarget=[10.3,-16.7,-1];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'sgACC_TRD_lh'
        mnitarget=[-4.6,26.1,-8.1];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'sgACC_TRD_rh'
        mnitarget=[5,26,-7.9];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'NAcc_TRD_lh'
        mnitarget=[-7.3,6.2,-4.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'NAcc_TRD_rh'
        mnitarget=[7.3,6.2,-4.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Thalamus_GTS_lh'
        mnitarget=[-6.3,-13.1,-0.3];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'Thalamus_GTS_rh'
        mnitarget=[6.3,-13.1,-0.3];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'GPi_GTS_lh'
        mnitarget=[-19.5,-3,-4.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
    case 'GPi_GTS_rh'
        mnitarget=[19.5,-3,-4.5];
        set(handles.MCP,'value',1); set(handles.AC,'value',0); set(handles.PC,'value',0);
        obj.planRelative=[2,1,1,1,3];
end
set(handles.right,'value',1); set(handles.left,'value',0);
set(handles.anterior,'value',1); set(handles.posterior,'value',0);
set(handles.ventral,'value',1); set(handles.dorsal,'value',0);


t.target=mnitarget;
t.entry=mnitarget+[15,15,+55].*([-1,1,1].^(1+double(mnitarget(1)>0)));

set(handles.targetX,'String',num2str(t.target(1)));
set(handles.targetY,'String',num2str(t.target(2)));
set(handles.targetZ,'String',num2str(t.target(3)));

set(handles.entryX,'String',num2str(t.entry(1)));
set(handles.entryY,'String',num2str(t.entry(2)));
set(handles.entryZ,'String',num2str(t.entry(3)));

obj.target=t;


function ea_plan2reconstruction(~,~,obj)
[fileName, path] = uiputfile(obj.options.subj.recon.recon,'Export Plan as Reconstruction');

elstruct = obj.plan2elstruct;
options = obj.options;

ea_save_reconstruction(elstruct(1).coords_mm,elstruct(1).trajectory,elstruct(1).markers,obj.plan2elstruct_model,1,options,[path, fileName]);


% --- Outputs from this function are returned to the command line.
function varargout = ea_trajectorycontrol_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showPlanning.
function showPlanning_Callback(hObject, eventdata, handles)
% hObject    handle to showPlanning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showPlanning
obj = getappdata(handles.trajectorycontrol,'obj');
obj.showPlanning = get(handles.showPlanning,'Value');
obj.togglestates(1) = get(handles.showPlanning,'Value');
obj.toggleH.State = get(handles.showPlanning,'Value');
ea_synctrajectoryhandles(handles,obj);

function targetX_Callback(hObject, eventdata, handles)
% hObject    handle to targetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetX as text
%        str2double(get(hObject,'String')) returns contents of targetX as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.target(1)=str2double(get(handles.targetX,'String'));


% --- Executes during object creation, after setting all properties.
function targetX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetY_Callback(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtY as text
%        str2double(get(hObject,'String')) returns contents of txtY as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.target(2)=str2double(get(handles.targetY,'String'));

% --- Executes during object creation, after setting all properties.
function txtY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetZ_Callback(hObject, eventdata, handles)
% hObject    handle to targetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetZ as text
%        str2double(get(hObject,'String')) returns contents of targetZ as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.target(3)=str2double(get(handles.targetZ,'String'));

% --- Executes during object creation, after setting all properties.
function targetZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entryX_Callback(hObject, eventdata, handles)
% hObject    handle to entryX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entryX as text
%        str2double(get(hObject,'String')) returns contents of entryX as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.entry(1)=str2double(get(handles.entryX,'String'));

% --- Executes during object creation, after setting all properties.
function entryX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entryX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entryY_Callback(hObject, eventdata, handles)
% hObject    handle to entryY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entryY as text
%        str2double(get(hObject,'String')) returns contents of entryY as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.entry(2)=str2double(get(handles.entryY,'String'));

% --- Executes during object creation, after setting all properties.
function entryY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entryY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entryZ_Callback(hObject, eventdata, handles)
% hObject    handle to entryZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entryZ as text
%        str2double(get(hObject,'String')) returns contents of entryZ as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.entry(3)=str2double(get(handles.entryZ,'String'));

% --- Executes during object creation, after setting all properties.
function entryZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entryZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in color.
function color_Callback(hObject, eventdata, handles)
% hObject    handle to color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c = ea_uisetcolor;
if any(c)
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.color=c;
    ea_synctrajectoryhandles(handles,obj);
end

% --- Executes on button press in AC.
function AC_Callback(hObject, eventdata, handles)
% hObject    handle to AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AC
if get(hObject,'Value')
    set(handles.MCP,'Value',0);
    set(handles.PC,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(1)=1;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in MCP.
function MCP_Callback(hObject, eventdata, handles)
% hObject    handle to MCP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MCP
if get(hObject,'Value')
    set(handles.AC,'Value',0);
    set(handles.PC,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(1)=2;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in PC.
function PC_Callback(hObject, eventdata, handles)
% hObject    handle to PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PC
if get(hObject,'Value')
    set(handles.MCP,'Value',0);
    set(handles.AC,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(1)=3;
else
    set(hObject,'Value',1);
end

% --- Executes on selection change in space.
function space_Callback(hObject, eventdata, handles)
% hObject    handle to space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns space contents as cell array
%        contents{get(hObject,'Value')} returns selected item from space

obj=getappdata(handles.trajectorycontrol,'obj');
obj.planRelative(5)=get(handles.space,'Value');
ea_synctrajectoryhandles(handles,obj);

% --- Executes during object creation, after setting all properties.
function space_CreateFcn(hObject, eventdata, handles)
% hObject    handle to space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of right
if get(hObject,'Value')
    set(handles.left,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(2)=1;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of left
if get(hObject,'Value')
    set(handles.right,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(2)=2;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in anterior.
function anterior_Callback(hObject, eventdata, handles)
% hObject    handle to anterior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anterior
if get(hObject,'Value')
    set(handles.posterior,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(3)=1;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in posterior.
function posterior_Callback(hObject, eventdata, handles)
% hObject    handle to posterior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of posterior
if get(hObject,'Value')
    set(handles.anterior,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(3)=2;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in ventral.
function ventral_Callback(hObject, eventdata, handles)
% hObject    handle to ventral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ventral
if get(hObject,'Value')
    set(handles.dorsal,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(4)=1;
else
    set(hObject,'Value',1);
end

% --- Executes on button press in dorsal.
function dorsal_Callback(hObject, eventdata, handles)
% hObject    handle to dorsal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dorsal
if get(hObject,'Value')
    set(handles.ventral,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(4)=2;
else
    set(hObject,'Value',1);
end


% --- Executes during object creation, after setting all properties.
function targetY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showMacro.
function showMacro_Callback(hObject, eventdata, handles)
% hObject    handle to showMacro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showMacro
obj=getappdata(handles.trajectorycontrol,'obj');
obj.showMacro=get(hObject,'Value');
obj.togglestates(2)=get(handles.showMacro,'Value');
ea_synctrajectoryhandles(handles,obj);

% --- Executes on selection change in plan_electrode_model.
function plan_electrode_model_Callback(hObject, eventdata, handles)
% hObject    handle to plan_electrode_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plan_electrode_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plan_electrode_model
obj=getappdata(handles.trajectorycontrol,'obj');
options.elmodeln = get(handles.plan_electrode_model,'Value');
string_list = get(handles.plan_electrode_model,'String');
obj.elmodel=string_list{options.elmodeln};
ea_synctrajectoryhandles(handles,obj);


% --- Executes during object creation, after setting all properties.
function plan_electrode_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plan_electrode_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showMicro.
function showMicro_Callback(hObject, eventdata, handles)
% hObject    handle to showMicro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showMicro
obj=getappdata(handles.trajectorycontrol,'obj');
obj.showMicro=get(hObject,'Value');
obj.togglestates(3)=get(handles.showMicro,'Value');
ea_synctrajectoryhandles(handles,obj);


% --- Executes on selection change in relateMicro.
function relateMicro_Callback(hObject, eventdata, handles)
% hObject    handle to relateMicro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns relateMicro contents as cell array
%        contents{get(hObject,'Value')} returns selected item from relateMicro
obj=getappdata(handles.trajectorycontrol,'obj');
switch get(hObject,'Value')
    case 1
        obj.relateMicro='macro';
    case 2
        obj.relateMicro='planning';
end
ea_synctrajectoryhandles(handles,obj);


% --- Executes during object creation, after setting all properties.
function relateMicro_CreateFcn(hObject, eventdata, handles)
% hObject    handle to relateMicro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addtraj.
function addtraj_Callback(hObject, eventdata, handles)
% hObject    handle to addtraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
obj=getappdata(handles.trajectorycontrol,'obj');
if obj.hasPlanning % -> This will add a new trajectory unrelated to the present one
    
    
else
    
    obj.target=ea_getstandardtarget(obj.side);
    obj.showPlanning=1;
    obj.hasPlanning=1;
    obj.showPlanning=1;
    ea_synctrajectoryhandles(handles,obj);
    
end


% --- Executes when user attempts to close trajectorycontrol.
function trajectorycontrol_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to trajectorycontrol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
obj=getappdata(handles.trajectorycontrol,'obj');
delete(hObject);

ea_save_electrode(obj);


% --- Executes on selection change in planningappearance.
function planningappearance_Callback(hObject, eventdata, handles)
% hObject    handle to planningappearance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns planningappearance contents as cell array
%        contents{get(hObject,'Value')} returns selected item from planningappearance
obj=getappdata(handles.trajectorycontrol,'obj');
switch get(hObject,'value')
    case 1
        obj.planningAppearance='line';
    case 2
        obj.planningAppearance='electrode';
end
ea_synctrajectoryhandles(handles,obj);


% --- Executes during object creation, after setting all properties.
function planningappearance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to planningappearance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in electrode_model_popup.
function electrode_model_popup_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrode_model_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_model_popup
obj=getappdata(handles.trajectorycontrol,'obj');
options.elmodeln = get(handles.electrode_model_popup,'Value');
string_list = get(handles.electrode_model_popup,'String');
obj.elmodel=string_list{options.elmodeln};
ea_synctrajectoryhandles(handles,obj);

% --- Executes during object creation, after setting all properties.
function electrode_model_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in electrode_model_plan.
function electrode_model_plan_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_model_plan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrode_model_plan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_model_plan
obj=getappdata(handles.trajectorycontrol,'obj');
options.elmodeln = get(handles.electrode_model_plan,'Value');
string_list = get(handles.electrode_model_plan,'String');
obj.plan2elstruct_model=string_list{options.elmodeln};
ea_synctrajectoryhandles(handles,obj);

% --- Executes during object creation, after setting all properties.
function electrode_model_plan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_model_plan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in electrode_relative_plan.
function electrode_relative_plan_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_relative_plan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrode_relative_plan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_relative_plan
obj=getappdata(handles.trajectorycontrol,'obj');
toptions.elmodel=obj.plan2elstruct_model;
toptions=ea_resolve_elspec(toptions);
if toptions.elspec.tipiscontact % tip should be coded as 0, contacts as 1,2,3...
    obj.electrodeRelativeToPlan = get(handles.electrode_relative_plan,'Value');
else
    obj.electrodeRelativeToPlan = get(handles.electrode_relative_plan,'Value')-1;
end
ea_synctrajectoryhandles(handles,obj);

% --- Executes during object creation, after setting all properties.
function electrode_relative_plan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_relative_plan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
