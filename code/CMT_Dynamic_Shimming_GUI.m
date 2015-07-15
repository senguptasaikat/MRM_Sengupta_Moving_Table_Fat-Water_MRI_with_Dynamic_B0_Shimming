function varargout = CMT_Dynamic_Shimming_GUI(varargin)
% COMBI_DYNAMIC_SHIMMING_GUI MATLAB code for COMBI_Dynamic_Shimming_GUI.fig
%      COMBI_DYNAMIC_SHIMMING_GUI, by itself, creates a new COMBI_DYNAMIC_SHIMMING_GUI or raises the existing
%      singleton*.
%
%      H = COMBI_DYNAMIC_SHIMMING_GUI returns the handle to a new COMBI_DYNAMIC_SHIMMING_GUI or the handle to
%      the existing singleton*.
%
%      COMBI_DYNAMIC_SHIMMING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMBI_DYNAMIC_SHIMMING_GUI.M with the given input arguments.
%
%      COMBI_DYNAMIC_SHIMMING_GUI('Property','Value',...) creates a new COMBI_DYNAMIC_SHIMMING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before COMBI_Dynamic_Shimming_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to COMBI_Dynamic_Shimming_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help COMBI_Dynamic_Shimming_GUI

% Last Modified by GUIDE v2.5 08-May-2015 15:16:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @COMBI_Dynamic_Shimming_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @COMBI_Dynamic_Shimming_GUI_OutputFcn, ...
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


% --- Executes just before COMBI_Dynamic_Shimming_GUI is made visible.
function COMBI_Dynamic_Shimming_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to COMBI_Dynamic_Shimming_GUI (see VARARGIN)

% Choose default command line output for COMBI_Dynamic_Shimming_GUI
handles.output = hObject;


if nargin<4 || ~isstruct(varargin{1})
      %%% Error !! Expecting outParams
else
      handles.outParams=varargin{1};
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes COMBI_Dynamic_Shimming_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = COMBI_Dynamic_Shimming_GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, ~)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, ~)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, ~, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


CMT_Dynamic_Shimming;


%%%% IF USING SPLINE INTERPOLATIONS :START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

 
load(['../data_output/','pcenter.mat']);
xx = 1: handles.py;

F0 = spline(double(pcenter)',handles.coeffs(1,:),xx);
Z = spline(double(pcenter)',handles.coeffs_1st(1,:),xx);
X = spline(double(pcenter)',handles.coeffs_1st(2,:),xx);
Y = spline(double(pcenter)',handles.coeffs_1st(3,:),xx);


%%%% scale values from -1 to 1 in mT/m to use in load_samples()

F0_max = max(abs(F0));
Z_max = max(abs(Z));
X_max = max(abs(X));
Y_max = max(abs(Y));

fid = fopen('../data_output/Shim_Table_MAXVALS.txt','wt');
     fprintf(fid,'%f \n',F0_max);
     fprintf(fid,'%f \n',Z_max);
     fprintf(fid,'%f \n',X_max);
     fprintf(fid,'%f \n',Y_max);
     fclose(fid);

F0 = F0/F0_max;
Z = Z/Z_max;
X = X/X_max;
Y = Y/Y_max;
Y = -Y;


%%%%%%%%%% WRITE OUT ALL SHIMS IN  FILE in mT/m %%%%%%%%%%%%%%%%%%

 fid = fopen('../data_output/Shim_Table_F0.txt','wt');
     fprintf(fid,'%f \n',F0);
     fclose(fid);
     
     fid = fopen('../data_output/Shim_Table_Z.txt','wt');
     fprintf(fid,'%f \n',Z);
     fclose(fid);
     
     fid = fopen('../data_output/Shim_Table_X.txt','wt');
     fprintf(fid,'%f \n',X);
     fclose(fid);
     
     
     fid = fopen('../data_output/Shim_Table_Y.txt','wt');
     fprintf(fid,'%f \n',Y);
     fclose(fid);


   
%%%% IF USING SPLINE INTERPOLATIONS :STOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  
% --- Executes on slider movement.
function slider1_Callback(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


value = round(get(hObject,'Value'));

set(handles.text14,'String',num2str(value));  

if value ~= 0
    
axes(handles.axes2);
imagesc(squeeze(handles.fieldmap_Thresholded(:,:,value)),[ handles.fmap_min handles.fmap_max]);

if isfield(handles,'A')
axes(handles.axes5);
imagesc(squeeze(handles.A(:,:,value)),[handles.fmap_min handles.fmap_max]);
end

if isfield(handles,'Freq_Offset') 
set(handles.text5,'String',num2str(handles.Freq_Offset(value))); % F0
end

if isfield(handles,'coeffs_1st') 
set(handles.text6,'String',num2str(handles.coeffs_1st(1,value)));  %%% Z
set(handles.text7,'String',num2str(handles.coeffs_1st(2,value)));  %%% X
set(handles.text8,'String',num2str(handles.coeffs_1st(3,value)));  %%% Y
end
     
end

guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, ~, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


value = round(get(hObject,'Value'));
set(handles.text15,'String',num2str(value));  


if value ~= 0
    
axes(handles.axes3);
imagesc(squeeze(handles.fieldmap_Thresholded(value,:,:)),[handles.fmap_min handles.fmap_max]);

if isfield(handles,'A')
axes(handles.axes6);
imagesc(squeeze(handles.A(value,:,:)),[handles.fmap_min handles.fmap_max]);
end

if isfield(handles,'volmag')
axes(handles.axes9);
imagesc(squeeze(handles.volmag(value,:,:))); 
end
     
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, ~, ~)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, ~, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

value = round(get(hObject,'Value'));
set(handles.text16,'String',num2str(value));  


if value ~= 0
    
axes(handles.axes4);
imagesc(squeeze(handles.fieldmap_Thresholded(:,value,:)),[ handles.fmap_min handles.fmap_max]);

if isfield(handles,'A')
axes(handles.axes7);
imagesc(squeeze(handles.A(:,value,:)),[handles.fmap_min handles.fmap_max]);
end
     
end

colorbar;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, ~, ~)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, ~, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

A =  get(handles.popupmenu1,'Value');
B = get(handles.popupmenu1,'String');
handles.threshold = B{A};


if strcmp(handles.threshold,'Automatic');
    handles.fieldmap_Thresholded  = Threshold(handles.vol(:,:,:,1,1),handles.fieldmap);
end



axes(handles.axes2);
imagesc(squeeze(handles.fieldmap_Thresholded(:,:,100)),[ handles.fmap_min handles.fmap_max]);

axes(handles.axes3);
imagesc(squeeze(handles.fieldmap_Thresholded(60,:,:,1,1)),[handles.fmap_min handles.fmap_max]);

axes(handles.axes4);
imagesc(squeeze(handles.fieldmap_Thresholded(:,50,:,1,1)),[ handles.fmap_min handles.fmap_max]);

axes(handles.axes9);
imagesc(squeeze(handles.vol(60,:,:,1,1)),[ handles.fmap_min handles.fmap_max]);


    
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
clc



function edit1_Callback(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double




% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double



% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double



% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double




% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double




% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, ~, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.fieldmap_Thresholded = handles.outParams.fieldmap.* handles.outParams.mask;
handles.py = handles.outParams.nPy;
handles.dz = handles.outParams.dz/10;
clear handles.outParams;

id = handles.fieldmap_Thresholded < -280;
handles.fieldmap_Thresholded(id) = 0;


handles.fmap_min = -250;
handles.fmap_max = 250;

axes(handles.axes2);
imagesc(squeeze(handles.fieldmap_Thresholded(:,:,50)),[ handles.fmap_min handles.fmap_max]);

axes(handles.axes3);
imagesc(squeeze(handles.fieldmap_Thresholded(60,:,:,1,1)),[handles.fmap_min handles.fmap_max]);

axes(handles.axes4);
imagesc(squeeze(handles.fieldmap_Thresholded(:,50,:,1,1)),[ handles.fmap_min handles.fmap_max]);

% axes(handles.axes9);
% imagesc(squeeze(handles.vol(60,:,:,1,1)),[ handles.fmap_min handles.fmap_max]);

set(handles.slider1,'Max',size(handles.fieldmap_Thresholded,3)); 
set(handles.slider1,'SliderStep',[1/size(handles.fieldmap_Thresholded,3), 1]); 

set(handles.slider2,'Max',size(handles.fieldmap_Thresholded,1)); 
set(handles.slider2,'SliderStep',[1/size(handles.fieldmap_Thresholded,1), 1]); 

set(handles.slider3,'Max',size(handles.fieldmap_Thresholded,2)); 
set(handles.slider3,'SliderStep',[1/size(handles.fieldmap_Thresholded,2), 1]); 
colorbar;


guidata(hObject, handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%% IF USING SPLINE INTERPOLATIONS :START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


load(['../data_output/','pcenter.mat']);
xx = 1: handles.py;

F0 = spline(double(pcenter)',handles.coeffs(1,:),xx);
Z = spline(double(pcenter)',handles.coeffs_1st(1,:),xx);
X = spline(double(pcenter)',handles.coeffs_1st(2,:),xx);
Y = spline(double(pcenter)',handles.coeffs_1st(3,:),xx);


%%%% scale values from -1 to 1 in mT/m to use in load_samples()

F0_max = max(abs(F0));
Z_max = max(abs(Z));
X_max = max(abs(X));
Y_max = max(abs(Y));

fid = fopen('../data_output/Shim_Table_MAXVALS.txt','wt');
     fprintf(fid,'%f \n',F0_max);
     fprintf(fid,'%f \n',Z_max);
     fprintf(fid,'%f \n',X_max);
     fprintf(fid,'%f \n',Y_max);
     fclose(fid);

F0 = F0/F0_max;
Z = Z/Z_max;
X = X/X_max;
Y = Y/Y_max;
% Y = -Y;


%%%%%%%%%% WRITE OUT ALL SHIMS IN  FILE in mT/m %%%%%%%%%%%%%%%%%%

 fid = fopen('../data_output/Shim_Table_F0.txt','wt');
     fprintf(fid,'%f \n',F0);
     fclose(fid);
     
     fid = fopen('../data_output/Shim_Table_Z.txt','wt');
     fprintf(fid,'%f \n',Z);
     fclose(fid);
     
     fid = fopen('../data_output/Shim_Table_X.txt','wt');
     fprintf(fid,'%f \n',X);
     fclose(fid);
     
     
     fid = fopen('../data_output/Shim_Table_Y.txt','wt'); 
     fprintf(fid,'%f \n',Y);
     fclose(fid);

     
    
     
     
       
