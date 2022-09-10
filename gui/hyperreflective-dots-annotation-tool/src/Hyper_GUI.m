function varargout = Hyper_GUI(varargin)
% HYPER_GUI MATLAB code for Hyper_GUI.fig
%      HYPER_GUI, by itself, creates a new HYPER_GUI or raises the existing
%      singleton*.
%
%      H = HYPER_GUI returns the handle to a new HYPER_GUI or the handle to
%      the existing singleton*.
%
%      HYPER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYPER_GUI.M with the given input arguments.
%
%      HYPER_GUI('Property','Value',...) creates a new HYPER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Hyper_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Hyper_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Hyper_GUI

% Last Modified by GUIDE v2.5 12-Jun-2020 16:12:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Hyper_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Hyper_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before Hyper_GUI is made visible.
function Hyper_GUI_OpeningFcn(hObject, eventdata, handles, varargin)

% OK
try
    % -----------  prepare environment -------------------------------
    % Get executable path macOS 
    if isdeployed & ismac
        pathRoot = ctfroot;
        indBar = strfind(pathRoot,'/');
        indApp = strfind(pathRoot,'.app');
        aux = sum(indBar < indApp);
        endPath = indBar(aux)-1;% -1 to avoid getting slash
        pathExe = pathRoot(1:endPath);
    else
        pathExe = pwd;
    end
    % clc;
    % addpath(genpath('Source/'));
    handles.output = hObject;
    
    %--------------------- Layout preparing -----------------------------
    set(handles.figure1,'CurrentAxes',handles.axes2); 
    [img, map, alphachannel] = imread('MU_posH.png'); %insert logo
    image(img, 'AlphaData', alphachannel);
    set(gca,'visible','off','XColor','none','YColor','none'); % hide axis

    set(handles.figure1,'CurrentAxes',handles.axes1);
    set(gca,'visible','off','XColor', 'none','YColor','none'); % hide axis

    setappdata(handles.Load_Image,'path_search',pathExe); % default search path = pwd
    setappdata(handles.Load_Image,'path_or',pathExe); % save absolute path = pwd

    DotSize = 5;
    handles.DotSizeChange.Value = DotSize;
    setappdata(handles.Load_Image,'DotSize',DotSize);
    handles.DotSizeText.String = ['Dot Size: ' num2str(DotSize)];

    handles.Enhance.Enable  ='inactive';
    handles.Mode.Enable = 'inactive';
    handles.Load_Image.Enable = 'inactive';
    handles.Export.Enable = 'inactive';
    handles.Layers.Enable = 'inactive';
    handles.ChangeBscan.Enable = 'off';
    handles.DotSizeChange.Enable = 'off';
    handles.AnnotType.Enable = 'inactive';
%--------------------------------------------------------------------
    guidata(hObject, handles);

    %---------- Initial checks of folders -----------------------
    if ~exist([pathExe '/Sessions'],'dir') % if folder doesn't exist create it
        mkdir([pathExe '/Sessions']);
    end

    if isfile([pathExe '/Sessions/SessionList.mat']) % if folder doesn't exist create it
        load([pathExe '/Sessions/SessionList.mat'],'SessionList');
    else
        SessionList = {};
        save([pathExe '/Sessions/SessionList.mat'],'SessionList');
    end

    setappdata(handles.figure1,'SessionList',SessionList);
  
    jedit = findjobj(handles.Enhance);
    jedit.setBorderPainted(false);

    jedit = findjobj(handles.Load_Image);
    jedit.setBorderPainted(false);

    jedit = findjobj(handles.Exit);
    jedit.setBorderPainted(false);

    jedit = findjobj(handles.Export);
    jedit.setBorderPainted(false);

    jedit = findjobj(handles.Layers);
    jedit.Border = [];

    jedit = findjobj(handles.pegotemobo);
    jedit.Border = [];
    
catch ME
%     errorMessage = sprintf('Path:%s Error in %s at line %d.\n\nErrorMessage:\n%s',...
%         pathExe,ME.stack(1).name,ME.stack(1).line,ME.message);
    errorMessage = sprintf('Path %s Error in %s at line %d.\n\nErrorMessage:\n%s',...
        ctfroot,ME.stack(1).name,ME.stack(1).line,ME.message); 
    uiwait(warndlg(errorMessage));
end

function varargout = Hyper_GUI_OutputFcn(hObject, eventdata, handles)

try
    pathExe = getappdata(handles.Load_Image,'path_or'); % save absolute path = pwd
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;

    % Center GUI
    movegui(gcf,'center');

    % --------  PREGUNTAR CARGAR O CREAR SESIÓN DE ANOTACIÓN ---------------
    answer = questdlg('Seleccione opción de anotación', ...
        'Menu inicial','Load previous session','New session','');
            
    % Prepare structure for msgbox 
    boxStruct.Interpreter = 'tex';
    boxStruct.WindowStyle = 'modal';

    switch answer
        case 'Load previous session'
            % Load Session List
            SessionList = getappdata(handles.figure1,'SessionList');

            % Check if there are previous sessions actually
            if isempty(SessionList)  
                msgbox('\fontsize{12} No previous sessions stored. Closing GUI.','Warning','warn',boxStruct);
                close all;
                return;
            end

            % Ask user to select session
            [indx,tf] = listdlg('ListString',SessionList);

            % Check user selects one session
            if tf == 0
                msgbox('\fontsize{12} No session selected. Closing GUI.','Error','error',boxStruct);
                close all;
                return;
            end

            % Check user doesn't select more than one session
            if length(indx) > 1
                msgbox('\fontsize{12} Two or more sessions selected. Closing GUI','Error','error',boxStruct);
                close all;
                return;
            end


            % Load images data
            Session = SessionList{indx};
            load([pathExe '/Sessions/' Session '/' Session '_data.mat'],'Images');

            % Inform user everything went ok
            msgbox(['\fontsize{12} Session \bf ' Session ' \rm correctly loaded'],'Message','none',boxStruct);

        case 'New session'
            SessionList = getappdata(handles.figure1,'SessionList');

            % Ask user to input new session name
            answer = inputdlg('Enter new session name','New Session');

            % Check if user cancels enters something
            if isempty(answer) 
                msgbox('\fontsize{12} Canceled by user','Warning','warn',boxStruct);
                close all;
                return;
            end

            % Check if clicks OK without entering anything
            if isempty(char(answer))
                msgbox('\fontsize{12} Empty name entered. Closing GUI.','Error','error',boxStruct);
                close all;
                return;
            end

            % Check if duplicate
            if sum(strcmp(answer,SessionList))>0
                msgbox('\fontsize{12} Session name already exists. Closing GUI.','Error','error',boxStruct);
                close all;
                return;
            end

            % Check session name doesn't contain blank spaces
            if sum(ismember(char(answer),' ')) > 0
                msgbox('\fontsize{12} Session name contains blank spaces. Closing GUI.','Error','error',boxStruct);
                close all;
                return;            
            end

            % Update SessionList
            Session = char(answer);
            SessionList{length(SessionList)+1} = Session;
            setappdata(handles.figure1,'SessionList',SessionList);
            save([pathExe '/Sessions/SessionList.mat'],'SessionList');

            % Make folders and files for new session
            mkdir([pathExe '/Sessions/' Session]);
            mkdir([pathExe '/Sessions/' Session '/Export']);

            Images = struct('name',{},'dots',{},'NFL',{},'GCL',{},'IPL',{},...
                'INL',{},'OPL',{},'HFL_ONL',{},'ELM_BM',{},'Total',{},'DDS',{},'Hypo',{},'Drusen',{});

            save([pathExe '/Sessions/' SessionList{end} '/' SessionList{end} '_data.mat'],'Images');

            % Inform user everything went ok
            msgbox(['\fontsize{12} Session \bf' Session ' \rm correctly created'],'Message','none',boxStruct);

        otherwise
            msgbox('\fontsize{12} No option selected. Closing GUI.','Warning','warn',boxStruct);
            close all; 
            return;
    end

    % Set Images and Session name
    setappdata(handles.figure1,'Images',Images);
    setappdata(handles.figure1,'Session',Session);

    % Set annot type to PHR
    handles.AnnotType.Value = 1;
    setappdata(handles.figure1,'AnnotType',1);

    % Update Metadata
    metadata = {['Session: ' Session],'Imagename: ','Format: ',...
        'PatientID: ','Eye: ','Date: '};
    handles.Metadata.String = metadata';

    % Activate Load Image
    handles.Load_Image.Enable = 'on';

catch ME
    errorMessage = sprintf('Error in %s at line %d.\n\nErrorMessage:\n%s',...
        ME.stack(1).name,ME.stack(1).line,ME.message);
    uiwait(warndlg(errorMessage));
end

function Load_Image_Callback(hObject, eventdata, handles)

% Save previous data 
save_data(handles);

% Get current session and search path and ask user to open an image 
Session = getappdata(handles.figure1,'Session');
path_search = getappdata(handles.Load_Image,'path_search'); % get search path
[filename,path_file] = uigetfile('*.*','Select an Image',path_search);  % ask user to select image

% Check that user opens something 
if filename == 0
    return;
end

% Check if it is .vol or .jpg
[~,imagename,ext] = fileparts([path_file filename]); 

% Attempt to read image
if strcmp(ext,'.vol')
%     options = 'visu';
        options = 'nodisp';
    try
        [header, ~, ~, BScans] = openVolFast_edit([path_file filename], options);
    catch error
        msgbox('Error loading .vol');
        return;
    end

    % Get important parameters
    nBscan = double(header.NumBScans);
    PatientID = char(header.PatientID);
    Eye = char(header.ScanPosition);
    Date = datestr(double(header.VisitDate+693960),'dd/mm/yyyy');

    % Get first B-Scan for plotting
    I = BScans(:,:,1).^0.25;

    % Set Max number of B-Scans and Slider step to select one by one
    handles.ChangeBscan.Max = nBscan;
    handles.ChangeBscan.Value = 1;
    
    if nBscan ==1
        handles.ChangeBscan.Enable = 'off';
    else
        handles.ChangeBscan.Enable = 'on';
        handles.ChangeBscan.SliderStep = [1/(nBscan-1) 1/(nBscan-1)];
    end

    % Save B-Scans and important parameters
    setappdata(handles.Load_Image,'BScans',BScans);
    
else
    try
        I = imread([path_file filename]); % load new image
        if size(I,3) ==3
            I = rgb2gray(I); % to grayscale for later on contrast enhancement
        end
        nBscan = 1;
    catch ME
        msgbox('Error loading image');
        return;
    end
    PatientID = 'N/A';
    Eye = 'N/A';
    Date = 'N/A';
    
    % Lock Slider (jpg has only 1-Bscan)
    handles.ChangeBscan.Enable = 'off';
end

% Update BScanNum text
handles.BScanNum.String = ['B-Scan:' '1/' num2str(nBscan)];

% Save new path and BScan index (later updated by guinput)
setappdata(handles.Load_Image,'path_search',path_file);
setappdata(handles.Load_Image,'Bind',1);

% Get Images data 
Images = getappdata(handles.figure1,'Images');


% Check if its a known image (all data is  already in Images)
try
    [inList,Id] = find(strcmp(imagename,{Images.name})==1);  % check if annotated
catch err 
    inList = []; % check if size(Images) == 0 and thus error 
end

% Create new subject data
if  isempty(inList) % new subject
    
    % Define Image Id and save it so Images
    Id = size(Images,2)+1;
    Images(Id).name = imagename;

    % Set initial image info
    if strcmp(ext,'.vol')
        Images(Id).dots = cell(1,nBscan);
        Images(Id).NFL = zeros(1,nBscan);
        Images(Id).GCL = zeros(1,nBscan);
        Images(Id).IPL = zeros(1,nBscan);
        Images(Id).INL = zeros(1,nBscan);
        Images(Id).OPL = zeros(1,nBscan);
        Images(Id).HFL_ONL = zeros(1,nBscan);
        Images(Id).ELM_BM = zeros(1,nBscan);
        Images(Id).Total = zeros(1,nBscan);
        Images(Id).DDS = zeros(1,nBscan);
        Images(Id).Hypo = zeros(1,nBscan);
        Images(Id).Drusen = zeros(1,nBscan);
    else   %.jpg .png or so     
        Images(Id).dots = cell(1);
        Images(Id).NFL = 0;
        Images(Id).GCL = 0;
        Images(Id).IPL = 0;
        Images(Id).INL = 0;
        Images(Id).OPL = 0;
        Images(Id).HFL_ONL = 0;
        Images(Id).ELM_BM = 0;
        Images(Id).Total = 0;
        Images(Id).DDS = 0;
        Images(Id).Hypo = 0;
        Images(Id).Drusen = 0;
    end
    
    setappdata(handles.figure1,'Images',Images);
end


% Set environment values 
setappdata(handles.figure1,'Id',Id);
setappdata(handles.Load_Image,'I_plot',I);
setappdata(handles.Load_Image,'I_or',I);
setappdata(handles.Layers,'LayerId',1); % set menú to 1
setappdata(handles.Mode,'mode',0);

% Clear axes
set(handles.figure1,'CurrentAxes',handles.axes1);  
cla reset;

% Set GUI interface values
handles.Enhance.Value = 0; % Set enhance button to off (default)
handles.Layers.Value = 1; % Set list to NFL
handles.Mode.Value = 0;
handles.Enhance.Enable = 'on'; % activate button
handles.Export.Enable  ='on';
handles.Mode.Enable = 'on';
handles.Layers.Enable = 'on';
handles.DotSizeChange.Enable = 'on';
handles.AnnotType.Enable = 'on';

% Update Metadata
metadata = {['Session: ' Session],['Imagename: ' imagename],['Format: ' ext],...
    ['PatientID: ' PatientID],['Eye: ' Eye],['Date: ' Date]};
handles.Metadata.String = metadata';

% Update Axes
updateLayout(handles);

function Exit_Callback(~, ~, handles)
save_data(handles);
close all;

function Enhance_Callback(hObject, ~, handles)

I = getappdata(handles.Load_Image,'I_or');

if (hObject.Value ==1)
    setappdata(handles.Load_Image,'I_plot',adapthisteq(I));
else
    setappdata(handles.Load_Image,'I_plot',I);
end

updateLayout(handles);

function Layers_Callback(hObject, ~, handles)
setappdata(handles.Layers,'LayerId',hObject.Value); % set layer number

function Layers_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CreateDot(~, eventdata, handles)
% Click handler
% Check that its in create dot mode
mode = getappdata(handles.Mode,'mode');
if mode == 1
   return; 
end

% Get Images, Id and B-scan index and annotation type
Images = getappdata(handles.figure1,'Images');
Id = getappdata(handles.figure1,'Id');
Bind = getappdata(handles.Load_Image,'Bind');
AnnotType = getappdata(handles.figure1,'AnnotType');

% f = fields(Images);
f = {'NFL','GCL','IPL','INL','OPL','HFL_ONL','ELM_BM','Total',...
    'DDS','Hypo','Drusen'};

%------- Left click (Create Dot) --------------------
if eventdata.Button ==1
    X = eventdata.IntersectionPoint(1);
    Y = eventdata.IntersectionPoint(2);

    switch AnnotType
        case 1 % PHR
            LayerId = getappdata(handles.Layers,'LayerId'); % get layer
            Images(Id).('Total')(Bind) = Images(Id).('Total')(Bind) + 1; % increase total PHR number
        case 2 % DDS
            LayerId = 9;
        case 3 % Hyporreflective
            LayerId = 10;
        case 4 % Drusen
            LayerId = 11;
    end
    
    Images(Id).dots{Bind}(size(Images(Id).dots{Bind},1)+1,:) = [X Y LayerId]; % save dot
    Images(Id).(f{LayerId})(Bind) = Images(Id).(f{LayerId})(Bind) + 1; % increase layer dot number
    
    setappdata(handles.figure1,'Images',Images);
    updateLayout(handles);
end

function DeleteDot(~, eventdata, handles)
% Check if delete mode selected
mode = getappdata(handles.Mode,'mode');
if mode == 0 
   return; 
end

% Get Images, Id and Bind
Images = getappdata(handles.figure1,'Images');
Id = getappdata(handles.figure1,'Id');
Bind = getappdata(handles.Load_Image,'Bind');

f = {'NFL','GCL','IPL','INL','OPL','HFL_ONL','ELM_BM','Total',...
    'DDS','Hypo','Drusen'};

%------ Left click (Erase dot) -------------------    
if eventdata.Button == 1 
    X = eventdata.IntersectionPoint(1);
    Y = eventdata.IntersectionPoint(2);
      
    Dist = sum(abs(Images(Id).dots{Bind}(:,1:2)-repmat([X Y],size(Images(Id).dots{Bind},1),1)),2);
    [~,ind] = min(Dist);
    LayerId = Images(Id).dots{Bind}(ind,3); % get layer
    Images(Id).dots{Bind}(ind,:)=[]; % delete dot
    if LayerId<=7
        Images(Id).('Total')(Bind) = Images(Id).('Total')(Bind) - 1; % decrease total dot number
    end
    Images(Id).(f{LayerId})(Bind) = Images(Id).(f{LayerId})(Bind) - 1; % decrease layer dot number
    
    setappdata(handles.figure1,'Images',Images);
    updateLayout(handles);
end

function updateLayout(handles)
% Get Images, Id and B-scan index data
Images = getappdata(handles.figure1,'Images');
Id = getappdata(handles.figure1,'Id');
Bind = getappdata(handles.Load_Image,'Bind');
DotSize = getappdata(handles.Load_Image,'DotSize');

% Clear figure
cla(handles.axes1); 

% Plot image
im = imshow(getappdata(handles.Load_Image,'I_plot')); % plot image
hold on;

% Plot annotations
if size(Images(Id).dots{Bind},1) > 0 
    % Plot PHR
    mask = Images(Id).dots{Bind}(:,3)<=7;
    PHR = Images(Id).dots{Bind}(mask,:);
    colors = [243 59 79; 250 167 5;22 247 5;5 247 240;247 5 239;247 247 5;59 81 243]./255;
    s1 = scatter(PHR(:,1)',PHR(:,2)',DotSize,colors(PHR(:,3),:),'MarkerFaceAlpha',0.8);
    
    % Plot DDS
    mask = Images(Id).dots{Bind}(:,3)==9;
    DDS = Images(Id).dots{Bind}(mask,:);
    s2 = scatter(DDS(:,1)',DDS(:,2)',DotSize,'xr','MarkerFaceAlpha',0.8);
    
    % Plot Hyporreflective
    mask = Images(Id).dots{Bind}(:,3)==10;
    Hypo = Images(Id).dots{Bind}(mask,:);
    s3 = scatter(Hypo(:,1)',Hypo(:,2)',DotSize,'sg','MarkerFaceAlpha',0.8);
    
    % Drusen
    mask = Images(Id).dots{Bind}(:,3)==11;
    Drusen = Images(Id).dots{Bind}(mask,:);
    s4 = scatter(Drusen(:,1)',Drusen(:,2)',DotSize,'db','MarkerFaceAlpha',0.8);
    
end

% Set callbacks for clicks 
im.ButtonDownFcn = {@CreateDot,handles}; % set callback for clicks
s1.ButtonDownFcn = {@DeleteDot,handles}; % set callback for clicks
s2.ButtonDownFcn = {@DeleteDot,handles}; % set callback for clicks
s3.ButtonDownFcn = {@DeleteDot,handles}; % set callback for clicks
s4.ButtonDownFcn = {@DeleteDot,handles}; % set callback for clicks

function save_data(handles) 

path_or = getappdata(handles.Load_Image,'path_or');
Session = getappdata(handles.figure1,'Session');
Images = getappdata(handles.figure1,'Images');
Id = getappdata(handles.figure1,'Id');

% Save mat 
save([path_or '/Sessions/' Session '/' Session '_data.mat'],'Images');

% Save png
% if ~isempty(Id) % check if it is empty
%     I = getappdata(handles.Load_Image,'I_or');
%     save_image(I,Images,Id,path_or);
% end

function save_image(I,Images,Id,path_or)
%---- create hidden new image to export data ----------
    f = figure('visible','off');
    imshow(I);
    hold on;
%------------- Plot Dots---------------------------
    if size(Images(Id).dots,1) > 0 % plot dots
        colors = [243 59 79; 250 167 5;22 247 5;5 247 240;247 5 239;247 247 5;59 81 243]./255;
        scatter(Images(Id).dots(:,1)',Images(Id).dots(:,2)',30,colors(Images(Id).dots(:,3),:),'x','MarkerFaceAlpha',0.8);
    end
    
    F = getframe(gca);
    Im = frame2im(F);
    imwrite(Im, [path_or '/Export/png/' Images(Id).name '_annot.png'])
  
    close(f);
      
% --- Executes on button press in Export.
function Export_Callback(hObject, eventdata, handles)

save_data(handles);

% Get Images path and B-Scan number
Images = getappdata(handles.figure1,'Images');
path_or = getappdata(handles.Load_Image,'path_or');
Session = getappdata(handles.figure1,'Session');

% Check if no figure annotated
if size(Images,1) == 0  
    return;
end

% Stack data with all image info
Imagename = {};
BScan = [];
NFL = [];
GCL = [];
IPL = [];
INL = [];
OPL = [];
HFL_ONL = [];
ELM_BM = [];
Total = [];
DDS = [];
Hypo = [];
Drusen = [];

for n=1:length(Images)
    nBscans  = length(Images(n).NFL);
    Imagename = [Imagename ; repmat({Images(n).name},nBscans,1)];
    BScan = [BScan ; (1:nBscans)'];
    NFL = [NFL; Images(n).NFL'];
    GCL = [GCL; Images(n).GCL'];    
    IPL = [IPL; Images(n).IPL'];
    INL = [INL; Images(n).INL'];    
    OPL = [OPL ;Images(n).OPL'];    
    HFL_ONL = [HFL_ONL; Images(n).HFL_ONL'];
    ELM_BM = [ELM_BM; Images(n).ELM_BM'];    
    Total = [Total; Images(n).Total'];    
    
    DDS = [DDS;Images(n).DDS'];
    Hypo = [Hypo;Images(n).Hypo'];
    Drusen = [Drusen;Images(n).Drusen'];
end

% Create table 
T_images_all = table(Imagename,BScan,NFL,GCL,IPL,INL,OPL,HFL_ONL,ELM_BM,Total,DDS,Hypo,Drusen,...
        'VariableNames',{'Image','BScan','NFL','GCL','IPL','INL','OPL','HFL_ONL','ELM_BM','Total_PHR','DDS','Hyporreflective','Drusen'});

annoted = sum(T_images_all{:,3:end},2)>0;
T_images_annot = T_images_all; 
% T_images_annot(T_images_annot.Total==0,:) = [];
T_images_annot(~annoted,:) = []; % delete not annoted lines

% Save table
xlsname = [path_or '/Sessions/' Session '/Export/Annotations_' datestr(now,'yyyy_mm_dd_HHMMSS') '.xls'];
writetable(T_images_annot,xlsname,'Sheet','Annotated Images'); % save fil
writetable(T_images_all,xlsname,'Sheet','All Images'); % save fil

% Prepare structure for msgbox 
boxStruct.Interpreter = 'tex';
boxStruct.WindowStyle = 'modal';

msgbox('\fontsize{12} Data correctly exported','Message','none',boxStruct);

% --- Executes on button press in Mode.
function Mode_Callback(hObject, eventdata, handles)

setappdata(handles.Mode,'mode',hObject.Value);

% --- Executes on slider movement.
function ChangeBscan_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeBscan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get and modify B-scan index
Bind = round(get(hObject,'Value'));
setappdata(handles.Load_Image,'Bind',Bind);

% Update text on image
nBscan = handles.ChangeBscan.Max;
handles.BScanNum.String = ['B-Scan:' num2str(Bind) '/' num2str(nBscan)];

% Update shown Image
BScans = getappdata(handles.Load_Image,'BScans');
setappdata(handles.Load_Image,'I_or',BScans(:,:,Bind).^0.25);
setappdata(handles.Load_Image,'I_plot',BScans(:,:,Bind).^0.25);

% Update Layout
updateLayout(handles);

% --- Executes during object creation, after setting all properties.
function ChangeBscan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChangeBscan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function DotSizeChange_Callback(hObject, eventdata, handles)
% Get and modify B-scan index
DotSize = get(hObject,'Value');
setappdata(handles.Load_Image,'DotSize',DotSize);

handles.DotSizeText.String = ['Dot Size:' num2str(round(DotSize))];
% Update Layout
updateLayout(handles);

% --- Executes during object creation, after setting all properties.
function DotSizeChange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DotSizeChange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in AnnotType.
function AnnotType_Callback(hObject, eventdata, handles)
setappdata(handles.figure1,'AnnotType',get(hObject,'Value'));

% --- Executes during object creation, after setting all properties.
function AnnotType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
