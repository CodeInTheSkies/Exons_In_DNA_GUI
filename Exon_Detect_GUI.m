function varargout = Exon_Detect_GUI(varargin)
%***************************************************************************************
% This software accompanies the paper "Filter-Based Methodology for the Location of
% Hot Spots in Proteins and Exons in DNA" by Ramachandran et al., "IEEE Transactions
% on Biomedical Engineering", Volume: 59, Issue: 6, June 2012.
% 
% Copyright (c) 2012 by P. Ramachandran.  All rights reserved.
% 
% This software may be used as is by individual researchers to carry out research but
% citation to the above paper would be expected.  For any other use, the permission
% of the author, Dr. P. Ramachandran (Email: rpara26@gmail.com), must be obtained.
% This is a demonstration software that comes with no warranty expressed or implied.
%***************************************************************************************
% EXON_DETECT_GUI M-file for Exon_Detect_GUI.fig
%      EXON_DETECT_GUI, by itself, creates a new EXON_DETECT_GUI or raises the existing
%      singleton*.
%
%      H = EXON_DETECT_GUI returns the handle to a new EXON_DETECT_GUI or the handle to
%      the existing singleton*.
%
%      EXON_DETECT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXON_DETECT_GUI.M with the given input arguments.
%
%      EXON_DETECT_GUI('Property','Value',...) creates a new EXON_DETECT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Hot_Spot_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Exon_Detect_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Exon_Detect_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Exon_Detect_GUI_OutputFcn, ...
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


%% --- OpeningFcn -- Executes just before Exon_Detect_GUI is made visible.
function Exon_Detect_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Exon_Detect_GUI (see VARARGIN)

% Choose default command line output for Exon_Detect_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Obtaining the name and path of this m-file so that the path can be added
% to the Matlab path
pcfn = mfilename('fullpath');
[nameofthisfile, pathofthisfile] = strtok(fliplr(pcfn), filesep);
nameofthisfile = fliplr(nameofthisfile);   pathofthisfile = fliplr(pathofthisfile);
if pathofthisfile(end) == filesep
    pathofthisfile = pathofthisfile(1:(length(pathofthisfile)-1));
end

% Checking if the path of this file is already in the Matlab search path.
% If it is already in path, then we leave everything as is, i.e., no
% changes to the current search path will be made. In case if it is already
% NOT in path, then we add it to the beginning of the current search path.
% However, when the GUI is closed, the added path will be automatically
% removed. Thus, effectively, there will be no permanent changes made to
% the Matlab search path. We will add it temporarily and then
% remove it when the GUI is closed. No worries! See the "CloseRequest_Fcn"
% callback at the very end of this m-file.
if ~isempty(strfind(lower(path), lower(pathofthisfile)))
    already_in_path = true;
else
    already_in_path = false;
    addpath(pathofthisfile, '-begin')
end
handles.already_in_path = already_in_path;
handles.pathofthisfile = pathofthisfile;

% Assigning default values for the optimization-based filter parameters
handles.R_default = 0.95;
handles.Eta_default = 1e-6;
handles.NN_default = 10000;
handles.K_default = 50;
handles.Tau_default = 0.03;
handles.ShowAll = true;

% Assigning the default values for the "no. of passes" fields
handles.Npass_C_default = 1;
handles.Npass_O_default = 1;

handles.charfrq = 2/3; % Assigning the char. freq. to be the period-3 freq. of N/3.
guidata(hObject, handles);
set(handles.text14,'String',num2str(handles.charfrq));

% Assigning values for the filter specs of the inverse-Chebyshev filter
%
handles.fstop1 = handles.charfrq - 8*1e-3;
set(handles.edit1,'string',num2str(handles.fstop1));
%
handles.fpass1 = handles.charfrq - 3*1e-3;
set(handles.edit2,'string',num2str(handles.fpass1));
%
handles.fpass2 = handles.charfrq + 3*1e-3;
set(handles.edit3,'string',num2str(handles.fpass2));
%
handles.fstop2 = handles.charfrq + 8*1e-3;
set(handles.edit4,'string',num2str(handles.fstop2));
%
handles.astop1 = 30;
set(handles.edit5,'string',num2str(handles.astop1));
%
handles.apass = 1;
set(handles.edit6,'string',num2str(handles.apass));
%
handles.astop2 = 30;
set(handles.edit7,'string',num2str(handles.astop2));

%Assigning values for the filter parameters of the optimization-based
%filter
handles.R    =  handles.R_default;
set(handles.edit9,'string',num2str(handles.R));
%
handles.w0   =  handles.charfrq;
set(handles.edit10,'string',num2str(handles.w0));
%
handles.Eta  =  handles.Eta_default;
set(handles.edit11,'string',num2str(handles.Eta));
%
handles.K    =  handles.K_default;
set(handles.edit12,'string',num2str(handles.K));
%
handles.Tau =  handles.Tau_default;
set(handles.edit13,'string',num2str(handles.Tau));
%
handles.NN   =  handles.NN_default;
set(handles.edit14,'string',num2str(handles.NN));
%
% Assigning the default values to the "no. of passes" fields
handles.Npass_Cheby = handles.Npass_C_default;
set(handles.edit16,'string',num2str(handles.Npass_Cheby));
%
handles.Npass_Opt   = handles.Npass_O_default;
set(handles.edit17,'string',num2str(handles.Npass_Opt));
% ****************************************
% designing the lowpass filter 
% ****************************************
hlp = fdesign.lowpass(.4,.5,1,80);
HdL = design(hlp, 'cheby2', 'MatchExactly', 'passband');
fcofL = coeffs(HdL);
[handles.bL, handles.aL] = sos2tf(fcofL.SOSMatrix, prod(fcofL.ScaleValues));

% Assigning the values for the default EIIP values, as per Cosic '94
handles.default_EIIPvals = [0.1260 0.1335 0.0806 0.1340];

% Updating the handles structure
guidata(hObject, handles);

ini_dir = pwd;

% Populate listbox1
load_listbox1(ini_dir,handles)

set(handles.text59,'String','Ready!');

% UIWAIT makes Exon_Detect_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% --- Outputs from this function are returned to the command line.
function varargout = Exon_Detect_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
 varargout{1} = handles.output;


%% --- Pushbutton2 Callback - "Predict (Inverse-Chebyshev)"
% Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'aIC') && isfield(handles,'Npass_Cheby') && isnumeric(handles.Npass_Cheby) &&...
        ~isnan(handles.Npass_Cheby) && isfield(handles,'selNumcLsqs')
    set(handles.text59,'String','Busy...'); pause(0.001),
    for q = 1:length(handles.selNumcLsqs)
        NumcLse = handles.selNumcLsqs{q};
        plot_ttl = ['Exon locations of gene with accession no. ' handles.selNumcLacc_nos{q}];
        b = handles.bIC;
        a = handles.aIC;
        y_in = NumcLse;
        for zeu = 1:handles.Npass_Cheby
            y = filtfilt(b, a, y_in);
            y_in = y;
        end
        ener = y.^2;
        ener = filtfilt(handles.bL, handles.aL, ener); % Lowpass filter application
        ener = ener/max(ener);
        figure, plot(ener), axis tight 'auto y'
        title({plot_ttl;...
            ['Inv.-Cheby. filter; Seq. Length: ' num2str(length(NumcLse)) '; '...
            'Avg. Mag.: ' num2str(sum(ener)/length(ener)) '.']}),
        xlabel('Nucleotides'), ylabel('Normalized signal power at frequency 2pi/3')
    end
    set(handles.text59,'String','Ready!');
    guidata(hObject, handles);
elseif ~isfield(handles,'aIC')
    beep
    errordlg(['First, design the filter by clicking on "Design"!'],'Bad Input','modal')
    return
else
    beep
    errordlg(['There is some other problem that could not be identified! '...
        'Recheck everything and try again!'],'Bad Input','modal')
    return
end


%%  --- Pushbutton3 Callback - "Freq. response (Inverse-Chebyshev)"
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text59,'String','Busy...'); pause(0.001),
if isfield(handles,'fHd')
    h = fvtool(handles.fHd);
    set(h,'MagnitudeDisplay','Magnitude','DesignMask','off','name','Inverse-Chebyshev Filter');
end
set(handles.text59,'String','Ready!');


%% --- Pushbutton4 Callback - "Filter Info (Inverse-Chebyshev)"
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text59,'String','Busy...'); pause(0.001),
if isfield(handles,'fHd')
    sizesos = size(handles.fHd.sosMatrix);
    f_ord = 2*sizesos(1);
    if datenum(version('-date')) > 732892 %732892 represents the release date of 
        % Matlab version (R2006b), i.e., August 03, 2006
        default_info = info(handles.fHd, 'long');
    else
        default_info = info(handles.fHd);
    end
    filinfo = strvcat(default_info(1:2,:), ['Filter Order : ' num2str(f_ord)], default_info(3:end,:));
    msgbox(filinfo,'Current Filter Information');
end
set(handles.text59,'String','Ready!');


%%  --- Pushbutton5 Callback - "About the Parameters (Optimized)"
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pinfo1 = ['"R" is the pole radius of the initial filter that is used as'...
    ' the starting point for the optimization. The default value assumed is '...
    num2str(handles.R_default) '.'];

pinfo2 = ['"w0" is the center frequency or, in other words, the notch frequency of interest'...
    ' forming the narrow passband of the filter. This corresponds to the period-3 frequency.'];

pinfo3 = ['"Eta" is a small constant. It is used to obtain the interval over which the'...
    ' discrete integration is performed. This interval is defined as the union of'...
    ' [0, w0-Eta] and [w0+Eta, pi]. The default value assumed is ' num2str(handles.Eta_default) '.'];

pinfo4 = ['"K" is the user-specified number of iterations of optimization. The default'...
    ' value assumed is ' num2str(handles.K_default) '.'];

pinfo5 = ['"Tau" is the stability margin. The stability of the filter designed'...
    ' is guaranteed if the poles are inside the unit circle. In order to ensure robust stability'...
    ' "Tau" imposes an additional constraint such that the poles will be forced to lie inside'...
    ' the circle with radius 1-Tau. The default value assumed for "Tau" is '...
    num2str(handles.Tau_default) '.'];

pinfo6 = ['"NN" is the number of divisions of the above interval of discrete integration.'...
    ' The larger this number, the greater would be the number of divisions, resulting in'...
    ' a denser sampling of the interval. The default value assumed is '...
    num2str(handles.NN_default) '.'];

pinfo = strvcat(pinfo1,' ',' ',pinfo2,' ',' ',pinfo3,' ',' ',pinfo4,' ',' ',pinfo5,' ',' ',pinfo6);
msgbox(pinfo,'About the filter parameters','help','replace');


%%  --- Pushbutton6 Callback - "Filter Info (Optimized)"
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'x0') && isfield(handles,'xkp1') &&  isfield(handles,'iter') &&  isfield(handles,'ndel')
    set(handles.text59,'String','Busy...'); pause(0.001),   
%     Assign the filter info from the handles structure
    x0 = handles.x0; xkp1 = handles.xkp1; iter = handles.iter; ndel = handles.ndel;
    
    B2 = [1 0 -1];
    A2 = [1 xkp1(2) xkp1(1)];
    h2 = (1 - xkp1(1))/2;
    
    filinfo = strvcat(...
        'Optimized Filter',...
        '==================',...
        ' ',...
        ['Gain constant : ' num2str(h2)],...
        ['Numerator coefficients : ' num2str(B2)],...
        ['Denominator coefficients : ' num2str(A2)]...
        );
    
    msgbox(filinfo,'Current Filter Information');
    set(handles.text59,'String','Ready!');
end

%%  --- Pushbutton7 Callback - "Freq. responses (Optimized)"
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'x0') && isfield(handles,'xkp1') &&  isfield(handles,'iter') &&  isfield(handles,'ndel')
    set(handles.text59,'String','Busy...'); pause(0.001),
%     Assign the filter info from the handles structure
    x0 = handles.x0; xkp1 = handles.xkp1; iter = handles.iter; ndel = handles.ndel;
    
    B2 = [1 0 -1];
    A2 = [1 xkp1(2) xkp1(1)];
    h2 = (1 - xkp1(1))/2;
    
    hn2 = fvtool(h2*B2,A2);
    set(hn2,'MagnitudeDisplay','Magnitude','DesignMask','off','name',...
        'Optimization Technique - Optimized Filter');
    set(handles.text59,'String','Ready!');
end


%%  --- Pushbutton8 Callback - "Predict (Optimized)"
% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'x0') && isfield(handles,'xkp1') && isfield(handles,'iter') && isfield(handles,'ndel') &&...
        isfield(handles,'Npass_Opt') && isnumeric(handles.Npass_Opt) && ~isnan(handles.Npass_Opt) &&...
        isfield(handles,'selNumcLsqs')
    set(handles.text59,'String','Busy...'); pause(0.001),
    for q = 1:length(handles.selNumcLsqs)
        NumcLse = handles.selNumcLsqs{q};
        plot_ttl = ['Exon locations of gene with accession no. ' handles.selNumcLacc_nos{q}];
        x0 = handles.x0; xkp1 = handles.xkp1; iter = handles.iter; ndel = handles.ndel;
        %
        B2 = [1 0 -1];
        A2 = [1 xkp1(2) xkp1(1)];
        h2 = (1 - xkp1(1))/2;
        %
        %       Apply the filter to the signal
        NumcLseqsset{q} = NumcLse;
        y_in = NumcLse;
        for zeu = 1:handles.Npass_Opt
            y = filtfilt(h2*B2, A2, y_in);
            y_in = y;
        end
        ener = y.^2;
        ener = filtfilt(handles.bL, handles.aL, ener); % LOWPASS FILTER APPLICATION   
        ener = ener/max(ener);
        avgenrval = sum(ener)/length(ener);
        figure, plot(ener); axis tight 'auto y'
        title({plot_ttl;...
            ['Opt. BPN filter; Seq. Length: ' num2str(length(NumcLse)) '; '...
            'Avg. Mag.: ' num2str(avgenrval) '.']}),
        xlabel('Nucleotides'), ylabel('Normalized signal power at frequency 2pi/3')
    end
    set(handles.text59,'String','Ready!');
    guidata(hObject, handles);
elseif ~isfield(handles,'ndel')
    beep
    errordlg('First, design the filter by clicking on "Design"!','Bad Input','modal')
    return
else
    beep
    errordlg(['There is some other problem that I cannot identify! '...
        'Recheck everything and try again!'],'Bad Input','modal')
    return
end


%%  --- Pushbutton9 Callback - "Show only '.txt' files"
% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ShowAll = false;
guidata(hObject, handles);
load_listbox1(pwd,handles)


%%  --- Pushbutton10 Callback - "Show all files"
% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ShowAll = true;
guidata(hObject, handles);
load_listbox1(pwd,handles)


%%  --- Pushbutton11 Callback - "Design Filter - Inverse-Chebyshev"
% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

everything_ok = inv_cheby_is_everything_fine(handles);

if everything_ok
    set(handles.text59,'String','Busy...'); pause(0.001),
    Hd = HSpot_DF_design_file(handles.charfrq,handles.fstop1,handles.fpass1,handles.fpass2,...
        handles.fstop2,handles.astop1,handles.apass,handles.astop2,handles.match);
    handles.fHd = Hd;
    fcof = coeffs(Hd);
    [b,a] = sos2tf(fcof.SOSMatrix, prod(fcof.ScaleValues)); 
    handles.bIC = b;
    handles.aIC = a;
    guidata(hObject, handles);
    set(handles.text59,'String','Designing...'); pause(0.1);
    set(handles.text59,'String','Ready!');
end
        

%%  --- Pushbutton12 Callback - "More Info about No. of Passes - Inverse-Chebyshev"
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pinfo = ['This field represents the no. of passes through the filter. The default value is 1.'...
    ' The greater the no. of passes, the more the effect of the filtering, equivalent to filtering'...
    ' with a higher order filter. For example, if the value in this field is 2, then the protein sequence'...
    ' is filtered twice with a second order filter, which brings in the effect of filtering with a filter'...
    ' of order 4.'];

msgbox(pinfo,'No. of Passes','help','replace');


%%  --- Pushbutton13 Callback - "Design Filter - Optimized"
% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

everything_ok = opt_is_everything_fine(handles);

if everything_ok
    set(handles.text59,'String','Busy...'); pause(0.001),
%     Switch off the "SwitchToMedScale" warning
    s = warning('off', 'optim:quadprog:SwitchToMedScale');

    [x0, xkp1, iter, ndel] = opt_anotch_with_stab_for_GUI(handles.R, handles.w0,...
        handles.Eta, handles.K, handles.Tau, handles.NN);
    
%     Reset all warnings to previous settings
    warning(s);
    
%     Assign the filter info to the handles structure
    handles.x0 = x0; 
    handles.xkp1 = xkp1; 
    handles.iter = iter; 
    handles.ndel = ndel;
    guidata(hObject, handles);
    set(handles.text59,'String','Ready!');
end

    
%%  --- Pushbutton14 Callback - "More Info about No. of Passes - Optimized"
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pinfo = ['This field represents the no. of passes through the filter. The default value is 1.'...
    ' The greater the no. of passes, the more the effect of the filtering, equivalent to filtering'...
    ' with a higher order filter. For example, if the value in this field is 2, then the protein sequence'...
    ' is filtered twice with a second order filter, which brings in the effect of filtering with a filter'...
    ' of order 4.'];

msgbox(pinfo,'No. of Passes','help','replace');


%%  --- Pushbutton25 Callback - "Load the file containing the FASTA DNA sequences"
% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text59,'String','Busy...'); pause(0.001),
handles.mainFASTAfile = handles.mainCHOSENfile;
fid = fopen(handles.mainFASTAfile{1});
accs_nums = {}; cc = 0;
while feof(fid) == 0
    tline = fgetl(fid);
    if tline(1) == '>'
        cc = cc+1;
        accs_nums{cc} = tline(2:end);        
    end
end
status = fclose(fid);
if status ~= 0
    disp('ERROR! File did not close properly!!')
end
handles.accs_nums = accs_nums;
guidata(hObject, handles);
load_listbox4(handles);
set(handles.text94,'String',num2str(length(handles.accs_nums)));
set(handles.text59,'String','Ready!');


%%  --- Pushbutton26 Callback - Read the selected character seqs., check for errors, and store in memory
% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text59,'String','Busy...'); pause(0.001),
fid = fopen(handles.mainFASTAfile{1});
handles.charsqs = {}; char_sq = []; cc = 1; beechmein = false; lencc = 0;
while feof(fid) == 0
    tline = fgetl(fid);
    if isempty(strfind(tline,'>')) && beechmein
        char_sq = strcat(char_sq, tline);
    end
    if ~isempty(strfind(tline,'>'))
        if beechmein
            handles.charsqs{cc} = char_sq;
            cc = cc+1;
            lencc = lencc + length(char_sq);
            beechmein = false;
        end
        if length(handles.selacc_nos) >= cc
            if ~isempty(strfind(tline,handles.selacc_nos{cc}))
                char_sq = [];
                beechmein = true;
            end
        end
    end
end
if beechmein
    handles.charsqs{cc} = char_sq;
    lencc = lencc + length(char_sq);
end
status = fclose(fid);
if status ~= 0
    disp('ERROR! File did not close properly!!')
end
errcellarry = {}; 
handles.errcellarry = errcellarry;
guidata(hObject, handles);
set(handles.text112,'String',num2str(length(handles.charsqs)));
set(handles.text59,'String','Ready!');


%%  --- Pushbutton27 Callback - Convert the character seqs. to numerical seqs.
% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text59,'String','Busy...'); pause(0.001),
if ~isfield(handles,'errcellarry') || ~isempty(handles.errcellarry)
    beep
    errordlg(['Error! Either you have not read the character sequences or they contain ambiguous nucleotides! '...
        'They must be first read and then checked before attempting to convert to numerical sequences!'],...
        'Bad Input','modal')
    set(handles.text59,'String','Ready!');
    return
else    
    if isfield(handles,'charsqs')
        handles.NumcLsqs = {};
        for jj = 1:length(handles.charsqs)
            handles.NumcLsqs{jj} = dnacharNumcL(handles.charsqs{jj});
        end
    else
        beep
        errordlg(['Error! No character sequences found in memory! Read them first!'],'Bad Input','modal')
        set(handles.text59,'String','Ready!');
        return
    end
end
guidata(hObject, handles);
load_listbox5(hObject, eventdata, handles);
set(handles.text59,'String','Ready!');

    
%% --- Listbox1 Callback
% Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
get(handles.figure1,'SelectionType');
if strcmp(get(handles.figure1,'SelectionType'),'open')
    filename = file_list{index_selected};
    if  handles.is_dir1(index_selected)
        cd (filename)
        load_listbox1(pwd,handles)
    end
else
    filename = file_list(index_selected);
    sel_type = handles.is_dir1(index_selected);
    if all(~sel_type)
        handles.mainCHOSENfile = filename;
        guidata(hObject, handles);
    end
end


%% --- Load Listbox1
function load_listbox1(dir_path,handles)
cd (dir_path)

if handles.ShowAll
    dir_struct = dir(dir_path); % If ShowAll is enabled, then get the whole list
else
    aaa = dir(fullfile(dir_path, '*.'));
    bbb = dir(fullfile(dir_path, '*.txt'));
    dir_struct = [aaa; bbb]; % Else, get only the list of directories and '.txt' files
end

handles.file_names1 = {dir_struct.name};
handles.is_dir1 = [dir_struct.isdir];
guidata(handles.figure1,handles);
set(handles.listbox1,'String',handles.file_names1,...
	'Value',1)
set(handles.text2,'String',pwd)


%% --- Listbox1 CreateFcn
% Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Listbox4 Callback
% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4
selacc_indx = get(handles.listbox4,'Value');
accs_nums = get(handles.listbox4,'String');
selacc_nos = accs_nums(selacc_indx);
handles.selacc_nos = selacc_nos;
guidata(hObject, handles);


%% --- Load Listbox4
function load_listbox4(handles)

set(handles.listbox4,'String',handles.accs_nums,...
	'Value',1)


%% --- Listbox4 CreateFcn
% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Listbox5 Callback
% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5

set(handles.text59,'String','Loading...');
seleacc_indx = get(handles.listbox5,'Value');
eaccs_nums = get(handles.listbox5,'String');
handles.selNumcLacc_nos = eaccs_nums(seleacc_indx);
handles.selNumcLsqs = handles.NumcLsqs(seleacc_indx);
handles.selcharsqs = handles.charsqs(seleacc_indx);
guidata(hObject, handles);
set(handles.text59,'String','Ready!');


%% --- Load Listbox5
function load_listbox5(hObject, eventdata, handles)
set(handles.listbox5,'String',handles.selacc_nos,'Value',1)
listbox5_Callback(hObject, eventdata, handles)


%% --- Listbox5 CreateFcn
% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit1 Callback - FStop1 Text Field
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fstop1 = user_entry;
guidata(hObject, handles);

%% --- Edit1 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit2 Callback - FPass1 Text Field
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fpass1 = user_entry;
guidata(hObject, handles);

%% --- Edit2 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit3 Callback - FPass2 Text Field
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fpass2 = user_entry;
guidata(hObject, handles);

%% --- Edit3 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit4 Callback - FStop2 Text Field
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fstop2 = user_entry;
guidata(hObject, handles);

%% --- Edit4 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit5 Callback - AStop1 Text Field
function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.astop1 = user_entry;
guidata(hObject, handles);

%% --- Edit5 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit6 Callback - APass Text Field
function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.apass = user_entry;
guidata(hObject, handles);

%% --- Edit6 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit7 Callback - AStop2 Text Field
function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.astop2 = user_entry;
guidata(hObject, handles);

%% --- Edit7 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit9 Callback - R (Pole Radius) Text Field
function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.R = user_entry;
guidata(hObject, handles);

%% --- Edit9 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit10 Callback - w0 (Center Freq.) Text Field
function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.w0 = user_entry;
guidata(hObject, handles);

%% --- Edit10 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit11 Callback - Eta (Small Constant) Text Field
function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Eta = user_entry;
guidata(hObject, handles);

%% --- Edit11 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit12 Callback - K (No. of Iterations) Text Field
function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.K = user_entry;
guidata(hObject, handles);

%% --- Edit12 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit13 Callback - Tau (Stability Margin) Text Field
function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Tau = user_entry;
guidata(hObject, handles);

%% --- Edit13 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit14 Callback - NN (No. of Divisions) Text Field
function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.NN = user_entry;
guidata(hObject, handles);

%% --- Edit14 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit16 Callback - (No. of Passes) - Inverse-Chebyshev
function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Npass_Cheby = user_entry;
guidata(hObject, handles);


%% --- Edit16 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit17 Callback - (No. of Passes) - Optimized
function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Npass_Opt = user_entry;
guidata(hObject, handles);


%% --- Edit17 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Popup Menu1 Callback
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupmenu1

val = get(hObject,'Value');
switch val
    case 1
        sel = 'passband';
    case 2
        sel = 'stopband';
end
handles.match = sel;
guidata(hObject, handles);

    
%% --- Popup Menu1 CreateFcn
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.match = 'passband';
guidata(hObject, handles);


%% --- Function to check if all inputs are proper for the Inverse-Chebyshev
%  based Technique
function everything_ok = inv_cheby_is_everything_fine(handles)

everything_ok = false;
allfreqs_zero_one = false;  allfreqs_sorted = false;
all_set_to_go = false; allvals_correct = false;

% Checking whether all fields in the handles structure exist
allfields_exist = isfield(handles,'fstop1') && isfield(handles,'fpass1') &&...
    isfield(handles,'fpass2') && isfield(handles,'fstop2') && isfield(handles,'astop1') &&...
    isfield(handles,'apass') && isfield(handles,'astop2');

% Checking if all values have been correctly entered
if allfields_exist
    allvals_correct = ~isnan(handles.fstop1) && ~isnan(handles.fpass1) && ~isnan(handles.fpass2) &&...
        ~isnan(handles.fstop2) && ~isnan(handles.astop1) && ~isnan(handles.apass) &&...
        ~isnan(handles.astop2);
end

% If all freq. values have been entered then checking if they are all in
% the correct range (between 0 and 1) and in the right
% order(fstop1 < fpass1 < fpass2 < fstop2)
if allfields_exist
    freqvec = [handles.fstop1 handles.fpass1 handles.fpass2 handles.fstop2];
    allfreqs_zero_one = (freqvec > 0) & (freqvec < 1);
    allfreqs_zero_one = all(allfreqs_zero_one);
    allfreqs_sorted = issorted(freqvec);
end

all_set_to_go = allfields_exist && allvals_correct && allfreqs_zero_one && allfreqs_sorted;

if ~(all_set_to_go) && isfield(handles,'charfrq')
    beep
    errordlg(['Error in the input values! All values must be correctly entered and all frequencies '...
        'must be in the open interval (0, 1) (normalized) and must obey the order '...
        'FStop1 < FPass1 < FPass2 < FStop2.'],'Bad Input','modal')
    return
elseif ~(all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg(['Something wrong with the preassigned char. freq. OR correct/enter input values! '...
        'All values must be correctly entered and all frequencies '...
        'must be in the open interval (0, 1) (normalized) and must obey the order '...
        'FStop1 < FPass1 < FPass2 < FStop2.'],'Bad Input','modal')
    return
elseif (all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg('Something wrong with the preassigned char. freq.','Bad Input','modal')
    return
else
    everything_ok = true;
end


%% --- Function to check if all inputs are proper for the Optimization
%  based Technique
function everything_ok = opt_is_everything_fine(handles)

everything_ok = false;
all_set_to_go = false; allvals_correct = false;
R_ok = false;  w0_ok = false; Tau_ok = false;

% Checking whether all fields in the handles structure exist
allfields_exist = isfield(handles,'R') && isfield(handles,'w0') &&...
    isfield(handles,'Eta') && isfield(handles,'K') && isfield(handles,'Tau') &&...
    isfield(handles,'NN');

% Checking if all values have been correctly entered
if allfields_exist
    allvals_correct = ~isnan(handles.R) && ~isnan(handles.w0) && ~isnan(handles.Eta) &&...
        ~isnan(handles.K) && ~isnan(handles.Tau) && ~isnan(handles.NN);
end

% If all parameter values have been entered then checking if they are all in
% the correct range (between 0 and 1) and in the right
% order(fstop1 < fpass1 < fpass2 < fstop2)
if allfields_exist && allvals_correct
    R_ok = (handles.R > 0) & (handles.R < 1);
    w0_ok = (handles.w0 > 0) & (handles.w0 < 1);
    Tau_ok = (handles.Tau > 0) & (handles.Tau < 1);
end

all_set_to_go = allfields_exist && R_ok && w0_ok && Tau_ok;

if ~(all_set_to_go) && isfield(handles,'charfrq')
    beep
    errordlg(['Error in the input values! All values must be correctly entered and R, w0, and Tau '...
        'must be in the open interval (0, 1).'],'Bad Input','modal')
    return
elseif ~(all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg(['Something wrong with the preassigned char. freq. OR correct/enter input values! '...
        'All values must be correctly entered and R,w0, and Tau '...
        'must be in the open interval (0, 1).'],'Bad Input','modal')
    return
elseif (all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg('Something wrong with the preassigned char. freq.','Bad Input','modal')
    return
else
    everything_ok = true;
end


%% --- HSpot_DF_design_file 
% 'charf' ---> inputs the characteristic frequency (the center frequency of
% the filter)
%HSPOT_DF_DESIGN_FILE Returns a discrete-time filter object.
function Hd = HSpot_DF_design_file(charf,Fstop1,Fpass1,Fpass2,Fstop2,Astop1,Apass,Astop2,match)
%
% Chebyshev Type II Bandpass filter designed using FDESIGN.BANDPASS.
% All frequency values are normalized to 1.
% IF desired to hardset values ...
% Fstop1 = charf - 8*1e-3;       % First Stopband Frequency
% Fpass1 = charf - 3*1e-3;       % First Passband Frequency
% Fpass2 = charf + 3*1e-3;       % Second Passband Frequency
% Fstop2 = charf + 8*1e-3;       % Second Stopband Frequency
% Astop1 = 30;                   % First Stopband Attenuation (dB)
% Apass  = 1;                    % Passband Ripple (dB)
% Astop2 = 30;                   % Second Stopband Attenuation (dB)
% match  = 'passband';           % Band to match exactly
% 
% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2);
Hd = design(h, 'cheby2', 'MatchExactly', match);

function [x0, xkp1, iter, ndel] = opt_anotch_with_stab_for_GUI(R,w0,eta,K,Tau,NN)
%% --- Opt. Filter Design File
% This m-file designs a second order antinotch digital filter by using an
% optimization technique
d00 = R^2;
d10 = -cos(pi*w0)*(1+d00);
x0 = [d00; d10];
divi = 1/NN;
dw = divi*pi; % dw is the infinitesimal size
omega = [0:divi:(w0-eta) (w0+eta):divi:1]*pi;
% Performing the optimization using QuadProg
xk = x0; iter = 0;
options = optimset('display','off');
hwai = waitbar(iter,'Performing Optimization ...');
while iter < K
    [H,f,A,b] = omgrad(xk, omega, Tau, dw);
    [del, fval, extflg, outpt] = quadprog(H,f,A,b,[],[],[],[],[],options);
    if extflg ~= 1
        outpt, extflg
    end
    xkp1 = xk + del;
    xk = xkp1;
    ndel = norm(del);
    iter = iter + 1;
    waitbar(iter/K)
end
close(hwai);


function [H,f,A,b] = omgrad(xk, omega, Tau, dw)
%% --- Function to compute gradient, H, f, A, and b for using in QuadProg
z = exp(1i*omega);
d0k = xk(1);  d1k = xk(2); 
gk1 = (z.^-2 - 1).*(1 + d1k*z.^-1 + z.^-2).*0.5./(1 + d1k*z.^-1 + d0k*z.^-2).^2;
gk2 = (d0k - 1).*(1 - z.^-2).*(z.^-1).*0.5./(1 + d1k*z.^-1 + d0k*z.^-2).^2;
gk = [gk1; gk2];
H = 2*real((gk*sqrt(dw))*(gk'*sqrt(dw))); 
Hak = 0.5.*(1 - d0k).*(1 - z.^-2)./(1 + d1k*z.^-1 + d0k*z.^-2);
Hak = Hak(:).';
f = 2*real((gk*sqrt(dw))*(Hak'*sqrt(dw))); % This '2' is part of the formula
% Constraints
Abar = [1 1; 1 -1; -1 0; 1 0];
tauu = [1-Tau; 1-Tau; 1-Tau; -Tau];
A = -Abar;
b = tauu + Abar*xk;


function [num_seq] = dnacharNumcL(char_seq)
%% --- DNACharNumcL Function
% This function receives a character DNA sequence as input, assigns the corresponding EIIP
% values for the nucleotides, and returns the numerical
% sequence(s) as output.
len = length(char_seq);
char_seq = upper(char_seq);
tmpeip = zeros(1,len);
% Default EIIP Values
tmpeip(strfind(char_seq,'A')) = 0.1260;
tmpeip(strfind(char_seq,'T')) = 0.1335;
tmpeip(strfind(char_seq,'G')) = 0.0806;
tmpeip(strfind(char_seq,'C')) = 0.1340;
% If you want to use the Optimized pseudo-EIIP Values (see paper)
% tmpeip(strfind(char_seq,'A')) = 0.1994;
% tmpeip(strfind(char_seq,'T')) = 0.1933;
% tmpeip(strfind(char_seq,'G')) = 0.0123;
% tmpeip(strfind(char_seq,'C')) = 0.0692;
num_seq = tmpeip(:)';


%% --- Figure1_CloseRequestFcn
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if ~handles.already_in_path && ~isempty(strfind(lower(path), lower(handles.pathofthisfile)))
    rmpath(handles.pathofthisfile);
end
delete(hObject);

