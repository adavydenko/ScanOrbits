function varargout = ScanOrbits(varargin)
% ScanOrbits M-file for ScanOrbits.fig
%      ScanOrbits, by itself, creates a new ScanOrbits or raises the existing
%      singleton*.
%
%      H = ScanOrbits returns the handle to a new ScanOrbits or the handle to
%      the existing singleton*.
%
%      ScanOrbits('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ScanOrbits.M with the given input arguments.
%
%      ScanOrbits('Property','Value',...) creates a new ScanOrbits or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ScanOrbits_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ScanOrbits_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help ScanOrbits

% Last Modified by GUIDE v2.5 05-Nov-2012 18:44:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ScanOrbits_OpeningFcn, ...
                   'gui_OutputFcn',  @ScanOrbits_OutputFcn, ...
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

%==========================================================================

% --- Executes just before ScanOrbits is made visible.
function ScanOrbits_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ScanOrbits (see VARARGIN)

% Choose default command line output for ScanOrbits
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ScanOrbits wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%==========================================================================

% --- Outputs from this function are returned to the command line.
function varargout = ScanOrbits_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%==========================================================================

% --- Executes on button press in makeLEButton.
function makeLEButton_Callback(~, ~, handles)
% hObject    handle to makeLEButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    currPointer = get(ScanOrbits, 'Pointer');
    set(ScanOrbits, 'Pointer', 'watch');
    
    a = str2double(get(handles.aEdit, 'String'));
    
    MakeLindbladEnvelop_Axes(handles.LEAxes, a);
%     
%     setappdata(ScanOrbits, 'I', I);
%     setappdata(ScanOrbits, 'H', H);
    
    pointsToScan = line('XData', [], 'YData', [], ...
                'ButtonDownFcn', 'ScanOrbits(''ButtonDown'',gcbo,[],guidata(gcbo))', ...
                'LineStyle', 'None', ...
                'Marker', '.', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
                'UIContextMenu', handles.ScanPointsCMenu);

    donePoints = line('XData', [], 'YData', [], ...
                'LineStyle', 'None', ...
                'Marker', '.', ...
                'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
            
    setappdata(ScanOrbits, 'PointsToScan', pointsToScan);
    setappdata(ScanOrbits, 'ScannedPoints', donePoints);
    
    set(ScanOrbits, 'Pointer', currPointer);
    
    setappdata(ScanOrbits, 'fHandle', 1);
  
%==========================================================================       

% --- Executes on button press in autoLEScanButton.
function autoLEScanButton_Callback(~, ~, handles)
% hObject    handle to autoLEScanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    currPointer = get(ScanOrbits, 'Pointer');
    
    a = str2double(get(handles.aEdit, 'String'));
    I = @(r) (sqrt((r^4)*(a^2+r^2)^(-3/2)));
    H = @(r) ((r^2/2)*(a^2+r^2)^(-3/2)-(a^2+r^2)^(-1/2));
    
    maxStep = 6;%str2double(get(handles.maxScanXEdit, 'String'));
    arraySize = int8((maxStep / 0.05) + 1);
    X = zeros(1, arraySize); Y = zeros(1, arraySize);
   
    k = 0;
    for r=0.0:0.05:maxStep
        k = k+1;
        X(k) = I(r);
        Y(k) = H(r);
    end
    P = polyfit(X, Y, 40);
    
        
    try
        set(ScanOrbits, 'Pointer', 'watch');

        dXScan = str2double(get(handles.scanXStepEdit, 'String'));
        dYScan = str2double(get(handles.scanYStepEdit, 'String'));
        
        maxIVal = str2double(get(handles.maxScanXEdit, 'String'));

        iVal = 0;
        
        %see AddNewPoint
        XData = get(getappdata(ScanOrbits, 'PointsToScan'), 'XData');
        YData = get(getappdata(ScanOrbits, 'PointsToScan'), 'YData');

        while iVal <= maxIVal
            val = polyval(P, iVal);

            hVal = 0;
            while hVal >= val              
                %AddNewPoint do GUI updated => less performance
                %AddNewPoint(-iVal, hVal);
                %AddNewPoint(iVal, hVal);
                XData = [XData iVal];
                YData = [YData hVal];
                XData = [XData -iVal];
                YData = [YData hVal];

                hVal = hVal - dYScan;
            end

            iVal = iVal + dXScan;
        end
        
        %see AddNewPoint
        set(getappdata(ScanOrbits, 'PointsToScan'), 'XData', XData, 'YData', YData);
    
        set(ScanOrbits, 'Pointer', currPointer);
        set(handles.ClearButton, 'Enable', 'on');
        set(handles.Start, 'Enable', 'on');
        set(handles.CalcCoverage, 'Enable', 'on');
        
    catch ME
        set(ScanOrbits, 'Pointer', currPointer);
        set(handles.ClearButton, 'Enable', 'on'); 
        rethrow(ME);
    end

%==========================================================================   
    
% --- Executes on button press in manualLEScanButton.
function manualLEScanButton_Callback(~, ~, handles)
% hObject    handle to manualLEScanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    button = 1;
    while (button == 1)
        [I, H, button] = ginput(1);
        
        if (button == 1)
            AddNewPoint(I, H)   
            
            set(handles.Start, 'Enable', 'on');
            set(handles.CalcCoverage, 'Enable', 'on');
            set(handles.ClearButton, 'Enable', 'on');
        end         
    end

%==========================================================================

function AddNewPoint (I, H)
    XData = get(getappdata(ScanOrbits, 'PointsToScan'), 'XData');
    YData = get(getappdata(ScanOrbits, 'PointsToScan'), 'YData');
    
    XData = [XData I];
    YData = [YData H];
    
    set(getappdata(ScanOrbits, 'PointsToScan'), 'XData', XData, 'YData', YData);

%==========================================================================

% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, ~, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(getappdata(ScanOrbits, 'PointsToScan'), 'XData', [], 'YData', []);
    set(getappdata(ScanOrbits, 'ScannedPoints'), 'XData', [], 'YData', []);
    
    set(handles.ReMarkButton, 'Enable', 'off');
    set(handles.Start, 'Enable', 'off');
    set(handles.CalcCoverage, 'Enable', 'off');
    set(hObject, 'Enable', 'off');

%==========================================================================

% --- Executes on button press in Start.
function Start_Callback(hObject, ~, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    currPointer = get(ScanOrbits, 'Pointer');
    set(ScanOrbits, 'Pointer', 'watch')
    
    try
        DoCalc(handles);
        set(ScanOrbits, 'Pointer', 'Arrow');

        set(ScanOrbits, 'Pointer', currPointer);
        set(hObject, 'Enable', 'off');
    catch ME
        set(ScanOrbits, 'Pointer', currPointer);
        set(hObject, 'Enable', 'off');
        rethrow(ME);
    end
    
%==========================================================================

function DoCalc(handles)

    PointsToScan = getappdata(ScanOrbits, 'PointsToScan');
    ScannedPoints = getappdata(ScanOrbits, 'ScannedPoints');

%     hwb = waitbar(0, 'Идет расчет орбит'); %при большом количестве точек
%     - начинает мерцать. Выключаю до лучших времен.
    
    f_i = getappdata(ScanOrbits, 'fHandle');
    f_start = f_i;
    f_count = length(get(PointsToScan, 'XData'));
    
    while ~isempty(get(PointsToScan, 'XData'))
        
        f_i = f_i + 1;
%         waitbar((f_i - f_start)/f_count, hwb, 'updated titile')
        
        XData = get(PointsToScan, 'XData');
        YData = get(PointsToScan, 'YData');

        XDataScanned = get(ScannedPoints, 'XData');
        YDataScanned = get(ScannedPoints, 'YData');
        
        h = XData(1);
        E = YData(1);
        
        hFigure = figure(f_i);
        set(hFigure, 'Name', num2str(f_i));
        set(0, 'CurrentFigure', hFigure);
        set(hFigure, 'Renderer', 'OpenGL');
        
        hold all  
        isShownOnFigure = DoOnePCalc(handles, h, E, f_i);
        if ~isShownOnFigure
            close(hFigure); %Если график не был использован, то закрываем это окно
        end
        hold off
        
%         figure(hwb) %Переводим окно прогресса на передний план
        
        %Помечаем, что мы рассчитали еще одно значение (точку на диаграмме)
        XDataScanned = [XDataScanned, XData(1)];%i
        XData(1) = [];

        YDataScanned = [YDataScanned, YData(1)];%h
        YData(1) = [];

        set(PointsToScan, 'XData', XData, 'YData', YData);
        set(ScannedPoints, 'XData', XDataScanned, 'YData', YDataScanned);
        
        setappdata(ScanOrbits, 'PointsToScan', PointsToScan);
        setappdata(ScanOrbits, 'ScannedPoints', ScannedPoints);
        
        %Кнопка "Пометить как нерассчитанные" становится теперь доступна
        set(handles.ReMarkButton, 'Enable', 'on');
    end
    
    setappdata(ScanOrbits, 'fHandle', f_i);
%     
%     close (hwb)

%==========================================================================
    
function result = DoOnePCalc(handles, h, E, f_num)
    isFigureUsed = false;   
    
    %Время интегрирования и параметры системы ур-ий
    t = [str2double(get(handles.tStartEdit, 'String')) str2double(get(handles.tEndEdit, 'String'))];
    a = str2double(get(handles.aEdit, 'String'));
    gamma = str2double(get(handles.gammaEdit, 'String'));
    
    %k_r = str2double(get(handles.krEdit, 'String'));
    %k_z = str2double(get(handles.kzEdit, 'String'));
    
    r = GetRho(E, h, a);
    
    %Подготавливаем место для графиков
    
    h_main = subplot(2,4,[1 2 5 6]); xlabel(h_main, 'x'); ylabel(h_main, 'y'); %Основной график движения XY
    h_xVx = subplot(2,4,3); xlabel(h_xVx, 'x'); ylabel(h_xVx, 'Vx'); %Отношение координаты к скоростям
    h_yVy = subplot(2,4,7); xlabel(h_yVy, 'y'); ylabel(h_yVy, 'Vy');
    h_rV = subplot(2,4,8); xlabel(h_rV, 'r'); ylabel(h_rV, 'V');
    h_p_yVy = subplot(2,4,4); xlabel(h_p_yVy, 'Пуанкаре - y'); ylabel(h_p_yVy, 'Пуанкаре - Vy'); %Сечение Пуанкаре при x=0
    
    hold all;
    [~, m] = size(r);
    if m ~= 0
        keepedPointNumber = 0;
        ecapedPointNumber = 0;
        for j = 1:m
             rho = r(j);
             w = 0;
             for w = 0:pi*0.05:2*pi
               [x_, y_, Vx_, Vy_] = getIC_full(rho, w, E, h, a); %Для случая изолированного скопления - getIC_isolated(rho, w, E, h)
               %Указываем функцию, описывающее уравнения движения
               odeFuncHandle = @myODE_full; % | @myODE_isolatedCluster | @my3dODE_full

               %проверим условие вылета, если оно не выполняется - переходим к следующей точке
               
               if isequal(get(handles.chbShowOnlyLost, 'Value'), 1)
                    if ~CheckOutCondition(gamma, a, x_, y_, Vx_, Vy_, f_num)
                        continue;
                    end
               end
                
                isFigureUsed = isFigureUsed || true; %Если нарисовали что-то на графике, то указываем это

                %рассчитываем ли мы одну звезду
                if isequal(get(handles.doGroupCalcChb, 'Value'), 0)
                    %Решаем уравнения движения
                    [t, Y] = solveMyODE45(odeFuncHandle, t, [x_ y_ Vx_ Vy_], gamma, a);

                    p_y = [];
                    p_Vy = [];                
                    for c = 1:(length(Y(:,1))-1)
                        if (Y(c,1)<=0 && Y(c+1,1)>0) || (Y(c,1)>=0 && Y(c+1,1)<0)
                            p_y = [p_y interp1([Y(c,1) Y(c+1,1)],[Y(c,2) Y(c+1,2)],0)];
                            p_Vy = [p_Vy interp1([Y(c,1) Y(c+1,1)],[Y(c,4) Y(c+1,4)],0)];
                        end
                    end 
                    
                    if isequal(get(handles.chbPrintCloud, 'Value'), 0)
                        set(gcf,'CurrentAxes',h_main); hold all; plot(Y(:,1), Y(:,2)); hold off;
                        set(gcf,'CurrentAxes',h_xVx); hold all; plot(Y(:,1), Y(:,3)); hold off; %FIXME              
                        set(gcf,'CurrentAxes',h_yVy); hold all; plot(h_yVy, Y(:,2), Y(:,4)); hold off; %FIXME
                        set(gcf,'CurrentAxes',h_rV); hold all; plot(h_rV, MySqrt(MySqr(Y(:,1)) + MySqr(Y(:,2))), MySqrt(MySqr(Y(:,3)) + MySqr(Y(:,4)))); hold off; %FIXME
                        plot(h_p_yVy, p_y, p_Vy, 'd');
                    else
                        if ~CheckOutCondition(gamma, a, x_, y_, Vx_, Vy_, f_num)
                            keepedPointNumber = keepedPointNumber + 1;
                            Z(:, :, 1, keepedPointNumber) = [Y(:,1), Y(:,2)];
                        else
                            ecapedPointNumber = ecapedPointNumber + 1;
                            Z(:, :, 2, ecapedPointNumber) = [Y(:,1), Y(:,2)];
                        end
                        
                    end


%                     if length(p_y) > 1
%                         spline_d = 1:length(p_y);
%                         spline_delta = 1:0.1:length(p_y);
%                         splined_p_y = spline (spline_d, p_y, spline_delta);
%                         splined_p_Vy = spline (spline_d, p_Vy, spline_delta);
% 
%                         plot(h_p_yVy, splined_p_y, splined_p_Vy);
%                     end

                %Или группу элементов
                else
                    n = str2double(get(handles.groupPopulationEdit, 'String'));
                    r_group = str2double(get(handles.groupRadiiEdit, 'String'));

                    for i=0:n-1
                       i_theta = (2*i*pi)/n;
                       [x_group, y_group] = pol2cart (i_theta, r_group);

                       x = x_ + x_group;
                       y = y_ + y_group;

                       %Для каждого элемента группы решаем уравнения движения
                       [t, Y] = solveMyODE45(odeFuncHandle, t, [x y Vx_ Vy_], gamma, a);
                       
                       set(gcf,'CurrentAxes',h_main); hold all; plot(Y(:,1), Y(:,2)); hold off;
                       
                    end %for...
                end %else %рассчитываем ли мы одну звезду               
             end %for w
        end %for j
    end %if m ~= 0
    
    if isequal(get(handles.chbPrintCloud, 'Value'), 1)
        Z = permute(Z, ...
                    [4, ... % количество рассчитанных орбит
                     2, ... % (x,y)
                     1, ... % рассчитанные значения для каждой звезды, соответствует t0..tn..t
                     3]);   % типа линии 1 - улетающая, 2 - остающаяся, (по идее можно разделить на 2 независимых массива)
        set(gcf,'CurrentAxes',h_main); hold all; 
        set(h_main, 'NextPlot', 'replaceChildren');
        axis equal;
        [~, ~, tCounts, lineCounts] = size(Z); % tCounts - соответствует количеству рассчитанных точек для i-й звезды
        for k = 1:tCounts
            if lineCounts == 1 % могут быть орбиты только одного типа, поэтому их нужно отображать по-разному
                marker = '';
                if keepedPointNumber > 0
                    marker = '.k';
                elseif ecapedPointNumber > 0
                    marker = '.r';
                end
                plot(Z(:,1,k), Z(:,2,k), marker);
            else
                plot(Z(:,1,k,1), Z(:,2,k,1), '.k', Z(:,1,k,2), Z(:,2,k,2), '.r');
            end            
            if k ~= tCounts
                %timeDiff = Z(1,1,3,k+1) - Z(1,1,3,k);
                timeDiff = t(k+1) - t(k);
                pause(timeDiff);
            end
        end
        hold off;
        
    end
    
    hold off;
    
    result = isFigureUsed;
    
%==========================================================================

function result = MySqr(V)
    r = zeros(1, length(V));
    for i=1:length(V)
        r(i) = V(i)^2;
    end    
    result = r;

%==========================================================================
    
function result = CheckOutCondition(gamma, a0, x, y, Vx, Vy, f_num)
%     g = gamma^2-1;
    g = sqrt(gamma^2-1);
    r = sqrt(x^2 + y^2);
    C1 = y - Vx*gamma/g^2;
    C2 = -x*gamma/g^2 - Vy/g^2;
    C3 = x*gamma/g^3 + Vy*gamma^2/g^3;
    C4 = Vx*gamma/g^2;
    
    X = -C2*gamma;
    Y = C1;
    a = sqrt((g^2/gamma^2)*(C3^2 + C4^2));
    
    R0 = abs(Y) - 2*abs(a);
    R02 = sqrt(X^2 + Y^2);

    A2 = -1/(2*a0^2); A4 = 3/(8*a0^4);
    EA = (2*A2/R0)-4*A4*R0^3;
    V = abs(X)/gamma - gamma*EA/g^3;
    
    %V = abs(X)/gamma - (1/a0^2)*(3*gamma/g^2 + 2/g)*R0*(1+R0^2/a0^2)^(-3/2);
    
    if (R0>0) && (X*Y<0) && (V>=0)
        msg = sprintf('{%d}: Условие выполняется (x = %d, y = %d, Vx = %d, Vy = %d, a0 = %d). \n R0 = %d\n X = %d; Y = %d; \n V = %d', f_num, x, y, Vx, Vy, a, R0, X, Y, V)
        %msgbox(msg);
        msg
        
        result = true;
    else
        msg = sprintf('{%d}: Условие НЕ выполняется (x = %d, y = %d, Vx = %d, Vy = %d, a0 = %d). \n R0 = %d\n X = %d; Y = %d; \n V = %d', f_num, x, y, Vx, Vy, a, R0, X, Y, V);
        %msgbox(msg;;;);
        msg
        
        result = false;
    end

    
%==========================================================================

function result = MySqrt(V)
    r = zeros(1, length(V));
    for i=1:length(V)
        r(i) = sqrt(V(i));
    end    
    result = r;
    
%==========================================================================

function result = GetJacobiValue(handles, x, y, Vx, Vy, a)
    result = 0.5*(Vx^2 + Vy^2) - F(handles, x, y, a) - x^2/2;
        
%==========================================================================

function result = F(~, x, y, a)
%     result = (1+(x^2+y^2)/r0^2)^(-1/2)
    result = (a^2+(x^2+y^2))^(-1/2);

%==========================================================================

% --- Executes on button press in ReMarkButton.
function ReMarkButton_Callback(hObject, ~, handles)
% hObject    handle to ReMarkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    pointsToScan = getappdata(ScanOrbits, 'PointsToScan');
    scannedPoints = getappdata(ScanOrbits, 'ScannedPoints');
    
    set(pointsToScan, 'XData', get(scannedPoints, 'XData'), 'YData', get(scannedPoints, 'YData'))
    set(scannedPoints, 'XData', [], 'YData', [])
    
    setappdata(ScanOrbits, 'PointsToScan', pointsToScan)
    setappdata(ScanOrbits, 'ScannedPoints', scannedPoints)

    set(handles.ClearButton, 'Enable', 'on');
    set(hObject, 'Enable', 'off');
    set(handles.Start, 'Enable', 'on');
    set(handles.CalcCoverage, 'Enable', 'on');

%==========================================================================

% --- Executes during object creation, after setting all properties.
function groupRadiiEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to groupRadiiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%==========================================================================

% --- Executes on button press in doGroupCalcChb.
function doGroupCalcChb_Callback(hObject, ~, handles)
% hObject    handle to doGroupCalcChb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doGroupCalcChb

    if isequal(get(hObject, 'Value'), 0)
        set(handles.groupPopulationEdit, 'Enable', 'Off')
        set(handles.groupRadiiEdit, 'Enable', 'Off')
    else
        set(handles.groupPopulationEdit, 'Enable', 'On')
        set(handles.groupRadiiEdit, 'Enable', 'On')
    end

%==========================================================================

% --- Executes on button press in addpointButton.
function addpointButton_Callback(hObject, ~, handles)
% hObject    handle to addpointButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    AddNewPoint(str2double(get(handles.addIEdit, 'String')), str2double(get(handles.addHEdit, 'String')))
    
    set(handles.Start, 'Enable', 'on');
    set(handles.CalcCoverage, 'Enable', 'on');
    set(handles.ClearButton, 'Enable', 'on');

%==========================================================================

% --- Executes on button press in CalcCoverage.
function CalcCoverage_Callback(hObject, ~, handles)
% hObject    handle to CalcCoverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    currPointer = get(ScanOrbits, 'Pointer');
    set(ScanOrbits, 'Pointer', 'watch')
    
    try
        DoCoverage(handles);
        set(ScanOrbits, 'Pointer', 'Arrow');

        set(ScanOrbits, 'Pointer', currPointer);
        set(hObject, 'Enable', 'off');
    catch ME
        set(ScanOrbits, 'Pointer', currPointer);
        set(hObject, 'Enable', 'off');
        rethrow(ME);
    end

    %==========================================================================

function DoCoverage(handles)

    PointsToScan = getappdata(ScanOrbits, 'PointsToScan');
    ScannedPoints = getappdata(ScanOrbits, 'ScannedPoints');
    
    totalOrbits = 0;
    lostOrbits = 0;
    
    count = length(get(PointsToScan, 'XData'));
    
    p_x = zeros(1, count);
    p_y = zeros(1, count);
    p_z = zeros(1, count);
    
    XDataScanned = zeros(1, count);
    YDataScanned = zeros(1, count);
    
    T_out = [str2double(get(handles.tStartEdit, 'String')) str2double(get(handles.tEndEdit, 'String'))];
    a_out = str2double(get(handles.aEdit, 'String'));
    gamma_out = str2double(get(handles.gammaEdit, 'String'));
       
    pointsCount = size(get(PointsToScan, 'XData'), 2);
    
    groupCalc = get(handles.doGroupCalcChb, 'Value');
    
    Points = getappdata(ScanOrbits, 'PointsToScan');
    XData = get(Points, 'XData');
    YData = get(Points, 'YData');
        
    parfor i = 1:pointsCount        
        T = T_out;
        a = a_out;
        gamma = gamma_out;
        
        h = XData(i);
        E = YData(i);
        
        p_x(i) = h;
        p_y(i) = E;
        
        if isequal(groupCalc, 0)
            [total, lost] = DoOneCoverage(h, E, T, a, gamma);
                        
            if (total > 0)
                p_z(i) = 1 - (lost/total);
            else
                p_z(i) = 0;
            end 
            
            totalOrbits = totalOrbits + total;
            lostOrbits = lostOrbits + lost;
        else
            %DoGroupCoverage...
        end
                
        XDataScanned(i) = XData(i);%i
        YDataScanned(i) = YData(i);%h        
    end
    
    set(PointsToScan, 'XData', XData, 'YData', YData);
    set(ScannedPoints, 'XData', XDataScanned, 'YData', YDataScanned);
        
    setappdata(ScanOrbits, 'PointsToScan', PointsToScan);
    setappdata(ScanOrbits, 'ScannedPoints', ScannedPoints);
    
    %Кнопка "Пометить как нерассчитанные" становится теперь доступна
    set(handles.ReMarkButton, 'Enable', 'on');
    
    hActive = get(0, 'CurrentFigure');
    hNewFigure = figure();
    set(0, 'CurrentFigure', hNewFigure);
    hAxes = axes();
    
    trisurf(delaunay(p_x, p_y), p_x, p_y, p_z);
    DrawLindblad(hAxes, a_out);
    titleText = sprintf('r_0:%0.5g | gamma:%0.5g | lost/total orbits: %d/%d', a_out, gamma_out, lostOrbits, totalOrbits);
    title(titleText);
    
    %switch back active figure
    set(0, 'CurrentFigure', hActive);
    

%==========================================================================
    
function [total, lost] = DoOneCoverage(h, E, T, a, gamma)
    calcsCount = 0;
    lostCount = 0;
    
    r = GetRho(E, h, a);
    
    [~, m] = size(r);
    if m ~= 0
        for j = 1:m
            rho = r(j);
            %w = 0;
             for w = 0:pi*0.1:2*pi                 
                [x_, y_, Vx_, Vy_] = getIC_full(rho, w, E, h, a); %Для случая изолированного скопления - getIC_isolated(rho, w, E, h)
                %Указываем функцию, описывающее уравнения движения
                odeFuncHandle = @myODE_full; % | @myODE_isolatedCluster | @my3dODE_full

                %Решаем уравнения движения
                [~, Y] = solveMyODE45(odeFuncHandle, T, [x_ y_ Vx_ Vy_], gamma, a);

                calcsCount = calcsCount + 1;

                p_x = Y(length(Y(:,1)),1);
                p_y = Y(length(Y(:,1)),2);

                if (p_x*p_x + p_y*p_y) > 10*10
                    lostCount = lostCount + 1;
                end %if
             end %for w
        end %for j
    end %if m ~= 0
    
    total = calcsCount;
    lost = lostCount;
    
%     if (calcsCount > 0)
%         result = 1 - (lostCount / calcsCount);
%     else
%         result = 0;
%     end %if
    
%==========================================================================

function result = DoGroupCoverage(handles, h, E)
    calcsCount = 0;
    lostCount = 0;
    
    t = [str2double(get(handles.tStartEdit, 'String')) str2double(get(handles.tEndEdit, 'String'))];
    a = str2double(get(handles.aEdit, 'String'));
    gamma = str2double(get(handles.gammaEdit, 'String'));
    
    r = GetRho(E, h, a);
    
    [~, m] = size(r);
    if m ~= 0
        for j = 1:m
            rho = r(j);
            %w = 0;
             for w = 0:pi*0.1:0
               [x_, y_, Vx_, Vy_] = getIC_full(rho, w, E, h, a); %Для случая изолированного скопления - getIC_isolated(rho, w, E, h)
               %Указываем функцию, описывающее уравнения движения
               odeFuncHandle = @myODE_full; % | @myODE_isolatedCluster | @my3dODE_full

               n = str2double(get(handles.groupPopulationEdit, 'String'));
               r_group = str2double(get(handles.groupRadiiEdit, 'String'));

                for i=0:n-1
                   i_theta = (2*i*pi)/n;
                   [x_group, y_group] = pol2cart (i_theta, r_group);

                   x = x_ + x_group;
                   y = y_ + y_group;

                   %Для каждого элемента группы решаем уравнения движения
                   [t, Y] = solveMyODE45(odeFuncHandle, t, [x y Vx_ Vy_], gamma, a);

                   calcsCount = calcsCount + 1;

                   p_x = Y(length(Y(:,1)),1);
                   p_y = Y(length(Y(:,1)),2);

                   if (p_x*p_x + p_y*p_y) > 10*10
                       lostCount = lostCount + 1;
                   end %if
                end %for...

             end %for w
        end %for j
    end %if m ~= 0
    
    if (calcsCount > 0)
        result = 1 - (lostCount / calcsCount);
    else
        result = 0;
    end %if
    
%==========================================================================
