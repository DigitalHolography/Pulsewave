% Vesselness_Filter_Tool
% Script to implement GUI for Vessleness_Filter_Tool. Allows user to load
% an image, select filter parameters, filter the image with a Frangi
% vesselness filter, and save the resulting image.
%
% MPHYGB24: Exercise Sheet #6 (Tasks 2 & 3)
% AUTHOR:   Edward James, Medical Imaging CDT MRes , UCL, 2017

function varargout = vesselness_filter_tool(varargin)
    % VESSELNESS_FILTER_TOOL MATLAB code for vesselness_filter_tool.fig
    %      VESSELNESS_FILTER_TOOL, by itself, creates a new VESSELNESS_FILTER_TOOL or raises the existing
    %      singleton*.
    %
    %      H = VESSELNESS_FILTER_TOOL returns the handle to a new VESSELNESS_FILTER_TOOL or the handle to
    %      the existing singleton*.
    %
    %      VESSELNESS_FILTER_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in VESSELNESS_FILTER_TOOL.M with the given input arguments.
    %
    %      VESSELNESS_FILTER_TOOL('Property','Value',...) creates a new VESSELNESS_FILTER_TOOL or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before vesselness_filter_tool_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to vesselness_filter_tool_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help vesselness_filter_tool

    % Last Modified by GUIDE v2.5 08-Dec-2017 09:49:13

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name', mfilename, ...
        'gui_Singleton', gui_Singleton, ...
        'gui_OpeningFcn', @vesselness_filter_tool_OpeningFcn, ...
        'gui_OutputFcn', @vesselness_filter_tool_OutputFcn, ...
        'gui_LayoutFcn', [], ...
        'gui_Callback', []);

    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end

    % End initialization code - DO NOT EDIT

    % --- Executes just before vesselness_filter_tool is made visible.
    function vesselness_filter_tool_OpeningFcn(hObject, eventdata, handles, varargin)
        % This function has no output args, see OutputFcn.
        % hObject    handle to figure
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % varargin   command line arguments to vesselness_filter_tool (see VARARGIN)

        %To initialise the GUI sigma and beta values:
        start_sigma = 2.7;
        start_beta = 0.5;
        %To put these values into the edit boxes and slider
        set(handles.Sigma_Slider, 'Value', start_sigma);
        set(handles.Sigma_EditBox, 'String', num2str(start_sigma));
        set(handles.Beta_Editbox, 'String', num2str(start_beta));

        % Choose default command line output for vesselness_filter_tool
        handles.output = hObject;

        % Update handles structure
        guidata(hObject, handles);

        % UIWAIT makes vesselness_filter_tool wait for user response (see UIRESUME)
        % uiwait(handles.figure1);

        % --- Outputs from this function are returned to the command line.
        function varargout = vesselness_filter_tool_OutputFcn(hObject, eventdata, handles)
            % varargout  cell array for returning output args (see VARARGOUT);
            % hObject    handle to figure
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)

            % Get default command line output from handles structure
            varargout{1} = handles.output;

            function Beta_Editbox_Callback(hObject, eventdata, handles)
                % hObject    handle to Beta_Editbox (see GCBO)
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)
                beta = str2double(get(hObject, 'String')); % get beta value
                guidata(hObject, handles); % update handles

                % --- Executes during object creation, after setting all properties.
                function Beta_Editbox_CreateFcn(hObject, eventdata, handles)
                    % hObject    handle to Beta_Editbox (see GCBO)
                    % eventdata  reserved - to be defined in a future version of MATLAB
                    % handles    empty - handles not created until after all CreateFcns called
                    if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
                        set(hObject, 'BackgroundColor', 'white');
                    end

                    % --- Executes on button press in Open_Pushbutton.
                    function Open_Pushbutton_Callback(hObject, eventdata, handles)
                        % hObject    handle to Open_Pushbutton (see GCBO)
                        % eventdata  reserved - to be defined in a future version of MATLAB
                        % handles    structure with handles and user data (see GUIDATA)

                        % Load image.
                        % File format drop-down box will show .png, .gif, .jpeg, or all file types.
                        [basefilename, folder] = uigetfile({'*.png'; '*.gif'; '*.jpeg'; '*.*'}, ...
                        'Select an image to show');

                        if basefilename == 0
                            % User has clicked the cancel button or selected a file that does not
                            % exist, for example
                            % Display an error message
                            errordlg('File not found', 'Upload Error', 'modal');
                            return; % GUI behaves as before (i.e. it will not crash)
                        end

                        fullfilename = fullfile(folder, basefilename);

                        % To catch an error if file selected is not a recognised image type
                        try
                            handles.original = imread(fullfilename);
                        catch
                            errordlg('File selected is not a recognised image type', ...
                                'Upload Error', 'modal');
                            return; % GUI behaves as before (i.e. it will not crash)
                        end

                        % Show image in first axes
                        % axes(handles.axes1);
                        imshow(handles.original);

                        % Update handles structure
                        guidata(hObject, handles);

                        % --- Executes on button press in Save_Pushbutton.
                        function Save_Pushbutton_Callback(hObject, eventdata, handles)
                            % hObject    handle to Save_Pushbutton (see GCBO)
                            % eventdata  reserved - to be defined in a future version of MATLAB
                            % handles    structure with handles and user data (see GUIDATA)

                            %To save filtered image as one of the below types
                            img = handles.filtered;
                            [filename, pathname] = uiputfile({'*.png'; '*.gif'; '*.jpeg'; '*.*'}, ...
                            'Save filtered image as');

                            if filename == 0
                                % User has clicked the cancel button, for example.
                                % Display an error message
                                errordlg('File not Saved', 'Save Error', 'modal');
                                return; % GUI behaves as before (i.e. it will not crash)
                            end

                            file = fullfile(pathname, filename);

                            % To catch an error if user attempts to save file as file type that is not
                            % recognised as an image
                            try
                                imwrite(img, file);
                            catch
                                errordlg('File not Saved - File format selected is not a recognised image type', ...
                                    'Save Error', 'modal');
                                return; % GUI behaves as before (i.e. it will not crash)
                            end

                            % --- Executes on button press in Filter_Pushbutton.
                            function Filter_Pushbutton_Callback(hObject, eventdata, handles)
                                % hObject    handle to Filter_Pushbutton (see GCBO)
                                % eventdata  reserved - to be defined in a future version of MATLAB
                                % handles    structure with handles and user data (see GUIDATA)

                                %Obtain values of beta and sigma
                                sigma = get(handles.Sigma_EditBox, 'String');
                                sigma = str2double(sigma);
                                beta = get(handles.Beta_Editbox, 'String');
                                beta = str2double(beta);

                                %Load original image
                                img = handles.original;
                                img = im2double(img); % convert image to double precision
                                img = imcomplement(img); % invert the image intensities
                                img = vesselness_filter(img, sigma, beta); % filter the image
                                axes(handles.axes2);
                                imshow(img); %show image in second axes

                                %Store and update handles structure.
                                handles.filtered = img;
                                guidata(hObject, handles);

                                % --- Executes on slider movement.
                                function Sigma_Slider_Callback(hObject, eventdata, handles)
                                    % hObject    handle to Sigma_Slider (see GCBO)
                                    % eventdata  reserved - to be defined in a future version of MATLAB
                                    % handles    structure with handles and user data (see GUIDATA)
                                    sigma = get(hObject, 'Value'); % returns position of slider
                                    set(handles.Sigma_EditBox, 'String', num2str(sigma)) % update edit box
                                    guidata(hObject, handles) % update handles

                                    % --- Executes during object creation, after setting all properties.
                                    function Sigma_Slider_CreateFcn(hObject, eventdata, handles)
                                        % hObject    handle to Sigma_Slider (see GCBO)
                                        % eventdata  reserved - to be defined in a future version of MATLAB
                                        % handles    empty - handles not created until after all CreateFcns called

                                        % Hint: slider controls usually have a light gray background.
                                        if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
                                            set(hObject, 'BackgroundColor', [.9 .9 .9]);
                                        end

                                        function Sigma_EditBox_Callback(hObject, eventdata, handles)
                                            % hObject    handle to Sigma_EditBox (see GCBO)
                                            % eventdata  reserved - to be defined in a future version of MATLAB
                                            % handles    structure with handles and user data (see GUIDATA)
                                            sigma = str2double(get(hObject, 'String')); % get value from box
                                            set(handles.Sigma_Slider, 'Value', sigma); % update slider
                                            guidata(hObject, handles); % update handles

                                            % --- Executes during object creation, after setting all properties.
                                            function Sigma_EditBox_CreateFcn(hObject, eventdata, handles)
                                                % hObject    handle to Sigma_EditBox (see GCBO)
                                                % eventdata  reserved - to be defined in a future version of MATLAB
                                                % handles    empty - handles not created until after all CreateFcns called
                                                if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
                                                    set(hObject, 'BackgroundColor', 'white');
                                                end

                                                % --- Executes during object creation, after setting all properties.
                                                function figure1_CreateFcn(hObject, eventdata, handles)
                                                    % hObject    handle to figure1 (see GCBO)
                                                    % eventdata  reserved - to be defined in a future version of MATLAB
                                                    % handles    empty - handles not created until after all CreateFcns called

                                                    % --- Executes during object creation, after setting all properties.
                                                    function axes2_CreateFcn(hObject, eventdata, handles)
                                                        % hObject    handle to axes2 (see GCBO)
                                                        % eventdata  reserved - to be defined in a future version of MATLAB
                                                        % handles    empty - handles not created until after all CreateFcns called

                                                        % --- Executes during object creation, after setting all properties.
                                                        function axes1_CreateFcn(hObject, eventdata, handles)
                                                            % hObject    handle to axes1 (see GCBO)
                                                            % eventdata  reserved - to be defined in a future version of MATLAB
                                                            % handles    empty - handles not created until after all CreateFcns called

                                                            % --- Executes during object deletion, before destroying properties.
                                                            function axes2_DeleteFcn(hObject, eventdata, handles)
                                                                % hObject    handle to axes2 (see GCBO)
                                                                % eventdata  reserved - to be defined in a future version of MATLAB
                                                                % handles    structure with handles and user data (see GUIDATA)
