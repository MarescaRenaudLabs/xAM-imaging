function varargout = saveDataMenu(varargin)
    % Menu to save data from Verasonics GUI. Select which parameters to save,
    % and customize parameters with the SaveDataMenuSettings config.
    %
    % See saveDataMenu_demo.m
    %
    % date:    30-01-2023
    % author:  R. Waasdorp (r.waasdorp@tudelft.nl)
    % ==========================================================================

    % variables
    persistent previous_state fig
    filename_field = [];
    folder_field = [];
    checkboxes = struct('vsx_pars', [], 'iqdata', [], 'rcvdata', [], 'figfig', [], ...
        'figpng', [], 'compression', []);

    % if called with one input, one output, and 'command', run cli_gateway
    if nargin == 1 && nargout == 1
        varargout{1} = cli_gateway(varargin{1});
        return;
    end

    % ==========================================================================
    % MAIN ROUTINE
    % ==========================================================================

    % check if VSX is freezed!
    if ~vsx_check_freeze()
        return
    end

    % Create a figure for the GUI
    if isempty(fig) || ~isgraphics(fig)
        fig = uifigure('Name', 'Save Data...', 'HandleVisibility', 'on', 'Tag', 'SaveDataMenu');
        create_components();
        fig.Position(4) = 360;
    else
        figure(get_handle())
    end

    % set defaults or previous state or custom defaults
    if ~isempty(previous_state) % check for previous call
        state = previous_state;
    elseif check_if_exist('base', 'SaveDataMenuSettings') % check if UIState exists in base workspace
        state = evalin_if_exist('base', 'SaveDataMenuSettings');
    else % just use default state
        state = get_default_settings();
    end
    restore_state(state);
    previous_state = get_state(state);

    % ==========================================================================
    % LOCAL METHODS
    % ==========================================================================

    function create_components()

        col_width = 80;
        row_height = 25;

        % Set the layout of the GUI
        layout = uigridlayout(fig, [2, 3]);
        layout.ColumnWidth = {col_width, '1x', col_width, col_width};
        layout.RowHeight = repelem({row_height}, 9); % {40, 40,'1x'};

        % Create a text label for the filename field
        filename_label = uilabel(layout, 'Text', 'Filename:');
        set_position(filename_label, 1, 1);
        filename_field = uieditfield(layout, 'text');
        set_position(filename_field, 1, [2 3]);

        datetime_button = uibutton(layout, 'Text', 'Add Time');
        datetime_button.ButtonPushedFcn = @add_datetime_filename;
        set_position(datetime_button, 1, 4);

        % Create a text label for the folder field
        folder_label = uilabel(layout, 'Text', 'Folder:');
        set_position(folder_label, 2, 1);
        folder_field = uieditfield(layout, 'Text');
        set_position(folder_field, 2, [2 3]);

        browse_button = uibutton(layout, 'Text', 'Browse');
        browse_button.ButtonPushedFcn = @browse_folder_callback;
        set_position(browse_button, 2, 4);

        % optional variables
        vars_label = uilabel(layout, 'Text', 'Variables:');
        set_position(vars_label, 3, 1);

        checkboxes.vsx_pars = uicheckbox(layout, 'Text', 'Verasonics structs');
        set_position(checkboxes.vsx_pars, 3, 2);

        checkboxes.iqdata = uicheckbox(layout, 'Text', 'IQData');
        set_position(checkboxes.iqdata, 4, 2);

        checkboxes.rcvdata = uicheckbox(layout, 'Text', 'RcvData');
        set_position(checkboxes.rcvdata, 5, 2);

        vars_label = uilabel(layout, 'Text', 'Figures:');
        set_position(vars_label, 6, 1);

        % optional figure saving
        checkboxes.figfig = uicheckbox(layout, 'Text', '.fig');
        set_position(checkboxes.figfig, 6, 2);

        checkboxes.figpng = uicheckbox(layout, 'Text', '.png');
        set_position(checkboxes.figpng, 7, 2);

        % create checkboxes for save options
        saveoptions_label = uilabel(layout, 'Text', 'Save options:');
        set_position(saveoptions_label, 8, 1);

        checkboxes.compression = uicheckbox(layout, 'Text', 'Compression');
        set_position(checkboxes.compression, 8, 2);

        % save and cancel
        cancel_button = uibutton(layout, 'Text', 'Reset');
        cancel_button.ButtonPushedFcn = @reset_callback;
        set_position(cancel_button, 7, 4);

        cancel_button = uibutton(layout, 'Text', 'Cancel');
        cancel_button.ButtonPushedFcn = @cancel_callback;
        set_position(cancel_button, 8, 4);

        save_button = uibutton(layout, 'Text', 'Save');
        save_button.ButtonPushedFcn = @save_data_callback;
        set_position(save_button, 9, 4);
    end

    function out = cli_gateway(cmd)
        switch cmd
            case 'defaults'
                out = get_default_settings();
            case 'variables'
                reset_callback();
                out = get_vars_to_save();
        end
    end
    % --------------------------------------------------------------------------
    % settings/state management
    % --------------------------------------------------------------------------
    function state = get_default_settings()
        state.filename = '';
        state.folder = '';
        state.checkboxes = struct( ...
            'vsx_pars', 1, ...
            'iqdata', 1, ...
            'rcvdata', 0, ...
            'figfig', 1, ...
            'figpng', 1, ...
            'compression', 1);
        state.custom_variables = {};
        state.keep_open = true;
    end

    function state = get_state(state)
        state.filename = filename_field.Value;
        state.folder = folder_field.Value;
        state.checkboxes = struct();
        for fn = fieldnames(checkboxes).'
            f = fn{1};
            state.checkboxes.(f) = checkboxes.(f).Value;
        end
    end

    function restore_state(state)
        filename_field.Value = state.filename;
        folder_field.Value = state.folder;
        for fn = fieldnames(checkboxes).'
            f = fn{1};
            checkboxes.(f).Value = state.checkboxes.(f);
        end
    end

    % --------------------------------------------------------------------------
    % button callbacks
    % --------------------------------------------------------------------------
    function h = get_handle()
        h = findobj('tag', 'SaveDataMenu');
    end

    function browse_folder_callback(~, ~)
        folder = uigetdir();
        folder_field.Value = folder;
        figure(get_handle()) % to pop back to front
    end

    function cancel_callback(~, ~)
        close_menu();
    end

    function reset_callback(~, ~)
        % restore_state(get_default_settings());
        state = evalin_if_exist('base', 'SaveDataMenuSettings');
        restore_state(state);
    end

    function set_position(widget, row, col)
        widget.Layout.Row = row;
        widget.Layout.Column = col;
    end

    function close_menu()
        close(get_handle())
    end

    function error_dialog(title, msg)
        fig_error = errordlg(msg, title);
        set(fig_error, 'WindowStyle', 'modal')
        uiwait(fig_error);
        return
    end

    function status_dialog(title, msg)
        fig_status = msgbox(msg, title);
        set(fig_status, 'WindowStyle', 'modal')
        uiwait(fig_status);
        return
    end

    function add_datetime_filename(~, ~)
        str = datestr(now, 'dd-mm-yy_HHMMSS');
        filename = filename_field.Value;
        [~, filename, ~] = fileparts(filename);
        try %#ok
            % will fail if its not a datestr
            datenum(filename(end - 14:end), 'dd-MM-yy_HHmmSS');
            filename(end - 15:end) = []; % remove '_date'
        end
        filename = [filename '_' str '.mat'];
        filename_field.Value = filename;
    end

    function save_data_callback(~, ~)
        folder = folder_field.Value;
        if ~isfolder(folder) && ~isempty(folder)
            mkdir(folder)
        end
        filename = filename_field.Value;
        if isempty(filename)
            error_dialog('Error', 'Filename cannot be empty!')
            return
        end

        if ~endsWith(filename, '.mat')
            filename = [filename '.mat'];
            filename_field.Value = filename;
        end
        out_path = fullfile(folder, filename);

        % check what to save
        if all(structfun(@(f) f.Value == 0, checkboxes))
            error_dialog('Error', 'Nothing to save. Select something to save.');
            return
        end

        % construct string to save stuff
        vars_to_save = get_vars_to_save();

        % save options
        save_options = {};
        if ~state.checkboxes.compression
            save_options = [save_options; '-nocompression'];
        end

        % call save in base workspace
        call_save_mat(out_path, vars_to_save, save_options);


        % save figures if requested
        if state.checkboxes.figfig || state.checkboxes.figpng
            save_plot(out_path);
        end

        previous_state = get_state(state); % save menu state for possible new image to save

        if state.keep_open
            status_dialog('Saved!', 'Data saved!');
            return
        else
            close_menu(); % close gui
        end

    end

    % --------------------------------------------------------------------------
    % save routines
    % --------------------------------------------------------------------------

    function call_save_mat(save_path, vars_to_save, save_options)
        varStr = char(join(compose('''%s'',', string(vars_to_save))));
        if exist('save_options', 'var') && ~isempty(save_options)
            optionStr = char(join(compose('''%s'',', string(save_options))));
            varStr = [varStr, optionStr];
        end
        saveStr = ['save(''' save_path ''', ' varStr(1:end - 1) ');'];
        disp(saveStr)
        evalin('base', saveStr);
    end

    function vars_to_save = get_vars_to_save()
        state = get_state(state);
        vsx_check_freeze();
        vsx_get_buffers(); % get verasonics buffers

        vars_to_save_vsx = {'Resource'; 'Trans'; 'TW'; 'TX'; 'Receive'; 'RcvProfile'; 'TPC'; 'TGC'; 'Media'};
        vars_to_save = {};

        if state.checkboxes.vsx_pars
            vars_to_save = [vars_to_save; vars_to_save_vsx];
            paramw = evalin('base', 'whos(''BFConfig'', ''P'',''Pm'');');
            vars_to_save = [vars_to_save; {paramw.name}.'];
        end

        if state.checkboxes.iqdata
            n = numel(vars_to_save);
            % search for IQ in base workspace
            iqw = evalin('base', 'whos(''IQ*'');');
            vars_to_save = [vars_to_save; {iqw.name}.'];

            % get verasonics buffers
            if check_if_exist('base', 'IData') % comes in pairs
                vars_to_save = [vars_to_save; {'IData'; 'QData'}];
            end
            if n == numel(vars_to_save)
                error_dialog('Error', 'No IQ data found!');
            end
        end

        if state.checkboxes.rcvdata
            n = numel(vars_to_save);
            if check_if_exist('base', 'RcvData') % comes in pairs
                vars_to_save = [vars_to_save; 'RcvData'];
            end
            rcv = evalin('base', 'whos(''RF*'');');
            vars_to_save = [vars_to_save; {rcv.name}.'];
            if n == numel(vars_to_save)
                error_dialog('Error', 'No RF data found!');
            end
        end

        % add potential custom variables
        if ~isempty(state.custom_variables)
            if ~iscell(state.custom_variables)
                state.custom_variables = {state.custom_variables};
            end
            if numel(state.custom_variables) >= 1 && ~isempty(state.custom_variables{1})
                vars_to_save = [vars_to_save; state.custom_variables(:)];
            end
        end

         % add potential custom variables
        if isfield(state, 'custom_figures') && ~isempty(state.custom_figures)
            if ~iscell(state.custom_figures)
                state.custom_figures= {state.custom_figures};
            end
        end
    end

    % --------------------------------------------------------------------------
    % save plots
    % --------------------------------------------------------------------------

    function save_plot(save_path)
        save_path_no_ext = save_path(1:end - 4);

        % check if custom figures need to be saved
        if isfield(state, 'custom_figures') && ~isempty(state.custom_figures)
            for k = 1:numel(state.custom_figures)
                figN = state.custom_figures{k};
                if ~ishandle(figN); continue; end
                figH = figure(figN);
                save_path = sprintf('%s_fig%i',save_path_no_ext, k);
                if state.checkboxes.figfig
                    savefig(figH, [save_path '.fig']);
                end
                if state.checkboxes.figpng
                    exportgraphics(figH, [save_path '.png'], 'Resolution', 300)
                end

            end
        end
        
        % get handle to plotwindow
        plotWindow = findobj('tag', 'plotWindow');
        if isempty(plotWindow) || ~isgraphics(plotWindow)
            return
        end
        if state.checkboxes.figfig
            savefig(plotWindow, [save_path_no_ext '.fig']);
        end
        if state.checkboxes.figpng
            % export_fig(plotWindow, [save_path_no_ext '.png'], '-r600');
            exportgraphics(plotWindow, [save_path_no_ext '.png'], 'Resolution', 600)
        end

    end

    % --------------------------------------------------------------------------
    % interface with verasonics
    % --------------------------------------------------------------------------
    function isfrozen = vsx_check_freeze()
        if ~isempty(findobj('tag', 'UI')) % running VSX
            if evalin('base', 'freeze') == 0 % no action if not in freeze
                msgbox('Please freeze VSX');
                isfrozen = 0;
                return
            end
        end
        isfrozen = 1;
    end

    function vsx_get_buffers()
        if ~isempty(findobj('tag', 'UI')) % running VSX
            Control.Command = 'copyBuffers';
            assignin('base', 'Control', Control)
            evalin('base', 'runAcq(Control)'); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
        end
    end
end
