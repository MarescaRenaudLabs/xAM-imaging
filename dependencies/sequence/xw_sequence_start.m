function xw_sequence_start(hObject)
    EN = evalin('base', 'EN');
    fprintf('Starting acquisition\n')
    eventNum = EN.StartEvent;

    % create control struct
    Control = evalin('base', 'Control');
    nc = numel(Control);
    % if ~isempty(Control(nc).Command)
    %     nc = nc + 1;
    % end
    Control(nc + 1).Command = 'set&Run';
    Control(nc + 1).Parameters = {'Parameters', 1, 'startEvent', eventNum}; % Event number that corresponds to start sequence
    assignin('base', 'Control', Control);

    % reset StatusDisplay
    %StatusDisplay(true);

    % unfreeze
    gui = hObject.Parent; % handle to gui
    freezeBut = findobj(gui.Children, 'String', 'Freeze'); % handle to button
    freezeBut.Value = 0;
    assignin('base', 'freeze', 0);
end
