function xw_change_dr(src, value)
    UISTATES = evalin('base', 'UISTATES');
    % update UI State
    switch src
        case 'bmode'
            UISTATES.dr_bmode = value;
        case 'xam'
            UISTATES.dr_xam = value;
    end
    % make persistent
    assignin('base', 'UISTATES', UISTATES);
    xw_plot(); % update plot
end
