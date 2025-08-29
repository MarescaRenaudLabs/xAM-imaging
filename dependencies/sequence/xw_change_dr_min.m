function xw_change_dr_min(UIValue)
    UISTATES = evalin('base', 'UISTATES');
    UISTATES.dr_min = UIValue;
    assignin('base', 'UISTATES', UISTATES);
    xw_plot(); % update plot
end
