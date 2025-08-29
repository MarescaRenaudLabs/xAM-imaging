function xw_change_dr_max(UIValue)
    UISTATES = evalin('base', 'UISTATES');
    UISTATES.dr_max = UIValue;
    assignin('base', 'UISTATES', UISTATES);
    xw_plot(); % update plot
end
