function xw_change_persistence(UIValue)
    UISTATES = evalin('base', 'UISTATES');
    UISTATES.persistence = UIValue;
    assignin('base', 'UISTATES', UISTATES);
end
