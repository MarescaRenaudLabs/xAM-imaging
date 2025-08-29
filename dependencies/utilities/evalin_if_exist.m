function out = evalin_if_exist(WS, expression, default_value)
    % out = evalin_if_exist(WS, expression, default_value)
    %   will return the value of expression if it exists in the provided
    %   workspace, otherwise will return the default value. default value = []
    %   by default.
    %
    % date:    03-04-2023
    % author:  R. Waasdorp (r.waasdorp@tudelft.nl)
    % ==============================================================================
    if ~exist('default_value', 'var') || isempty(default_value)
        default_value = [];
    end
    if evalin(WS, ['exist(''' expression ''',''var'')'])
        out = evalin(WS, expression);
    else
        out = default_value;
    end
end
