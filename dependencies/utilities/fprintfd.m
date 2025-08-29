function fprintfd(varargin)
    % print the inputs if the DEBUG variable is set to 1 in the BASE workspace.
    % 31-01-2024 - R.Waasdorp
    if evalin_if_exist('base','DEBUG',0)
        fprintf(varargin{:});
    end
end
