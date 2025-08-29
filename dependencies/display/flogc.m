function varargout = flogc(varargin)
    % FLOGC performs log compression on the input data
    if nargin == 1
        in = varargin{1};
        out = 20 * log10(abs(in) ./ max(abs(in(:))));
        varargout{1} = out;
    elseif nargin > 1 && nargout > 1
        tmp = cellfun(@(x) 20 * log10(abs(x)), varargin, 'UniformOutput', 0);
        m = max(cellfun(@(x) double(max(x(:))), tmp));
        varargout = cellfun(@(x) x - m, tmp, 'UniformOutput', 0);
    elseif nargin > 1 && nargout == 1
        tmp = cellfun(@(x) 20 * log10(abs(x)), varargin, 'UniformOutput', 0);
        m = max(cellfun(@(x) double(max(x(:))), tmp));
        varargout{1} = cellfun(@(x) x - m, tmp, 'UniformOutput', 0);
    else
        error('dont understand calling convention')
    end
end
