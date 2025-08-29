function tt_S = TimeTag2Sec(tt_RF)
    % Convert a verasonics time tag to seconds.
    % Input can be multiple time tags, need to be stored as [2xn] array
    timeTags = double(reshape(tt_RF, 2, [], 1));
    timeTags(1, timeTags(1, :) < 0) = timeTags(1, timeTags(1, :) < 0) + 65536; % fix int roll over
    tt_S = sum(timeTags .* [1; 65536]) / 4e4; % and convert to seconds
end
