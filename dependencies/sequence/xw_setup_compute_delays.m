function [DelayX, DelayL, DelayR, channelsX, channelsL, channelsR] = xw_setup_compute_delays(BisectorIdx, Aperture, AngleRad, Trans, P)
    % Function to calculate delays for xBmode/xAM
    
    % Initialize values
    DelayX = zeros(1, 128); DelayL = DelayX; DelayR = DelayX;
    theta = AngleRad;
    c0 = P.speed_of_sound;

    if ~strcmp(Trans.units, 'mm')
        warning('off', 'backtrace');
        warning('BuildDelayVector not implemented in wavelength, please set Trans.units to mm');
        warning('on', 'backtrace');
        return
    end

    % Calculate Delay once in wavelength! (first in seconds then scaled to match Verasonics wavelength)
    if floor(BisectorIdx) == BisectorIdx
        Delay = Trans.ElementPos(1:Aperture/2, 1) .* 1E-3 .* sin(theta) / c0; % Delay in seconds of the center frequency (1/f)
        Delay = Delay * P.transmit_frequency * 1e6;
        Delay= [Delay; 0;flip(Delay)];
        Delay = Delay - min(Delay(:));
        
        try
            Delay(Aperture / 2 + 1) = 0; % Put the delay at SubAperture/2+1 to 0 (not necessary)
        catch
            keyboard
        end

        % Then put this Delay at the right spot
        % Construct channel indices
        channelsX = BisectorIdx + [- (Aperture / 2):Aperture / 2];
        channelsL = channelsX(1:Aperture / 2);
        channelsR = channelsX((Aperture / 2 + 2):end);

        % Fill in DelayX, DelayL, DelayR
        DelayX(channelsX) = Delay; % Delay X are the delays for the X-transmission
        DelayL(channelsL) = Delay(1:Aperture / 2); % Delay L are the delays for the /-transmission
        DelayR(channelsR) = Delay(1 + [Aperture / 2 + 1:Aperture]); % Delay R are the delays for the \-transmission
    else
        % In that case there will be no dead element in the middle, we transmit in between elements
        Delay = Trans.ElementPos(1:Aperture/2, 1) .* 1E-3 .* sin(theta) / c0; % Delay in seconds of the center frequency (1/f)
        Delay = Delay * P.transmit_frequency * 1e6;
        Delay = [Delay; flip(Delay)];        
        Delay = Delay -  min(Delay(:));
        
        % Then put this Delay at the right spot
        % Construct channel indices and remove the last one to be sure
        % you don't have a dead element
        ListElements = [-(Aperture / 2):Aperture / 2]; ListElements(end) = [];
        channelsX = BisectorIdx + 0.5 + ListElements;
        channelsL = channelsX(1:Aperture / 2);
        channelsR = channelsX((Aperture / 2 + 1):end);

        % Fill in DelayX, DelayL, DelayR
        DelayX(channelsX) = Delay; % Delay X are the delays for the X-transmission
        DelayL(channelsL) = Delay(1:Aperture / 2); % Delay L are the delays for the /-transmission
        DelayR(channelsR) = Delay([Aperture / 2 + 1:Aperture]); % Delay R are the delays for the \-transmission
    end

end
