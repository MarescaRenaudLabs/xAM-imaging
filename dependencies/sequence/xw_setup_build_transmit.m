function TX = xw_setup_build_transmit(TX, idx_tx, idx_ray, idx_angle, Trans, P)
    % This function will build the transmit structure for xAM. The transmit will be
    % organized first with the all the X subaperture transmits and then left and then right.
    % It will calculate different subapertures based on the bisector position with the added
    % possibility to receive between elements (non integer bisectors)
    
    idx_bisector = P.ray_positions(idx_ray); % Fetch the bisector

    if idx_bisector - P.aperture_size_max / 2 < 1
        % Left reduced aperture
        Aperture = 2 * ceil(idx_bisector - 1);
    elseif idx_bisector + P.aperture_size_max / 2 > 128
        % Right reduced aperture
        Aperture = 2 * ceil(128 - idx_bisector);
    else
        % Normal aperture
        Aperture = P.aperture_size_max;
    end

    % calculate the angle for the selected Aperture
    if P.use_adaptive_xwave_angle
        P.ray_angle(idx_angle, idx_ray) = acotd((P.max_image_depth * 2) / (Aperture * Trans.spacingMm * 1e-3)); % Changing TX plane wave angle
        P.ray_angle(idx_angle, idx_ray) = min(P.ray_angle(idx_angle, idx_ray), P.xwave_angles(idx_angle));
    else
        P.ray_angle(idx_angle, idx_ray) = P.xwave_angles(idx_angle);
    end

    AngleRad = P.ray_angle(idx_angle, idx_ray) .* pi / 180; % Fetch the angle
    % get the size of the half aperture
    P.ray_half_aperture(idx_angle, idx_ray) = Aperture / 2 + rem(idx_bisector + 0.5, ceil(idx_bisector)); % +rem(32.5+1,32) to find out where is the middle of the aperture!

    % determine max imaging depth for angle and aperture size
    P.ray_max_image_depth(idx_angle, idx_ray) = (Aperture * Trans.spacingMm * 1e-3) / 2 * cotd(P.ray_angle(idx_angle, idx_ray));

    assignin('base', 'P', P);

    % If you want to make sure the end depth is maximized per
    % angle you can change the Angle value
    % AngleRad = min(acot(P.ImagingWindow(2)./P.Wavelength./(Aperture/2)),AngleRad);% The angle has to be the minimum between the
    % fixed angle you give it and the acot(Depth_wanted/Aperture)

    % Build Delay X, L and R based on the bisector number
    [DelayX, DelayL, DelayR, channelsX, channelsL, channelsR] = ...
        xw_setup_compute_delays(idx_bisector, Aperture, AngleRad, Trans, P);
     
    % Calculate apodization if needed
    if (1) % ((Aperture + 5) > P.aperture_size_max) % 5 is arbitrary
        switch P.transmit_apodization
            case 'none'
                ApodList = ones(1, size(channelsL, 2));
            case 'tukey'
                ApodList = tukeywin(size(channelsL, 2), 0.2)';
            case 'hamming'
                ApodList = hamming(size(channelsL, 2))';
            case 'kaiser'
                ApodList = kaiser(size(channelsL, 2))';
        end
    else
        ApodList = ones(1, size(channelsL, 2)); % No apod for small apertures
    end

    % Set the bisector element to 0
    if (size(channelsX, 2) ./ 2) ~= floor(size(channelsX, 2) ./ 2) % here you have a dead element in the middle
        ApodList = [ApodList, 0, ApodList];
    else
        ApodList = [ApodList, ApodList];
    end

    % adaptive TX apod
    % P.minApodVal = 0.70;
    % fapodweight = @(a) (a - P.aperture_size_min) * (1 - P.minApodVal) / (P.aperture_size_max - P.aperture_size_min) + P.minApodVal;
    % ApodList = ApodList * fapodweight(Aperture);

    % 1st transmit is X
    TX(idx_tx).waveform = 1;
    TX(idx_tx).Apod = zeros(1, 128);
    TX(idx_tx).Apod(channelsX) = ApodList; % This will transmit only on 65 elements around ReconChannel
    TX(idx_tx).Delay = DelayX; % X

    % 2nd transmit is /
    TX(idx_tx + P.num_rays * P.num_angles).waveform = 1;
    TX(idx_tx + P.num_rays * P.num_angles).Apod = zeros(1, 128);
    TX(idx_tx + P.num_rays * P.num_angles).Apod(channelsL) = ApodList(1:length(channelsL)); % This will transmit only on 32 elements before ReconChannel
    TX(idx_tx + P.num_rays * P.num_angles).Delay = DelayL; % /

    % 3rd transmit is \
    TX(idx_tx + 2 * P.num_rays * P.num_angles).waveform = 1;
    TX(idx_tx + 2 * P.num_rays * P.num_angles).Apod = zeros(1, 128);
    TX(idx_tx + 2 * P.num_rays * P.num_angles).Apod(channelsR) = ApodList(end:-1:(end - length(channelsR) + 1)); % This will transmit only on 32 elements after ReconChannel
    TX(idx_tx + 2 * P.num_rays * P.num_angles).Delay = DelayR;


end
