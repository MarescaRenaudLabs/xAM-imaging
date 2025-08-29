function LUTStruct = xw_beamform_compute_LUT(bf_modes, z_vector, P, param, Trans, TW)
    % LUTStruct = bfXwaveComputeLuts(modes,z_vector, P, param, Trans, TW)
    % Computes LUT for beamforming for x-wave imaging, but also for left and right
    % plane waves.
    %
    %
    % Input arguments:  (see call signature in bfXwaveWideFOVconstDepth_v3)
    %
    %   modes:      specify the requested modes for which to generate the lookup
    %               table, cell array. If all modes are needed,
    %                   modes = {'PW', 'xBmode', 'xAM'}
    %
    % Output arguments:
    %   LUTStruct:  struct that contains all the LUTs, if all modes requested,
    %               you'll end up with the fields, 'Tx_PW_L', 'Tx_PW_R',
    %               'Tx_XW' and 'Rcv'.
    %               The Tx_* luts dimensions are [nz, nrays, nangles]
    %               The Rcv (Receive) lut is [nz, nrays, nelements]
    %
    %
    % date:    02-11-2023
    % author:  R. Waasdorp (r.waasdorp@tudelft.nl)
    % ==============================================================================

    % size_LUT_T = [length(z_vector), P.num_rays, P.num_angles]; % size LUTS
    % size_LUT_R = [length(z_vector), P.num_rays, Trans.numelements]; % size LUTS

    ListRaysPos_x = interp1(1:Trans.numelements, Trans.ElementPos(:, 1), P.ray_positions) * 1e-3; % in m
    ElementPos_x = Trans.ElementPos(:, 1) * 1e-3; % in m

    % Experiment with lens correction
    CorrectionForTW = TW(1).peak / (Trans.frequency * 1e6);
    CorrectionForLens = Trans.lensCorrection * 1e-3 ./ param.c_us;
    toffset = CorrectionForTW + 2 * CorrectionForLens;

    % if ismember('PW_L', modes)
    %     LUT_T = zeros(length(z_vector), P.num_rays, P.num_angles);
    %     for i_Ray = 1:P.num_rays
    %         for i_z = 1:length(z_vector)
    %             for i_Angle = 1:P.num_angles
    %                 % forward delay
    %                 Da = P.ray_half_aperture(i_Angle, i_Ray) * param.p; % size of active aperture [m]
    %                 tTx = (z_vector(i_z) * cosd(P.ray_angle(i_Angle, i_Ray)) + Da * sind(P.ray_angle(i_Angle, i_Ray))) / param.c_us;
    %                 LUT_T(i_z, i_Ray, i_Angle) = tTx;
    %             end
    %         end
    %     end
    %     LUTStruct.Tx_PW_L = LUT_T;
    % end

    % if ismember('PW_R', modes)
    %     LUT_T = zeros(length(z_vector), P.num_rays, P.num_angles);
    %
    %     for i_Ray = 1:P.num_rays
    %         for i_z = 1:length(z_vector)
    %             for i_Angle = 1:P.num_angles
    %                 % forward delay
    %                 Da = P.ray_half_aperture(i_Angle, i_Ray) * param.p; % size of active aperture [m]
    %                 tTx = (z_vector(i_z) * cosd(P.ray_angle(i_Angle, i_Ray)) + Da * sind(P.ray_angle(i_Angle, i_Ray))) / param.c_us;
    %                 LUT_T(i_z, i_Ray, i_Angle) = tTx;
    %             end
    %         end
    %     end
    %     LUTStruct.Tx_PW_R = LUT_T;
    % end

    % if ismember('xAM', modes) || ismember % Not needed anymore.
    LUT_T = zeros(length(z_vector), P.num_rays, P.num_angles);
    for i_Ray = 1:P.num_rays
        for i_z = 1:length(z_vector)
            for i_Angle = 1:P.num_angles
                % forward delay
                Da = P.ray_half_aperture(i_Angle, i_Ray) * param.p; % size of active aperture [m]
                tTx = (Da * sind(P.ray_angle(i_Angle, i_Ray)) + z_vector(i_z) ...
                    * cosd(P.ray_angle(i_Angle, i_Ray))) / param.c_us;
                LUT_T(i_z, i_Ray, i_Angle) = tTx;
            end
        end
    end
    LUTStruct.TX = LUT_T;
    %end

    % always compute Rcv LUT - this contains the offsets for TW and lens
    LUT_R = zeros(length(z_vector), P.num_rays, Trans.numelements);

    for i_Ray = 1:P.num_rays
        xrecon = ListRaysPos_x(i_Ray);
        for i_elem = 1:Trans.numelements
            xpiezo = ElementPos_x(i_elem);
            DistPiezo = abs(xrecon - xpiezo);
            for i_z = 1:length(z_vector)
                % compute aperture based on depth and f#
                aperture = z_vector(i_z) / param.f;
                % if position satisfies f#, determine delay
                if (DistPiezo < (aperture / 2))
                    % computing tof from point A to receiver transducer element
                    tRx = sqrt(z_vector(i_z) .^ 2 + DistPiezo .^ 2) / param.c_us;
                    LUT_R(i_z, i_Ray, i_elem) = tRx + toffset;
                end
            end
        end
    end
    LUTStruct.RX = LUT_R;

end
