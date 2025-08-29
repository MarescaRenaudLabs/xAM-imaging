function BFData = xw_beamform_run(RFData)

    fprintfd('Starting beamforming\n')

    P = evalin('base', 'P');
    Receive = evalin('base', 'Receive');
    Trans = evalin('base', 'Trans');
    TW = evalin('base', 'TW');

    % defining parameters structure
    param.p = Trans.spacingMm * 1e-3; % in [m]
    param.lambda = P.scale_wvl2mm * 1e-3; % 1 wvl in [m]
    param.c_us = P.speed_of_sound; % [m/s]
    param.Fs = Receive(1).decimSampleRate * 1e6; % [Hz]
    param.Fc = Receive(1).demodFrequency * 1e6; % [Hz] value was: P.fSI;
    if isfield(P, 'fNumber')
        param.f = P.fNumber; % f number
    else
        param.f = 1.28; % f number
    end

    % set beamforming modes, either run all modes, or select only xAM or
    % xBmode, or subset.
    if isfield(P, 'beamform_image_modes')
        bf_modes = P.beamform_image_modes;
    else
        bf_modes = {'PW', 'xBmode', 'xAM'}; % choose from: or all {'PW', 'xBmode', 'xAM'}
    end
    checkmode = @(mode) ismember(mode, bf_modes); % function returns true if mode is turned ON.

    %% =========================================================================
    % Preprocess RF Data
    % ==========================================================================

    % look at time tag:
    TimeTag = TimeTag2Sec(RFData(1:2, 1));
  
    % define min and max image depth and set z vector
    %z_min =  (P.aperture_size_min - 1) * param.p / 2 * sind(max(P.ray_angle(:))); % [m]
    z_min = P.image_start_depth_mm * 1e-3;
    z_max = P.image_end_depth_mm * 1e-3;
    z_vector = z_min:param.lambda / 4:z_max; % depth axis [m]
    P.x_axis = interp1(1:Trans.numelements, Trans.ElementPos(:, 1), P.ray_positions) * 1e-3;
    P.z_axis = z_vector;

    % Separate RF, and reorder if needed
    Nz_RF = Receive(1).endSample; % Offset for RF data
    ntot = (Nz_RF * P.num_rays * P.num_angles);
    inds = 1:ntot;
    NFrames = size(RFData, 3);
    RFStruct.PW_L = RFData(inds + 1 * ntot, :, :);
    RFStruct.PW_R = RFData(inds + 2 * ntot, :, :);
    RFStruct.XW = RFData(inds + 0 * ntot, :, :);

    if NFrames > 1
        RFStruct.PW_L = reshaper3to2D(RFStruct.PW_L);
        RFStruct.PW_R = reshaper3to2D(RFStruct.PW_R);
        RFStruct.XW = reshaper3to2D(RFStruct.XW);
    end

    % SUM RF for NONLINEAR IMAGING
    RFStruct.XW_AM = RFStruct.XW - RFStruct.PW_R - RFStruct.PW_L;

    % setup LUT, if not in function Workspace
    persistent LUTStruct
    if isempty(LUTStruct)
        LUTStruct = xw_beamform_compute_LUT(bf_modes, z_vector, P, param, Trans, TW);
        LUTStruct = structfun(@single, LUTStruct, 'UniformOutput', false); % cast to single.
    end

    % setup image reconstruction configuration.
    persistent BFC
    if isempty(BFC)
        warning OFF ParseStruct:notExistOrEmpty
        if ~exist('MLBFXWave', 'file') == 3
            error('MLBFXWave:notfound', 'Could not find MLBFXWave, which is the C++ beamformer. Contact r.waasdorp@tudelft.nl')
        end
        BFC.XPiezo = Trans.ElementPos(:, 1) * 1e-3; % in m
        BFC.ZPiezo = 0 * BFC.XPiezo;
        BFC.Origin = [P.z_axis(1); P.x_axis(1)];
        BFC.ScaleZ = P.z_axis(2) - P.z_axis(1);
        BFC.ScaleX = P.x_axis(2) - P.x_axis(1);
        BFC.Nx = numel(P.x_axis);
        BFC.Nz = numel(P.z_axis);
        BFC.Nz_RF = Nz_RF;
        BFC.fNumber = param.f;
        BFC.NTx = P.num_angles;
        BFC.Angles = P.xwave_angles;
        BFC.NAngles = P.num_angles;
        BFC.NRays = P.num_rays;
        BFC.NFrames = NFrames;
        BFC.decimSampleRate = Receive(1).decimSampleRate * 1e6;
        BFC.demodFrequency = Receive(1).demodFrequency * 1e6;
        BFC.speedOfSound = P.speed_of_sound;
        BFC.CorrectionForTW = 0; % already added to LUTs
        BFC.CorrectionForLens = 0;
        BFC.CorrectionForRcvData = 0;
        BFC.connectorNeedOrder = 0;
        % BFC.InterpMode = 'lagrange';
        BFC.InterpMode = 'linear';
        BFC.DemodMode = 'IQ'; % 'IQ' / 'none' / 'hilbert'
        BFC.BFMode = 'XW';
    end

    % Run image reconstruction
    if checkmode('PW')
        PW_L = MLBFXWave(BFC, RFStruct.PW_L, LUTStruct.TX, LUTStruct.RX);
        PW_R = MLBFXWave(BFC, RFStruct.PW_R, LUTStruct.TX, LUTStruct.RX);
        BFData.PW = PW_R + PW_L; % sum left and right plane waves for compound image
    end
    if checkmode('xBmode'); BFData.xBmode = MLBFXWave(BFC, RFStruct.XW, LUTStruct.TX, LUTStruct.RX); end
    if checkmode('xAM'); BFData.xAM = MLBFXWave(BFC, RFStruct.XW_AM, LUTStruct.TX, LUTStruct.RX); end
    BFData.TimeTag = TimeTag; % keep time tag with bf'd data.

    % apply persistence if set.
    persistence = evalin('base', 'UISTATES.persistence');
    if persistence < 1
        assignin('base', 'BFData_raw', BFData); % assigin raw BFData without persistence
        BFData_prev = evalin('base', 'BFData');
        pfun = @(b, b_prev, f, p) p * b.(f) + (1 - p) * b_prev.(f);
        if checkmode('PW'); BFData.PW = pfun(BFData, BFData_prev, 'PW', persistence); end
        if checkmode('xBmode'); BFData.xBmode = pfun(BFData, BFData_prev, 'xBmode', persistence); end
        if checkmode('xAM'); BFData.xAM = pfun(BFData, BFData_prev, 'xAM', persistence); end
    end

    % assign to base workspace
    assignin('base', 'BFData', BFData);
    assignin('base', 'LUTStruct', LUTStruct)
    assignin('base', 'P', P);
    assignin('base', 'BFC', BFC);

end

function BFC = createEmptyBFConfig(P)
    BFC.NFrames = 0;
    BFC.NAngles = 0;
    BFC.NTx = 0;
    BFC.Nx = 0;
    BFC.Nz = 0;
    BFC.IQdemod = 0;
    BFC.IQ_demodFreq = 0;
    BFC.Nz_RF = 0;
    BFC.InterpMode = '';
    BFC.mode = '';
    BFC.NXPiezo = 0;
    BFC.decimSampleRate = 0;

    BFC.ScaleZ = 0;
    BFC.ScaleX = 0;
    BFC.fNumber = 0;
    BFC.ptSourceZ0 = 0;
    BFC.samplesPerWave = 0;
    BFC.demodFrequency = 0;
    BFC.speedOfSound = 0;
    BFC.CorrectionForTW = 0;
    BFC.CorrectionForLens = 0;
    BFC.CorrectionForRcvData = 0;
    BFC.Origin = [0; 0];
    BFC.ZPiezo = zeros(128, 1);
    BFC.XPiezo = zeros(128, 1);
    BFC.Angles = zeros(P.num_rays, 1);
    BFC.addToDelay = zeros(P.num_rays, 1);
    BFC.ptSourceX = zeros(P.num_rays, 1);
    BFC.ptSourceZ = zeros(P.num_rays, 1);
end
