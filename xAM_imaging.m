% =========================================================================
% Cross-Amplitude Modulation (xAM) Ultrasound Imaging
% =========================================================================
% 
% This script implements the xAM imaging paradigm for artifact-free detection
% of gas vesicles (GVs). By transmitting cross-propagating plane waves,
% xAM suppresses nonlinear propagation artifacts while preserving nonlinear
% contrast from GVs, as described in:
%
%   Maresca et al., "Nonlinear Ultrasound Imaging of Nanoscale Acoustic Biomolecules"
%   Phys. Rev. X 8, 041002 (2018). DOI: 10.1103/PhysRevX.8.041002
%
% Authors:
%   Rick Waasdorp (r.waasdorp@tudelft.nl)
%   David Maresca (d.maresca@tudelft.nl)
%
% -------------------------------------------------------------------------
% Parameter Reference (struct P)
% -------------------------------------------------------------------------
%
% Basic parameters:
%   P.image_start_depth_mm   : Start depth of imaging region [mm]
%   P.image_end_depth_mm     : End depth of imaging region [mm]
%   P.xwave_angles           : xWave angle [degrees]
%   P.speed_of_sound         : Speed of sound in medium
%                               (1480 m/s in water/agar, 1540 m/s in tissue)
%   P.image_voltage          : Transmit voltage [V]. Set safely to avoid GV collapse
%   P.aperture_size_min      : Minimum number of elements in active aperture
%                               (e.g. 48 for wide FOV)
%   P.aperture_size_max      : Maximum number of elements in active aperture
%   P.transmit_apodization   : Apodization function to reduce edge waves.
%                               Options: 'none', 'kaiser', 'hamming', 'tukey'
%   P.fps                    : Acquisition frame rate [Hz]
%   P.save_path              : Default path for saving data
%
% Advanced options:
%   P.num_accumulations      : Number of RF accumulations. Increases SNR but
%                               lowers FPS and may cause clipping
%   P.use_adaptive_xwave_angle : If true, xWave angles are adapted based on max
%                                 requested imaging depth
%   P.use_half_pitch_scanning  : Enables half-pitch scanning for finer sampling,
%                                 at cost of max frame rate.
%   P.transmit_frequency       : Transducer transmit frequency [MHz]
%
% Dependent parameters (computed from above):
%   P.half_aperture_size     : Half the number of active elements in transmit aperture
%   P.ray_positions          : Transmit ray locations (element indices).
%                               With half-pitch scanning enabled, indices fall
%                               between neighboring elements.
%   P.ray_positions_mm       : Physical positions [mm] of transmit rays in element coordinates
%
% -------------------------------------------------------------------------
% Output Data (struct BFData)
% -------------------------------------------------------------------------
% Data is beamformed on a fixed image grid of P.x_axis and P.z_axis. 
% After beamforming, data is returned in the struct BFData with fields:
%
%   BFData.PW       : B-mode image (complex IQ) reconstructed from plane-wave
%                     data [pixels: depth × lateral]
%
%   BFData.xAM      : xAM image (complex IQ) highlighting GV-specific contrast
%                     [pixels: depth × lateral]
%   
%   BFData.TimeTag  : Time stamp of recorded data (numeric)
%
% =========================================================================

clear all %#ok<CLALL>
format compact
close all

% add paths
addpath('dependencies/display')
addpath('dependencies/reconstruct')
addpath('dependencies/sequence')
addpath('dependencies/utilities')
addpath('dependencies/utilities/save')

%% =============================================================================
% SETUP ACQUISITION PARAMETERS
% ==============================================================================
Trans.name = 'L22-14vX'; % tested transducers are L22-14vX / L11-5v

P.image_start_depth_mm = 0; % start depth [mm]
P.image_end_depth_mm = 10; % end depth [mm]

P.xwave_angles = 17; % xWave angle in degrees

P.speed_of_sound = 1480; % agar/water 1480 m/s, tissue 1540 m/s
P.image_voltage = 2.5; % set to safe number to avoid collapse

P.aperture_size_min = 48; % min number of elements in active aperture (48 for wide FOV)
P.aperture_size_max = 64; % max number of elements in active aperture

P.transmit_apodization = 'tukey'; % Apodization to avoid edge waves. Options: none/kaiser/hamming/tukey
P.fps = 2; % acquisition frame rate in [Hz]

P.save_path = 'data'; % default path for data saving

%% =============================================================================
% TRANSDUCER PARAMETERS
% ==============================================================================

% Define Trans structure array
Trans.units = 'mm';
Trans = computeTrans(Trans);

%% =============================================================================
% SYSTEM PARAMETERS
% ==============================================================================
% Define system parameters
Resource.Parameters.numTransmit = 128; % number of transmit channels
Resource.Parameters.numRcvChannels = 128; % number of receive channels
Resource.Parameters.speedOfSound = P.speed_of_sound; % speed of sound [m/s]
Resource.Parameters.verbose = 1; % controls output warnings
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0; % 0; simulation onn [1] simulation off [0]
Resource.Parameters.fakeScanhead = 0;
Resource.VDAS.dmaTimeout = 1000 * 4; % 4 seconds
Resource.Parameters.waitForProcessing = 1;

% Compute and set dependent parameters

% Advanced options
P.num_accumulations = 1; % if >1, will accumulate RF data on the device. Keep low to avoid clipping! Lowers max FPS!
P.use_adaptive_xwave_angle = true; % if true updates the xWave Angles dependent on the  imaging depth
P.use_half_pitch_scanning = true;
P.transmit_frequency = Trans.frequency; % transducer transmit frequency [MHz]

P.num_pulses = 3; % number of pulses in xAM pulse sequence, hardcoded
P.timetag_enabled = 1; % debug option, see time tags of frames
P.num_rf_frames = 4; % must be >4 or even number, don't change
P.beamform_image_modes = {'PW', 'xAM'}; % PW = bf(\)+bf(/), xBmode = bf(X), xAM = bf(X-\-/)

% choose rays for pitch or pitch/2 scanning
step = 1 / (P.use_half_pitch_scanning + 1);
P.half_aperture_size = (P.aperture_size_min - rem(P.aperture_size_min, 2)) / 2; % number of active elements in half transmit aperture
P.ray_positions = ceil(P.aperture_size_min / 2) + 1:step:Resource.Parameters.numTransmit - floor(P.aperture_size_min / 2); % pitch scanning
P.ray_positions_mm = interp1(1:Trans.numelements, Trans.ElementPos(:, 1), P.ray_positions);
P.num_rays = numel(P.ray_positions);

% warn user about max angle used
if P.use_adaptive_xwave_angle
    P.max_image_depth = (P.aperture_size_max * Trans.spacingMm * 1e-3) / 2 * cotd(max(P.xwave_angles));
    if P.max_image_depth < P.image_end_depth_mm * 1e-3
        warning('Imaging depth too high for combination of angles. Max depth of intersecting xWave is %.2fmm. Either\n - Decrease imaging depth\n - Set P.use_adaptive_xwave_angle to false\n - Dont trust xAM results below %.2fmm\n - Decrease the max angle', P.max_image_depth * 1e3, P.max_image_depth * 1e3)
    end
end
assert(P.num_rf_frames >= 4, 'num_rf_frames must be larger than 4');

% some additional dependent parameters
P.scale_wvl2mm = (Resource.Parameters.speedOfSound / Trans.frequency) * 1e-3; % length of 1 wavelength in mm
P.scale_mm2wvl = Trans.frequency / (Resource.Parameters.speedOfSound / 1000); %wavelength per mm
P.num_angles = numel(P.xwave_angles);

% specify Media object
pt1;
Media.function = 'movePoints';

%% =============================================================================
% INPUT VALIDATION
% ==============================================================================
assert(P.image_end_depth_mm >= P.image_start_depth_mm, 'Image end depth must be larger than start depth. Recommended start = 0, end = 10mm')
assert(P.aperture_size_min < P.aperture_size_max, 'Min aperture size must be smaller than max aperture size. Recommended min = 33, max = 64')

%% =============================================================================
% RESOURE BUFFERS
% ==============================================================================

% determine demod frequency, and compute Nz_RF
TxFreq = 250 ./ (2 .* ((6:197).'));
TxFreqValid = TxFreq(2 .* TxFreq < 3.5 & TxFreq > 1.5); % depends on Trans.Bandwidth [1.5 3.5] for P4-1
demodFreqsValid = TxFreq(rem(250, TxFreq * 4) == 0);
[~, idx] = min(abs(Trans.frequency - demodFreqsValid));
demodFrequency = demodFreqsValid(idx);
% fprintf('\tTX    FREQ: %.3f MHz\n', Trans.frequency)
% fprintf('\tDEMOD FREQ: %.3f MHz\n', demodFrequency)

maxAcqLength = ceil(sqrt(P.image_end_depth_mm ^ 2 + ((Trans.numelements - 1) * Trans.spacingMm) ^ 2)); % [mm]
maxAcqLengthWvl = maxAcqLength * P.scale_mm2wvl;
P.Nz_RF = (ceil(maxAcqLengthWvl * 2 * 4 * demodFrequency / Trans.frequency / 128) * 128); %2*4*enddepth*factor (beam forming factor)

% allocate resources
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.num_rays * P.num_angles * P.num_pulses * P.Nz_RF;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.num_rf_frames;

% RcvProfile.LnaZinSel = 31; % turn on high Z-state, increase sensitivity

%% Timing parameters
P.num_tx_per_frame = P.num_accumulations * P.num_angles * P.num_pulses * P.num_rays;
P.ttna = round(maxAcqLength * 5 * 1e-3 / Resource.Parameters.speedOfSound * 1e6); % in microseconds!
P.time_per_frame = P.ttna * P.num_tx_per_frame;
P.max_framerate = 1e6 / P.time_per_frame;

%% =============================================================================
% TRANSMIT WAVEFORM
% ==============================================================================
% Specify Transmit Waveform (TW) structure array
TW(1).type = 'parametric';
TW(1).Parameters = [P.transmit_frequency, 0.67, 2, 1];

%% =============================================================================
% TRANSMIT AND RECEIVE STRUCTURES
% ==============================================================================

% Specify TX structure array
% Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0, 0.0, 0.0], ...
    'Steer', [0.0, 0.0], ...
    'Apod', zeros(1, Trans.numelements), ...
    'focus', 0, ...
    'Delay', zeros(1, Trans.numelements)), ...
    1, P.num_pulses * P.num_angles * P.num_rays);

TPC.hv = P.image_voltage;

% Specify Receive structure arrays.
P.use_accumulation = P.num_accumulations > 1;
nRcvObj_p_accum = P.num_pulses * P.num_rays * P.num_angles;
nRcvObj_p_frame = nRcvObj_p_accum * (P.use_accumulation + 1);
Receive = repmat(struct('Apod', zeros(1, Trans.numelements), ...
    'startDepth', 0, ... % P.image_start_depth_mm * P.scale_mm2wvl, ...
    'endDepth', maxAcqLengthWvl, ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ... % newly acquired data is written to the storage memory, replacing any previous data
    'callMediaFunc', 0), ...
    1, nRcvObj_p_frame * P.num_rf_frames);

% Set event specific Receive and TX attributes.
% Compute transmit delays
P.ray_angle = zeros(P.num_angles, P.num_rays); % xWave angles per ray
for idx_ray = 1:P.num_rays
    for idx_angle = 1:P.num_angles
        idx_tx = idx_angle + (idx_ray - 1) * P.num_angles;
        TX = xw_setup_build_transmit(TX, idx_tx, idx_ray, idx_angle, Trans, P);

        for idx_frame = 1:P.num_rf_frames
            idx = idx_angle + (idx_ray - 1) * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum;

            if (idx_ray == 1 && idx_angle == 1)
                Receive(idx).callMediaFunc = 1; % 1 for all 3 pulses
            end

            Receive(idx).Apod(:) = 1.0; % activate all elements
            Receive(idx).mode = 0;
            Receive(idx).framenum = idx_frame;
            Receive(idx).acqNum = idx_angle + (idx_ray - 1) * P.num_angles; % for X pulse

            Receive(idx + P.num_angles * P.num_rays).Apod = Receive(idx).Apod;
            Receive(idx + P.num_angles * P.num_rays).mode = 0;
            Receive(idx + P.num_angles * P.num_rays).framenum = idx_frame;
            Receive(idx + P.num_angles * P.num_rays).acqNum = idx_angle + (idx_ray - 1) * P.num_angles + P.num_angles * P.num_rays; % for right active aperture

            Receive(idx + 2 * P.num_angles * P.num_rays).Apod = Receive(idx).Apod;
            Receive(idx + 2 * P.num_angles * P.num_rays).mode = 0;
            Receive(idx + 2 * P.num_angles * P.num_rays).framenum = idx_frame;
            Receive(idx + 2 * P.num_angles * P.num_rays).acqNum = idx_angle + (idx_ray - 1) * P.num_angles + P.num_angles * 2 * P.num_rays; % for left active aperture

            if P.use_accumulation
                idx = idx_angle + (idx_ray - 1) * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum + nRcvObj_p_accum * P.num_rf_frames;

                Receive(idx).mode = 1; % activate all elements
                Receive(idx).Apod(:) = 1.0; % activate all elements
                Receive(idx).framenum = idx_frame;
                Receive(idx).acqNum = idx_angle + (idx_ray - 1) * P.num_angles; % for X pulse

                Receive(idx + P.num_angles * P.num_rays).mode = 1;
                Receive(idx + P.num_angles * P.num_rays).Apod = Receive(idx).Apod;
                Receive(idx + P.num_angles * P.num_rays).framenum = idx_frame;
                Receive(idx + P.num_angles * P.num_rays).acqNum = idx_angle + (idx_ray - 1) * P.num_angles + P.num_angles * P.num_rays; % for right active aperture

                Receive(idx + 2 * P.num_angles * P.num_rays).mode = 1;
                Receive(idx + 2 * P.num_angles * P.num_rays).Apod = Receive(idx).Apod;
                Receive(idx + 2 * P.num_angles * P.num_rays).framenum = idx_frame;
                Receive(idx + 2 * P.num_angles * P.num_rays).acqNum = idx_angle + (idx_ray - 1) * P.num_angles + P.num_angles * 2 * P.num_rays; % for left active aperture
            end
        end
    end
end

%% =============================================================================
% TIME GAIN CONTROL
% ==============================================================================
% Specify TGC Waveform structure (time-gain-compensation)
TGC(1).CntrlPts = 800 * ones(1, 8); % 1023 is max gain
TGC(1).rangeMax = P.image_end_depth_mm * P.scale_mm2wvl; % in wavelengths
TGC(1).Waveform = computeTGCWaveform(TGC(1));

%% =============================================================================
% PROCESS FUNCTIONS
% ==============================================================================

% specify external processing event
np = 1;
PN.image_reconstruct = np;
Process(np).classname = 'External';
Process(np).method = 'xw_beamform_run';
Process(np).Parameters = {'srcbuffer', 'receive', ...
                              'srcbufnum', 1, ...
                              'srcframenum', -1, ...
                              'dstbuffer', 'none'};
np = np + 1;

PN.image_reconstruct_init = np;
Process(np).classname = 'External';
Process(np).method = 'xw_beamform_run';
Process(np).Parameters = {'srcbuffer', 'receive', ...
                              'srcbufnum', 1, ...
                              'srcframenum', 1, ...
                              'dstbuffer', 'none'};
np = np + 1;

PN.image_display = np;
Process(np).classname = 'External';
Process(np).method = 'xw_plot';
Process(np).Parameters = {'srcbuffer', 'none', ...
                              'dstbuffer', 'none'};
np = np + 1;

PN.saveAcqPars = np;
Process(np).classname = 'External'; % process structure for 1st Doppler ensemble
Process(np).method = 'xw_sequence_save_acquisition_parameters';
Process(np).Parameters = {'srcbuffer', 'none', ... % name of buffer to process.
                              'dstbuffer', 'none'};
np = np + 1;

%% =============================================================================
% SEQUENCE CONTROL
% ==============================================================================

% Specify SeqControl structure arrays.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = P.ttna;
SC.ttna = 1;
SeqControl(2).command = 'returnToMatlab';
SC.matlab = 2;
SeqControl(3).command = 'jump';
SeqControl(3).argument = []; % set later
SC.jump = 3;
SeqControl(4).command = 'stop'; % time between frames
SC.stop = 4;
SeqControl(5).command = 'loopCnt';
SeqControl(5).argument = (P.num_accumulations - 1);
SC.loopCnt = 5;
SeqControl(6).command = 'timeToNextAcq';
P.ttna_frame = 1e6 / P.fps - P.ttna * (P.num_tx_per_frame);
assert(P.ttna_frame > P.ttna, 'frame rate too high! Max is %.2f Hz', P.max_framerate)
SeqControl(6).argument = P.ttna_frame; %  4e6; % 4 seconds max! in microseconds!
SC.ttna_frame = 6;

SeqControl(7).command = 'sync';
SeqControl(7).argument = 1000 * 4;
SC.sync = 7;

nsc = length(SeqControl) + 1;

%% =============================================================================
% EVENT STRUCTURES
% ==============================================================================

% Specify Event structure arrays.
n = 1;

Event(n).info = 'initialize beamformer';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = [PN.image_reconstruct_init, PN.image_display];
Event(n).seqControl = SC.matlab;
n = n + 1;

% save acquisition paramters
Event(n).info = 'saveAcqPars';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = PN.saveAcqPars;
Event(n).seqControl = SC.matlab;
n = n + 1;

Event(n).info = 'enter freeze and stop hardware';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = [SC.matlab, SC.stop];
n = n + 1;

% EN.startDop = n;
EN.StartEvent = n;
SeqControl(SC.jump).argument = EN.StartEvent; % set here

for idx_frame = 1:P.num_rf_frames
    for idx_ray = 1:P.num_rays % 1st set of acquisitions
        for idx_angle = 1:P.num_angles
            idx_TX = idx_angle + (idx_ray - 1) * P.num_angles;
            Event(n).info = sprintf('Acquire X R%03i A%02i F%02i', idx_ray, idx_angle, idx_frame);
            Event(n).tx = idx_TX;
            Event(n).rcv = idx_TX + (idx_frame - 1) * nRcvObj_p_accum;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = SC.ttna; % Time between acquisitions
            n = n + 1;

            Event(n).info = sprintf('Acquire R R%03i A%02i F%02i', idx_ray, idx_angle, idx_frame);
            Event(n).tx = idx_TX + P.num_rays * P.num_angles;
            Event(n).rcv = idx_TX + P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = SC.ttna; % Time between acquisitions
            n = n + 1;

            Event(n).info = sprintf('Acquire L R%03i A%02i F%02i', idx_ray, idx_angle, idx_frame);
            Event(n).tx = idx_TX + 2 * P.num_rays * P.num_angles;
            Event(n).rcv = idx_TX + 2 * P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = SC.ttna; % Time between acquisitions
            n = n + 1;
        end
    end
    if ~P.use_accumulation
        Event(n - 1).seqControl = SC.ttna_frame; % Time between acquisitions
    end

    if P.use_accumulation
        Event(n).info = 'Set loop count for number of accumulates.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = SC.loopCnt;
        n = n + 1;

        Event(n).info = 'Jump to end of accumulate events for loop count test.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'jump'; % Argument set below.
        nsc = nsc + 1;
        n = n + 1;

        nstart = n;

        % and accumulation acqs
        for idx_ray = 1:P.num_rays % 1st set of acquisitions
            for idx_angle = 1:P.num_angles
                idx_TX = idx_angle + (idx_ray - 1) * P.num_angles;
                Event(n).info = sprintf('Acquire X R%03i A%02i F%02i accum', idx_ray, idx_angle, idx_frame);
                Event(n).tx = idx_TX;
                Event(n).rcv = idx_TX + nRcvObj_p_accum * P.num_rf_frames + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = SC.ttna; % Time between acquisitions
                n = n + 1;

                Event(n).info = sprintf('Acquire R R%03i A%02i F%02i accum', idx_ray, idx_angle, idx_frame);
                Event(n).tx = idx_TX + P.num_rays * P.num_angles;
                Event(n).rcv = idx_TX + nRcvObj_p_accum * P.num_rf_frames + P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = SC.ttna; % Time between acquisitions
                n = n + 1;

                Event(n).info = sprintf('Acquire L R%03i A%02i F%02i accum', idx_ray, idx_angle, idx_frame);
                Event(n).tx = idx_TX + 2 * P.num_rays * P.num_angles;
                Event(n).rcv = idx_TX + nRcvObj_p_accum * P.num_rf_frames + 2 * P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = SC.ttna; % Time between acquisitions
                n = n + 1;
            end
        end

        SeqControl(nsc - 1).argument = n;
        Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'loopTst';
        SeqControl(nsc).argument = nstart;
        nsc = nsc + 1;
        n = n + 1;

        Event(n).info = 'dummy tx';
        Event(n).tx = 1;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = SC.ttna_frame; % Time between acquisitions
        n = n + 1;
    end

    Event(n).info = 'Transfer frame to host buffer';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc + 1;
    n = n + 1;

    Event(n).info = 'waitForTransferComplete';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc; % waitForTransferComplete
    n = n + 1;
    % set waitTransfer to reference last SeqControl
    SeqControl(nsc).command = 'waitForTransferComplete'; % wait for transfer complete, needed for large superframes
    SeqControl(nsc).argument = nsc - 1;
    nsc = nsc + 1;

    Event(n).info = 'Image Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = PN.image_reconstruct;
    Event(n).seqControl = 0; %
    n = n + 1;

    Event(n).info = 'markTransferProcessed';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc; % markTransferProcessed
    n = n + 1;
    % set waitTransfer to reference last SeqControl
    SeqControl(nsc).command = 'markTransferProcessed'; % wait for transfer complete, needed for large superframes
    SeqControl(nsc).argument = nsc - 2;
    nsc = nsc + 1;

    Event(n).info = 'Image Display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = PN.image_display;
    Event(n).seqControl = SC.matlab; %
    n = n + 1;

end

Event(n).info = 'Jump back to AcquisitionLoop.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = [SC.matlab, SC.jump]; % jump back to Event 1
n = n + 1;

%% =============================================================================
% User specified UI Control Elements
% =============================================================================
if P.timetag_enabled == 1
    UI(1).Statement = 'xw_setup_timetag()';
end

n_ui = 1;
% button to restart sequence
UI(n_ui).Control = {'UserB1', 'Style', 'VsPushButton', 'Label', 'Start Seq.'};
UI(n_ui).Callback = {'xw_sequence_start(hObject)'};
n_ui = n_ui + 1;

UISTATES = struct;

% Change persistence
P.image_persistence = 1;
UISTATES.persistence = P.image_persistence;
UI(n_ui).Control = {'UserB6', 'Style', 'VsSlider', 'Label', 'Persistence', ...
                        'SliderMinMaxVal', [0.1, 1, P.image_persistence], ...
                        'SliderStep', [0.01, 0.1], 'ValueFormat', '%1.2f'};
UI(n_ui).Callback = {'xw_change_persistence(UIValue)'};
n_ui = n_ui + 1;

UISTATES.dr_bmode = 40; % dynamic range B-mode
UI(n_ui).Control = {'UserC8', 'Style', 'VsSlider', 'Label', 'DR B-Mode', ...
                        'SliderMinMaxVal', [0, 60, UISTATES.dr_bmode], ...
                        'SliderStep', [0.01, 0.1], 'ValueFormat', '%2.1fdB'};
UI(n_ui).Callback = {'xw_change_dr(''bmode'', UIValue)'};
n_ui = n_ui + 1;

UISTATES.dr_xam = 20; % dynamic range minimum
UI(n_ui).Control = {'UserC7', 'Style', 'VsSlider', 'Label', 'DR xAM', ...
                        'SliderMinMaxVal', [0, 30, UISTATES.dr_xam], ...
                        'SliderStep', [0.01, 0.1], 'ValueFormat', '%2.1fdB'};
UI(n_ui).Callback = {'xw_change_dr(''xam'', UIValue)'};
n_ui = n_ui + 1;

UI(n_ui).Control = {'UserC1', 'Style', 'VsPushButton', 'Label', 'Save Data'};
UI(n_ui).Callback = {'saveDataMenu()'};
n_ui = n_ui + 1;

% configure saveDataMenu
P.matfile_file_name = 'xAM_imaging.mat'; filename = P.matfile_file_name;
SaveDataMenuSettings = saveDataMenu('defaults'); % call with 1 output arg and one input arg to obtain defaults
SaveDataMenuSettings.filename = P.matfile_file_name;
SaveDataMenuSettings.custom_figures = {2};
SaveDataMenuSettings.folder = P.save_path;
SaveDataMenuSettings.custom_variables = {'UISTATES', 'BFData', 'BFData_raw'};
BFData_raw = []; BFData = [];

% save sequence definition for VSX
if ~isfolder('MatFiles'); mkdir('MatFiles'); end
addpath('MatFiles')
save(['MatFiles/', P.matfile_file_name]);

% Automatically Launch Verasonics
VSX
