% Description:
%   L22-14v(X) Linear array
%   xWave wide field of view script with constant depth
%   96 ray transmits and 96 receive acquisitions
%   128 receive channels are active for each acquisition of the 96 transmits
%   62.5 MHz sampling for a 15.625 MHz processing center frequency
%
%   Dual mode, and continuous
%
% Last update:
%   03/11/2023 - Rick Waasdorp (r.waasdorp@tudelft.nl) (Beamforming xBmode)
%   25/03/2022 - Baptiste Heiles (Angles)

clear all % need for persistence
format compact
close all

% add paths
addpath('dependencies/display')
addpath('dependencies/reconstruct')
addpath('dependencies/sequence')
addpath('dependencies/utilities')
addpath('dependencies/utilities/save')

DEBUG = 0;

%% =============================================================================
% SETUP ACQUISITION PARAMETERS
% ==============================================================================
P.image_start_depth_mm = 0; % start depth [mm]
P.image_end_depth_mm = 10; % end depth [mm]

P.num_accumulations = 1; % if >1, will accumulate RF data on the device. Keep low to avoid clipping!
P.xwave_angles = [7.5 15.5]; % Ensemble of xWave angles in degrees

P.speed_of_sound = 1480; % default 1480 m/s, CIRS phantom 1540 m/s
P.image_start_voltage = 3; % set to safe number to avoid collapse

P.aperture_size_min = 63; % min number of elements in active aperture (33 for wide FOV)
P.aperture_size_max = 64; % max number of elements in active aperture

P.save_path = 'data'; % default path for data saving

P.transmit_apodization = 'tukey'; % Apodization to avoid edge waves. Options: none/kaiser/hamming/tukey
P.beamform_image_modes = {'PW', 'xBmode', 'xAM'}; % PW = bf(\)+bf(/), xBmode = bf(X), xAM = bf(X-\-/)
P.transmit_frequency = 15.6250; % transducer transmit frequency [MHz]
P.fps = 2; % acquisition frame rate in [Hz]

%% =============================================================================
% TRANSDUCER PARAMETERS
% ==============================================================================

% Define Trans structure array
Trans.name = 'L22-14vX';
Trans.units = 'mm';
Trans = computeTrans(Trans);

%% =============================================================================
% SYSTEM PARAMETERS
% ==============================================================================
% Define system parameters
Resource.Parameters.numTransmit = Trans.numelements; % number of transmit channels
Resource.Parameters.numRcvChannels = Trans.numelements; % number of receive channels
Resource.Parameters.speedOfSound = P.speed_of_sound; % speed of sound [m/s]
Resource.Parameters.verbose = 1; % controls output warnings
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0; %0; simulation onn [1] simulation off [0]
Resource.Parameters.fakeScanhead = 0;
Resource.VDAS.dmaTimeout = 1000 * 4; % 4 seconds
Resource.Parameters.waitForProcessing = 1;

% Compute and set dependent parameters

% some additional dependent parameters
P.scale_wvl2mm = (Resource.Parameters.speedOfSound / Trans.frequency) * 1e-3; % length of 1 wavelength in mm
P.scale_mm2wvl = Trans.frequency / (Resource.Parameters.speedOfSound / 1000); %wavelength per mm
P.num_angles = numel(P.xwave_angles);

P.samples_per_wavelength = 4; % 4 samples per wavelength, hardcoded
P.num_pulses = 3; % number of pulses in xAM pulse sequence, hardcoded
P.timetag_enabled = 1; % debug option, see time tags of frames
P.num_rf_frames = 4; % must be >4 or even number, don't change

% Advanced options
P.use_adaptive_xwave_angle = true; % if true updates the xWave Angles dependent on the  imaging depth
P.use_half_pitch_scanning = true; %

% choose rays for pitch or pitch/2 scanning
step = 1 / (P.use_half_pitch_scanning + 1);
P.half_aperture_size = (P.aperture_size_min - rem(P.aperture_size_min, 2)) / 2; % number of active elements in half transmit aperture
P.ray_positions = ceil(P.aperture_size_min / 2) + 1:step:Resource.Parameters.numTransmit - floor(P.aperture_size_min / 2); % pitch scanning
P.num_rays = numel(P.ray_positions);

% warn user about max angle used
if P.use_adaptive_xwave_angle
    P.max_image_depth = (P.aperture_size_max * Trans.spacingMm * 1e-3) / 2 * cotd(max(P.xwave_angles)) *1e3;
    if P.max_image_depth < P.image_end_depth_mm * 1e-3
        warning('Imaging depth too high for combination of angles. Max depth of intersecting xWave is %.2fmm. Either\n - Decrease imaging depth\n - Set P.use_adaptive_xwave_angle to false\n - Dont trust xAM results below %.2fmm\n - Decrease the max angle', P.max_image_depth * 1e3, P.max_image_depth * 1e3)
    end
end
if P.num_rf_frames < 4; error('num_rf_frames must be larger than 4'); end

% specify Media object
pt1;
Media.function = 'movePoints';


%% =============================================================================
% INPUT VALIDATION
% ==============================================================================

assert(P.image_end_depth_mm >= P.image_start_depth_mm, 'Image end depth must be larger than start depth. Recommended start = 0, end = 10mm')



%% =============================================================================
% RESOURE BUFFERS
% ==============================================================================

% Specify Resources
% Buf len = distance from edge of aperture to max depth
maxAcqLength = ceil(sqrt(P.image_end_depth_mm ^ 2 + ((Trans.numelements - 1) * Trans.spacingMm) ^ 2)); % [mm]
samples_per_mm = (Trans.frequency * P.samples_per_wavelength) / Resource.Parameters.speedOfSound * 1e3;
wavelength_mm = Resource.Parameters.speedOfSound / Trans.frequency * 1e-3; % [mm]
z_delta_mm = wavelength_mm * 0.5; % space between rows of beamform image grid [mm]
maxAcqLengthWvl = round(maxAcqLength * P.scale_mm2wvl);
P.Nz_RF = (ceil(maxAcqLengthWvl * 2 * P.samples_per_wavelength / 128) * 128);

Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.num_rays * P.num_angles * P.num_pulses * P.Nz_RF;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.num_rf_frames;

RcvProfile.LnaZinSel = 31;

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

TPC.hv = P.image_start_voltage;

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
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % reconstruction
Event(n).process = [PN.image_reconstruct_init, PN.image_display]; % process
Event(n).seqControl = SC.matlab; %4;
n = n + 1;

% save acquisition paramters
Event(n).info = 'saveAcqPars';
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % reconstruction
Event(n).process = PN.saveAcqPars; % process
Event(n).seqControl = SC.matlab; %4;
n = n + 1;

Event(n).info = 'enter freeze and stop hardware';
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % reconstruction
Event(n).process = 0; %PN.startFreeze; % process
Event(n).seqControl = [SC.matlab, SC.stop]; %4;
n = n + 1;

% EN.startDop = n;
EN.StartEvent = n;
SeqControl(SC.jump).argument = EN.StartEvent; % set here

for idx_frame = 1:P.num_rf_frames
    for idx_ray = 1:P.num_rays % 1st set of acquisitions
        for idx_angle = 1:P.num_angles
            idx_TX = idx_angle + (idx_ray - 1) * P.num_angles;
            Event(n).info = ['Acquire cross ray line ', num2str(idx_ray), ' angle ', num2str(idx_angle), ' frame ' num2str(idx_frame)];
            Event(n).tx = idx_TX; % use next TX structure.
            Event(n).rcv = idx_TX + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
            Event(n).recon = 0; % no reconstruction.
            Event(n).process = 0; % no processing
            Event(n).seqControl = SC.ttna; % Time between acquisitions
            n = n + 1;

            Event(n).info = ['Acquire right ray line ', num2str(idx_ray), ' angle ', num2str(idx_angle), ' frame ' num2str(idx_frame)];
            Event(n).tx = idx_TX + P.num_rays * P.num_angles; % use next TX structure.
            Event(n).rcv = idx_TX + P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
            Event(n).recon = 0; % no reconstruction.
            Event(n).process = 0; % no processing
            Event(n).seqControl = SC.ttna; % Time between acquisitions
            n = n + 1;

            Event(n).info = ['Acquire left ray line ', num2str(idx_ray), ' angle ', num2str(idx_angle), ' frame ' num2str(idx_frame)];
            Event(n).tx = idx_TX + 2 * P.num_rays * P.num_angles; % use next TX structure.
            Event(n).rcv = idx_TX + 2 * P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
            Event(n).recon = 0; % no reconstruction.
            Event(n).process = 0; % no processing
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
        % and accumulate acqs
        for idx_ray = 1:P.num_rays % 1st set of acquisitions
            for idx_angle = 1:P.num_angles
                idx_TX = idx_angle + (idx_ray - 1) * P.num_angles;
                Event(n).info = ['Acquire cross ray line ', num2str(idx_ray), ' angle ', num2str(idx_angle), ' frame ' num2str(idx_frame), ' accum'];
                Event(n).tx = idx_TX; % use next TX structure.
                Event(n).rcv = idx_TX + nRcvObj_p_accum * P.num_rf_frames + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
                Event(n).recon = 0; % no reconstruction.
                Event(n).process = 0; % no processing
                Event(n).seqControl = SC.ttna; % Time between acquisitions
                n = n + 1;

                Event(n).info = ['Acquire right ray line ', num2str(idx_ray), ' angle ', num2str(idx_angle), ' frame ' num2str(idx_frame), ' accum'];
                Event(n).tx = idx_TX + P.num_rays * P.num_angles; % use next TX structure.
                Event(n).rcv = idx_TX + nRcvObj_p_accum * P.num_rf_frames + P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
                Event(n).recon = 0; % no reconstruction.
                Event(n).process = 0; % no processing
                Event(n).seqControl = SC.ttna; % Time between acquisitions
                n = n + 1;

                Event(n).info = ['Acquire left ray line ', num2str(idx_ray), ' angle ', num2str(idx_angle), ' frame ' num2str(idx_frame), ' accum'];
                Event(n).tx = idx_TX + 2 * P.num_rays * P.num_angles; % use next TX structure.
                Event(n).rcv = idx_TX + nRcvObj_p_accum * P.num_rf_frames + 2 * P.num_rays * P.num_angles + (idx_frame - 1) * nRcvObj_p_accum; % use next receive structure
                Event(n).recon = 0; % no reconstruction.
                Event(n).process = 0; % no processing
                Event(n).seqControl = SC.ttna; % Time between acquisitions
                n = n + 1;
            end
        end
        %         Event(n - 1).seqControl = SC.ttna_frame; % Time between acquisitions

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
        Event(n).tx = 1; % use next TX structure.
        Event(n).rcv = 0; % use next receive structure
        Event(n).recon = 0; % no reconstruction.
        Event(n).process = 0; % no processing
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
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 0; % reconstruction
    Event(n).process = 0; % process
    Event(n).seqControl = nsc; % waitForTransferComplete
    n = n + 1;
    % set waitTransfer to reference last SeqControl
    SeqControl(nsc).command = 'waitForTransferComplete'; % wait for transfer complete, needed for large superframes
    SeqControl(nsc).argument = nsc - 1;
    nsc = nsc + 1;

    Event(n).info = 'Image Reconstruct';
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 0; % no reconstruction
    Event(n).process = PN.image_reconstruct; % external processing function
    Event(n).seqControl = 0; %
    n = n + 1;

    Event(n).info = 'markTransferProcessed';
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 0; % reconstruction
    Event(n).process = 0; % process
    Event(n).seqControl = nsc; % markTransferProcessed
    n = n + 1;
    % set waitTransfer to reference last SeqControl
    SeqControl(nsc).command = 'markTransferProcessed'; % wait for transfer complete, needed for large superframes
    SeqControl(nsc).argument = nsc - 2;
    nsc = nsc + 1;

    %         Event(n).info = 'sync';
    %         Event(n).tx = 0; % no transmit
    %         Event(n).rcv = 0; % no rcv
    %         Event(n).recon = 0; % reconstruction
    %         Event(n).process = 0; % process
    %         Event(n).seqControl = SC.sync; % markTransferProcessed
    %         n=n+1;

    Event(n).info = 'Image Display';
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = 0; % no reconstruction
    Event(n).process = PN.image_display; % external processing function
    Event(n).seqControl = SC.matlab; %
    n = n + 1;

    %     Event(n).info = 'matlab';
    %     Event(n).tx = 0; % no transmit
    %     Event(n).rcv = 0; % no rcv
    %     Event(n).recon = 0; % reconstruction
    %     Event(n).process = 0; % process
    %     Event(n).seqControl = [SC.matlab, SC.sync]; %

end

Event(n).info = 'Jump back to AcquisitionLoop.';
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % no reconstruction
Event(n).process = 0; % no process
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

UISTATES.dr_max = 0; % dynamic range minimum
UI(n_ui).Control = {'UserC8', 'Style', 'VsSlider', 'Label', 'Dyn. Range Max', ...
                        'SliderMinMaxVal', [-18, 0, UISTATES.dr_max], ...
                        'SliderStep', [0.01, 0.1], 'ValueFormat', '%2.1fdB'};
UI(n_ui).Callback = {'xw_change_dr_max(UIValue)'};
n_ui = n_ui + 1;

UISTATES.dr_min = -40; % dynamic range minimum
UI(n_ui).Control = {'UserC7', 'Style', 'VsSlider', 'Label', 'Dyn. Range Min', ...
                        'SliderMinMaxVal', [-60, -22, UISTATES.dr_min], ...
                        'SliderStep', [0.01, 0.1], 'ValueFormat', '%2.1fdB'};
UI(n_ui).Callback = {'xw_change_dr_min(UIValue)'};
n_ui = n_ui + 1;

% button to save data
% VarsToSave = {
%               'P';
%               'Resource';
%               'Trans';
%               'TW';
%               'TX';
%               'Receive';
%               'RcvProfile';
%               'TPC';
%               'TGC';
%               'Media';
%               'UISTATES';
%               'TT_array_movmean';
%               };
% % VarsToSave = [VarsToSave; 'RcvData'];
% VarsToSave = [VarsToSave; 'IQbf_xam'; 'IQbf_xbmode'];

UI(n_ui).Control = {'UserC1', 'Style', 'VsPushButton', 'Label', 'Save Data'};
UI(n_ui).Callback = {'saveDataMenu()'};
n_ui = n_ui + 1;

% configure saveDataMenu
P.matfile_file_name = 'xAM_sequence.mat'; filename = P.matfile_file_name;
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
