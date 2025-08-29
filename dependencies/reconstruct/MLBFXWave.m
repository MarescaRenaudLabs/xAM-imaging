% ==============================================================================
% BFIQ = MLBFXWave(BFConfig, RF, LUT_T, LUT_R)          to use LUT Beamforming
% ------------------------------------------------------------------------------
% GPU Beamforming implementation that uses LUTs. Requires input struct BFConfig
% and RF data. Supports 2D/3D plane wave, and 2D X-wave. Uses Demodulator if
% needed. If supported, can also compute traveltime LUTS, otherwise, supply
% externally computed LUTs.
%
%   Input arguments
%       - BFConfig      defines the beamform parameters.
%       - RF            RF Data. RF must be int16.
%       - LUT_T         Traveltime lookuptable forward delay
%       - LUT_R         Traveltime lookuptable return delay
%   Output arguments
%       - BFIQ          Beamformed image. Complex floats.
%
% ------------------------------------------------------------------------------
%
% BFConfig field information:
%                 InterpMode    Interpolation mode, 'linear/nearest/lagrange's
%                     BFMode    Beamform mode, 'PW2D'/'PW3D'/'XW'. If 2D mode is
%                                 used, Y fields are optional in BFConfig struct.
%                  DemodMode    Demodulation mode, 'none/IQ/hilbert'. If DemodMode
%                                 is set to 'none', the data is real, the
%                                 interpolation algortihm will assume that the
%                                 data is sampled at 4 points per wavelength,
%                                 to produce an intermediate analyitical signal.
%                        NTx    Number of transmissions per frame  (redundant probably)
%                    NAngles    Number of angles per frame         (mostly used instead of NTx)
%                    NFrames    Number of frames to reconstuct
%            decimSampleRate    The sampling frequency in MHz
%             demodFrequency    The IQ demodulation frequency in MHz (Fs/4 most likely)
%            CorrectionForTW    Correction for time to peak (s)
%          CorrectionForLens    Correction for lens travel time (s)
%       CorrectionForRcvData    Correction for start depth (not correcty implemented yet)
%             CorrectionFor*    Note: CorrectionFor* fields can be left 0 when
%                                 using LUT, since already added to T or R LUT.
%               speedOfSound    Speed of sound (m/s)
%                    fNumber    fNumber. When using LUTs, f number mask must be
%                                 applied to LUT_R. Set LUT_R to 0 (or actually
%                                 < 1e-15s) to not use that element in BF.
%                      Nz_RF    Number of samples per channel per acquisition.
%                                 Calculate using Receive.endSample-Receive.startSample
%                     Origin    Array with the Origin of the reconstruction grid.
%                                 For 2D must be in format [z,x]. For 3D [x,y,z].
%                 Nx, Ny, Nz    Number of grid points to reconstruct in x,y,z
%     ScaleX, ScaleY, ScaleZ    Scale of the reconstruction grid in x,y,z
%     XPiezo, YPiezo, ZPiezo    Array of Piezo element positions in x,y,z
%      ptSourceX, ...Y, ...Z    Array of Point source location
%                 addToDelay    Delay to add due to traveltime in virtual medium behind
%                                 transducer. Point source sorcery.
%
% date:    06-11-2023
% author:  R. Waasdorp (r.waasdorp@tudelft.nl)
% ==============================================================================
