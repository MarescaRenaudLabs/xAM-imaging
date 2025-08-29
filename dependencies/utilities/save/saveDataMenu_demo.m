% Short demo of saveDataMenu
%
% date:    30-01-2023
% author:  R. Waasdorp (r.waasdorp@tudelft.nl)
% ==============================================================================

clear
clear saveDataMenu % clear persistence

% get default ui state
SaveDataMenuSettings = saveDataMenu('defaults'); % call with 1 output arg and one input arg to obtain defaults
% and edit defaults
SaveDataMenuSettings.filename = 'test123.mat';
SaveDataMenuSettings.checkboxes.rcvdata = 1;
SaveDataMenuSettings.checkboxes.rf_to_bin = 1;
SaveDataMenuSettings.checkboxes.compression = 0;
SaveDataMenuSettings.custom_variables = {'TEST', 'BLA'};

%% create some empty verasonics fields and data for demo
P = struct();
BFConfig = struct();

Resource = struct();
Trans = struct();
TW = struct();
TX = struct();
Receive = struct();
RcvProfile = struct();
TPC = struct();
TGC = struct();
Media = struct();

IQData = cell(4, 1);
RcvData = {zeros(100, 'int16')};
RcvData = repmat(RcvData, 2, 1);

% set some custom variables
BLA = 13;
TEST = 26;

f = figure(1); clf;
f.Tag = 'plotWindow';
imagesc

%% launch gui to save workspace data
saveDataMenu()
