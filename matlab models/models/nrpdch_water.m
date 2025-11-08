%% NR PDSCH Throughput 
% This reference simulation shows how to measure the physical downlink
% shared channel (PDSCH) throughput of a 5G New Radio (NR) link, as defined
% by the 3GPP NR standard. The example implements the PDSCH and downlink
% shared channel (DL-SCH). The transmitter model includes PDSCH
% demodulation reference signals (DM-RS) and PDSCH phase tracking reference
% signals (PT-RS). The example supports both clustered delay line (CDL) and
% tapped delay line (TDL) propagation channels. You can perform perfect or
% practical synchronization and channel estimation. To reduce the total
% simulation time, you can execute the SNR points in the SNR loop in
% parallel by using the Parallel Computing Toolbox(TM).

% Copyright 2017-2023 The MathWorks, Inc.

%% Introduction
% This example measures the PDSCH throughput of a 5G link, as defined by
% the 3GPP NR standard [ <#13 1> ], [ <#13 2> ], [ <#13 3> ], [ <#13 4> ].
% 
% The example models these 5G NR features:
% 
% * DL-SCH transport channel coding
% * Multiple codewords, dependent on the number of layers
% * PDSCH, PDSCH DM-RS, and PDSCH PT-RS generation
% * Variable subcarrier spacing and frame numerologies (2^n * 15 kHz)
% * Normal and extended cyclic prefix
% * TDL and CDL propagation channel models
% 
% Other features of the simulation are:
% 
% * PDSCH subband precoding using SVD
% * CP-OFDM modulation
% * Slot wise and non slot wise PDSCH and DM-RS mapping
% * Perfect or practical synchronization and channel estimation
% * HARQ operation with 16 processes
% * The example uses a single bandwidth part across the whole carrier
% 
% The figure shows the implemented processing chain. For clarity, the 
% DM-RS and PT-RS generation are omitted.
% 
% <<../PDSCHLinkExampleProcessingChain.png>>
%
% For a more detailed explanation of the steps implemented in this example,
% see <docid:5g_gs#mw_b41487f2-4e53-412e-83ff-b77bd3984819 Model 5G NR
% Communication Links> and
% <docid:5g_gs#mw_5441cb0f-86ed-4038-8735-d0a50ea68f60 DL-SCH and PDSCH
% Transmit and Receive Processing Chain>.
%
% This example supports both wideband and subband precoding. The precoding
% matrix is determined using SVD by averaging the channel estimate across
% all PDSCH PRBs in the allocation (wideband case) or in the subband.
% 
% To reduce the total simulation time, you can use the Parallel Computing
% Toolbox to execute the SNR points of the SNR loop in parallel.
load('csi_1.mat');
% antenna_idx = 1; user_idx = 1; subcarrier_idx = 10;
% channel_amplitude = abs(squeeze(Channel_result_f(antenna_idx, user_idx, :, subcarrier_idx)));
% 
% figure;
% plot(channel_amplitude);
% title('Channel Magnitude vs. Time');
% xlabel('Time Index');
% ylabel('Magnitude');
% 
% correlation = corrcoef(channel_amplitude(1:end-1), channel_amplitude(2:end));
% disp('Correlation between adjacent frames:');
% disp(correlation(1,2));

%% 

% 参数设置
selected_antennas = 1:4; % 选择 4 根天线
selected_user = 1;       % 选择用户（假设用户索引为 1）
num_subcarriers = 624;   % 目标子载波数
num_symbols = 14;        % OFDM 符号数
num_traces = size(Channel_result_f, 3); % 100 个信道 trace

% 初始化信道存储 (624 x 14 x 4)
h = zeros(num_subcarriers, num_symbols, numel(selected_antennas));

% 循环处理 100 个信道 trace 中的第一个
trace_idx = 50; % 假设使用第一个 trace（可按需调整）

% 提取对应天线、用户和 trace 的信道
selected_channel = Channel_result_f(selected_antennas, selected_user, trace_idx, :); % 4x1x1024

% 插值到 624 子载波
for ant = 1:numel(selected_antennas) % 遍历天线
    for sym = 1:num_symbols % 遍历符号
        % 提取原始 1024 子载波信道
        original_h = squeeze(selected_channel(ant, :, :)); % 维度: 1x1024
        % 插值到 624 子载波
        interpolated_h = interp1(1:1024, original_h, linspace(1, 1024, num_subcarriers), 'linear');
        % 存储到目标信道格式
        h(:, sym, ant) = interpolated_h; % 存储为 624x14x4
    end
end


%绘制Argos信道响应
antenna_idx = 1; % 选择绘制的天线索引
symbol_idx = 1;  % 选择绘制的 OFDM 符号索引

% 提取信道响应
h_magnitude = abs(h(:, symbol_idx, antenna_idx)); % 幅度
h_phase = angle(h(:, symbol_idx, antenna_idx));   % 相位

% 子载波索引
subcarrier_idx = 1:size(h, 1);

%绘制幅度响应
figure;
plot(subcarrier_idx, h_magnitude, 'LineWidth', 2);
grid on;
xlabel('Subcarrier Index');
ylabel('Channel Magnitude');
title(sprintf('Channel Magnitude Response (Antenna %d, Symbol %d)', antenna_idx, symbol_idx));

figure;
hold on;
for antenna_idx = 1:size(h, 3) % 遍历所有天线
    h_magnitude = abs(h(:, symbol_idx, antenna_idx));
    plot(subcarrier_idx, h_magnitude, 'LineWidth', 2, 'DisplayName', sprintf('Antenna %d', antenna_idx));
end
hold off;
grid on;
xlabel('Subcarrier Index');
ylabel('Channel Magnitude');
title('Channel Magnitude Response for Multiple Antennas');
legend;

%% Simulation Length and SNR Points
% Set the length of the simulation in terms of the number of 10ms frames. A
% large number of NFrames should be used to produce meaningful throughput
% results. Set the SNR points to simulate. The SNR for each layer is
% defined per RE, and it includes the effect of signal and noise across all
% antennas. For an explanation of the SNR definition that this example
% uses, see <docid:5g_ug#mw_37cef3ca-2f4b-433d-8d68-117a881ca5fd SNR
% Definition used in Link Simulations>.

simParameters = struct();       % Clear simParameters variable to contain all key simulation parameters 
simParameters.NFrames = 2;      % Number of 10 ms frames
simParameters.SNRIn = [-9, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]; % SNR range (dB)
%throughputResults = zeros(numel(simParameters.SNRIn), 3); % 3 for lower, current, and upper CQI

%% Channel Estimator Configuration
% The logical variable |PerfectChannelEstimator| controls channel
% estimation and synchronization behavior. When set to |true|, perfect
% channel estimation and synchronization is used. Otherwise, practical
% channel estimation and synchronization is used, based on the values of
% the received PDSCH DM-RS.

simParameters.PerfectChannelEstimator = false;

%% Simulation Diagnostics
% The variable |DisplaySimulationInformation| controls the display of
% simulation information such as the HARQ process ID used for each
% subframe. In case of CRC error, the value of the index to the RV sequence
% is also displayed.

simParameters.DisplaySimulationInformation = true;

%%
% The |DisplayDiagnostics| flag enables the plotting of the EVM per layer.
% This plot monitors the quality of the received signal after equalization.
% The EVM per layer figure shows:
%
% * The EVM per layer per slot, which shows the EVM evolving with time.
% * The EVM per layer per resource block, which shows the EVM in frequency.
%
% This figure evolves with the simulation and is updated with each slot. 
% Typically, low SNR or channel fades can result in decreased signal 
% quality (high EVM). The channel affects each layer differently,
% therefore, the EVM values may differ across layers.
%
% In some cases, some layers can have a much higher EVM than others. These
% low-quality layers can result in CRC errors. This behavior may be caused
% by low SNR or by using too many layers for the channel conditions. You
% can avoid this situation by a combination of higher SNR, lower number
% of layers, higher number of antennas, and more robust transmission 
% (lower modulation scheme and target code rate).

simParameters.DisplayDiagnostics = false;

%% Carrier and PDSCH Configuration
% Set the key parameters of the simulation. These include:
% 
% * The bandwidth in resource blocks (12 subcarriers per resource block).
% * Subcarrier spacing: 15, 30, 60, 120 (kHz)
% * Cyclic prefix length: normal or extended
% * Cell ID
% * Number of transmit and receive antennas
% 
% A substructure containing the DL-SCH and PDSCH parameters is also
% specified. This includes:
% 
% * Target code rate
% * Allocated resource blocks (PRBSet)
% * Modulation scheme: 'QPSK', '16QAM', '64QAM', '256QAM'
% * Number of layers
% * PDSCH mapping type
% * DM-RS configuration parameters
% * PT-RS configuration parameters
% 
% Other simulation wide parameters are:
% 
% * Propagation channel model delay profile (TDL or CDL)

% Set waveform type and PDSCH numerology (SCS and CP type)
simParameters.Carrier = nrCarrierConfig;         % Carrier resource grid configuration
simParameters.Carrier.NSizeGrid = 52;            % Bandwidth in number of resource blocks (51 RBs at 30 kHz SCS for 20 MHz BW)
simParameters.Carrier.SubcarrierSpacing = 30;    % 15, 30, 60, 120 (kHz)
simParameters.Carrier.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
simParameters.Carrier.NCellID = 1;               % Cell identity

% PDSCH/DL-SCH parameters
simParameters.PDSCH = nrPDSCHConfig;      % This PDSCH definition is the basis for all PDSCH transmissions in the BLER simulation
simParameters.PDSCHExtension = struct();  % This structure is to hold additional simulation parameters for the DL-SCH and PDSCH

% Define PDSCH time-frequency resource allocation per slot to be full grid (single full grid BWP)
simParameters.PDSCH.PRBSet = 0:simParameters.Carrier.NSizeGrid-1;                 % PDSCH PRB allocation
simParameters.PDSCH.SymbolAllocation = [0,simParameters.Carrier.SymbolsPerSlot];  % Starting symbol and number of symbols of each PDSCH allocation
simParameters.PDSCH.MappingType = 'A';     % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))

% Scrambling identifiers
simParameters.PDSCH.NID = simParameters.Carrier.NCellID;
simParameters.PDSCH.RNTI = 1;

% PDSCH resource block mapping (TS 38.211 Section 7.3.1.6)
simParameters.PDSCH.VRBToPRBInterleaving = 0; % Disable interleaved resource mapping
simParameters.PDSCH.VRBBundleSize = 4;

% Define the number of transmission layers to be used
simParameters.PDSCH.NumLayers = 1;            % Number of PDSCH transmission layers
simParameters.PDSCHExtension.MCSTable      = 'Table 1'; % 'Table1',...,'Table4'
% Define codeword modulation and target coding rate
% The number of codewords is directly dependent on the number of layers so ensure that 
% layers are set first before getting the codeword number
if simParameters.PDSCH.NumCodewords > 1                             % Multicodeword transmission (when number of layers being > 4)
    simParameters.PDSCH.Modulation = {'16QAM','16QAM'};             % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PDSCHExtension.TargetCodeRate = [490 490]/1024;   % Code rate used to calculate transport block sizes
else
    simParameters.PDSCH.Modulation = 'QPSK';                       % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PDSCHExtension.TargetCodeRate = 78/1024;         % Code rate used to calculate transport block sizes
end

% DM-RS and antenna port configuration (TS 38.211 Section 7.4.1.1)
simParameters.PDSCH.DMRS.DMRSPortSet = 0:simParameters.PDSCH.NumLayers-1; % DM-RS ports to use for the layers
simParameters.PDSCH.DMRS.DMRSTypeAPosition = 2;      % Mapping type A only. First DM-RS symbol position (2,3)
simParameters.PDSCH.DMRS.DMRSLength = 1;             % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
simParameters.PDSCH.DMRS.DMRSAdditionalPosition = 2; % Additional DM-RS symbol positions (max range 0...3)
simParameters.PDSCH.DMRS.DMRSConfigurationType = 2;  % DM-RS configuration type (1,2)
simParameters.PDSCH.DMRS.NumCDMGroupsWithoutData = 1;% Number of CDM groups without data
simParameters.PDSCH.DMRS.NIDNSCID = 1;               % Scrambling identity (0...65535)
simParameters.PDSCH.DMRS.NSCID = 0;                  % Scrambling initialization (0,1)

% PT-RS configuration (TS 38.211 Section 7.4.1.2)
simParameters.PDSCH.EnablePTRS = 0;                  % Enable or disable PT-RS (1 or 0)
simParameters.PDSCH.PTRS.TimeDensity = 1;            % PT-RS time density (L_PT-RS) (1, 2, 4)
simParameters.PDSCH.PTRS.FrequencyDensity = 2;       % PT-RS frequency density (K_PT-RS) (2 or 4)
simParameters.PDSCH.PTRS.REOffset = '00';            % PT-RS resource element offset ('00', '01', '10', '11')
simParameters.PDSCH.PTRS.PTRSPortSet = [];           % PT-RS antenna port, subset of DM-RS port set. Empty corresponds to lower DM-RS port number

% Reserved PRB patterns, if required (for CORESETs, forward compatibility etc)
simParameters.PDSCH.ReservedPRB{1}.SymbolSet = [];   % Reserved PDSCH symbols
simParameters.PDSCH.ReservedPRB{1}.PRBSet = [];      % Reserved PDSCH PRBs
simParameters.PDSCH.ReservedPRB{1}.Period = [];      % Periodicity of reserved resources

% Additional simulation and DL-SCH related parameters
%
% PDSCH PRB bundling (TS 38.214 Section 5.1.2.3)
simParameters.PDSCHExtension.PRGBundleSize = [];     % 2, 4, or [] to signify "wideband"
%
% HARQ process and rate matching/TBS parameters
simParameters.PDSCHExtension.XOverhead = 6*simParameters.PDSCH.EnablePTRS; % Set PDSCH rate matching overhead for TBS (Xoh) to 6 when PT-RS is enabled, otherwise 0
simParameters.PDSCHExtension.NHARQProcesses = 16;    % Number of parallel HARQ processes to use
simParameters.PDSCHExtension.EnableHARQ = false;      % Enable retransmissions for each process, using RV sequence [0,2,3,1]

% LDPC decoder parameters
% Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
simParameters.PDSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
simParameters.PDSCHExtension.MaximumLDPCIterationCount = 6;

% Define the overall transmission antenna geometry at end-points
% If using a CDL propagation channel then the integer number of antenna elements is
% turned into an antenna panel configured when the channel model object is created
simParameters.NTxAnts = 4;                        % Number of PDSCH transmission antennas (1,2,4,8,16,32,64,128,256,512,1024) >= NumLayers
if simParameters.PDSCH.NumCodewords > 1           % Multi-codeword transmission
    simParameters.NRxAnts = 8;                    % Number of UE receive antennas (even number >= NumLayers)
else
    simParameters.NRxAnts = 1;                    % Number of UE receive antennas (1 or even number >= NumLayers)
end

% Define data type ('single' or 'double') for resource grids and waveforms
simParameters.DataType = 'single';

% Define the general CDL/TDL propagation channel parameters
simParameters.DelayProfile = 'CDL-A';   % Use CDL-C model (Urban macrocell model)
simParameters.DelaySpread = 300e-9;
%simParameters.MaximumDopplerShift = 5; 

% Cross-check the PDSCH layering against the channel geometry 
validateNumLayers(simParameters);
%% csi-rs
% Add this section for CSI reporting configuration


% Bandwidth part (BWP) configuration
NStartBWP = 0;
NSizeBWP = 52;

% CSI-RS configuration
csirs = nrCSIRSConfig;
csirs.CSIRSType = {'nzp','nzp','nzp'};  % Non-zero-power CSI-RS
csirs.RowNumber = [4 4 4];
csirs.NumRB = 51;                       % Resource blocks for CSI-RS
csirs.RBOffset = 0;
csirs.CSIRSPeriod = [4 0];
csirs.SymbolLocations = {0, 0, 0};
csirs.SubcarrierLocations = {0, 4, 8};
csirs.Density = {'one','one','one'};

% Get number of CSI-RS ports
csirsPorts = csirs.NumCSIRSPorts(1);

% Get CDM lengths corresponding to configured CSI-RS resources
cdmLengths = getCDMLengths(csirs);

% Initialize the practical timing offset as zero. It is updated in the
% slots when the correlation is strong for practical synchronization.
offsetPractical = 0;

% Configure the number of antennas for CSI-RS
nTxAnts = csirs.NumCSIRSPorts(1);  % CSI-RS ports correspond to transmit antennas
nRxAnts = 1;                       % Number of receive antennas

% Validate the CSI-RS configuration
validateCSIRSConfig(carrier, csirs, nTxAnts);

% CSI reporting configuration (CQI, PMI, RI computation)
reportConfig.NStartBWP = NStartBWP;
reportConfig.NSizeBWP = NSizeBWP;
reportConfig.CQITable = 'table1';
reportConfig.CodebookType = 'Type1SinglePanel';
reportConfig.PanelDimensions = [2 1];
reportConfig.CQIMode = 'Subband';
reportConfig.PMIMode = 'Subband';
reportConfig.SubbandSize = 4;

% Initialize storage for CSI feedback
cqiPracticalPerSlot = [];
subbandCQIPractical = [];
pmiPracticalPerSlot = struct('i1',[],'i2',[]);
SINRPerSubbandPerCWPractical = [];
cqiPerfectPerSlot = [];
subbandCQIPerfect = [];
pmiPerfectPerSlot = struct('i1',[],'i2',[]);
SINRPerSubbandPerCWPerfect = [];
riPracticalPerSlot = [];
riPerfectPerSlot = [];

% CQI to MCS and SNR mapping from the table
CQI_MCS_Table = struct( ...
    'CQI',   num2cell([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]), ...
    'Modulation', {'QPSK', 'QPSK', 'QPSK', 'QPSK', 'QPSK', 'QPSK', '16QAM', '16QAM', '16QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM'}, ...
    'CodeRate', num2cell([78/1024, 193/1024, 449/1024, 378/1024, 490/1024, 616/1024, 466/1024, 567/1024, 666/1024, 772/1024, 873/1024, 711/1024, 797/1024, 885/1024, 948/1024]), ...
    'SNRdB', num2cell([-6.3, -5.6, -5.2, -3.6, -1.4, 1.3, 3.9, 5.3, 8.1, 9.8, 11.7, 13.6, 15.8, 18.8, 22.5]));
snrValues = simParameters.SNRIn;
numSNR = numel(snrValues);
lowerCQIThroughput = zeros(numSNR, 1);  
currentCQIThroughput = zeros(numSNR, 1); 
upperCQIThroughput = zeros(numSNR, 1); 
%% 
% The simulation relies on various pieces of information about the baseband 
% waveform, such as sample rate.

waveformInfo = nrOFDMInfo(simParameters.Carrier); % Get information about the baseband waveform after OFDM modulation step

%% Propagation Channel Model Construction
% Create the channel model object for the simulation. Both CDL and TDL channel 
% models are supported [ <#13 5> ].

% Constructed the CDL or TDL channel model object
if contains(simParameters.DelayProfile,'CDL','IgnoreCase',true)
    
    channel = nrCDLChannel; % CDL channel object
    
    % Turn the number of antennas into antenna panel array layouts. If
    % NTxAnts is not one of (1,2,4,8,16,32,64,128,256,512,1024), its value
    % is rounded up to the nearest value in the set. If NRxAnts is not 1 or
    % even, its value is rounded up to the nearest even number.
    channel = hArrayGeometry(channel,simParameters.NTxAnts,simParameters.NRxAnts);
    simParameters.NTxAnts = prod(channel.TransmitAntennaArray.Size);
    simParameters.NRxAnts = prod(channel.ReceiveAntennaArray.Size);
else
    channel = nrTDLChannel; % TDL channel object
    
    % Set the channel geometry
    channel.NumTransmitAntennas = simParameters.NTxAnts;
    channel.NumReceiveAntennas = simParameters.NRxAnts;
end

% Assign simulation channel parameters and waveform sample rate to the object
channel.DelayProfile = simParameters.DelayProfile;
channel.DelaySpread = simParameters.DelaySpread;
%channel.MaximumDopplerShift = simParameters.MaximumDopplerShift;
channel.SampleRate = waveformInfo.SampleRate;

%% 
% Get the maximum number of delayed samples by a channel multipath
% component. This is calculated from the channel path with the largest
% delay and the implementation delay of the channel filter. This is
% required later to flush the channel filter to obtain the received signal.

chInfo = info(channel);
maxChDelay = chInfo.MaximumChannelDelay;

%% Processing Loop
% To determine the throughput at each SNR point, analyze the PDSCH data 
% per transmission instance using the following steps:
% 
% * _Update current HARQ process._ Check the transmission status for the
% given HARQ process to determine whether a retransmission is required. If
% that is not the case then generate new data.
% * _Resource grid generation._ Perform channel coding by calling the
% <docid:5g_ref#mw_sysobj_nrDLSCH nrDLSCH> System object. The object
% operates on the input transport block and keeps an internal copy of the
% transport block in case a retransmission is required. Modulate the coded
% bits on the PDSCH by using the <docid:5g_ref#mw_function_nrPDSCH nrPDSCH>
% function. Then apply precoding to the resulting signal.
% * _Waveform generation._ OFDM modulate the generated grid.
% * _Noisy channel modeling._ Pass the waveform through a CDL or TDL
% fading channel. Add AWGN. For an explanation of the SNR definition that
% this example uses, see
% <docid:5g_ug#mw_37cef3ca-2f4b-433d-8d68-117a881ca5fd SNR Definition used
% in Link Simulations>.
% * _Perform synchronization and OFDM demodulation._ For perfect
% synchronization, reconstruct the channel impulse response to synchronize
% the received waveform. For practical synchronization, correlate the
% received waveform with the PDSCH DM-RS. Then OFDM demodulate the
% synchronized signal.
% * _Perform channel estimation._ For perfect channel estimation,
% reconstruct the channel impulse response and perform OFDM demodulation.
% For practical channel estimation, use the PDSCH DM-RS.
% * _Perform equalization and CPE compensation._ MMSE equalize the
% estimated channel. Estimate the common phase error (CPE) by using the
% PT-RS symbols, then correct the error in each OFDM symbol within the
% range of reference PT-RS OFDM symbols.
% * _Precoding matrix calculation._ Generate the precoding matrix W for the
% next transmission by using singular value decomposition (SVD).
% * _Decode the PDSCH._ To obtain an estimate of the received codewords,
% demodulate and descramble the recovered PDSCH symbols for all transmit
% and receive antenna pairs, along with a noise estimate, by using the
% <docid:5g_ref#mw_function_nrPDSCHDecode nrPDSCHDecode> function.
% * _Decode the downlink shared channel (DL-SCH) and update HARQ process
% with the block CRC error._ Pass the vector of decoded soft bits to the
% <docid:5g_ref#mw_sysobj_nrDLSCHDecoder nrDLSCHDecoder> System object. The
% object decodes the codeword and returns the block CRC error used to
% determine the throughput of the system.

% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(simParameters.SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(simParameters.SNRIn),1);

% Set up redundancy version (RV) sequence for all HARQ processes
if simParameters.PDSCHExtension.EnableHARQ
    % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
    % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
    % should be used. If self-decodability is the priority, it should be
    % taken into account that the upper limit of the code rate at which
    % each RV is self-decodable is in the following order: 0>3>2>1
    rvSeq = [0 2 3 1];
else
    % HARQ disabled - single transmission with RV=0, no retransmissions
    rvSeq = 0; 
end

% Create DL-SCH encoder system object to perform transport channel encoding
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;

% Create DL-SCH decoder system object to perform transport channel decoding
% Use layered belief propagation for LDPC decoding, with half the number of
% iterations as compared to the default for belief propagation decoding
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;
decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

% Initialize structures to hold the data
snrData = [];
channelEstData = [];
noiseEstData = [];
berData = [];
throughputResults = zeros(numSNR, 10);

for snrIdx = 1:numel(simParameters.SNRIn)      % comment out for parallel computing
% parfor snrIdx = 1:numel(simParameters.SNRIn) % uncomment for parallel computing
% To reduce the total simulation time, you can execute this loop in
% parallel by using the Parallel Computing Toolbox. Comment out the 'for'
% statement and uncomment the 'parfor' statement. If the Parallel Computing
% Toolbox is not installed, 'parfor' defaults to normal 'for' statement.
% Because parfor-loop iterations are executed in parallel in a
% nondeterministic order, the simulation information displayed for each SNR
% point can be intertwined. To switch off simulation information display,
% set the 'displaySimulationInformation' variable above to false

    % Reset the random number generator so that each SNR point will
    % experience the same noise realization
    rng('default');
    
    % Take full copies of the simulation-level parameter structures so that they are not 
    % PCT broadcast variables when using parfor 
    simLocal = simParameters;
    waveinfoLocal = waveformInfo;
    
    % Take copies of channel-level parameters to simplify subsequent parameter referencing 
    carrier = simLocal.Carrier;
    pdsch = simLocal.PDSCH;
    pdschextra = simLocal.PDSCHExtension;
    decodeDLSCHLocal = decodeDLSCH;  % Copy of the decoder handle to help PCT classification of variable
    decodeDLSCHLocal.reset();        % Reset decoder at the start of each SNR point
    pathFilters = [];
     
    % Prepare simulation for new SNR point
    SNRdB = simLocal.SNRIn(snrIdx);
    fprintf('\nSimulating transmission scheme 1 (%dx%d) and SCS=%dkHz with %s channel at %gdB SNR for %d 10ms frame(s)\n', ...
        simLocal.NTxAnts,simLocal.NRxAnts,carrier.SubcarrierSpacing, ...
        simLocal.DelayProfile,SNRdB,simLocal.NFrames);

    [lowerCQI, currentCQI, upperCQI] = getAdjacentCQIs(SNRdB, CQI_MCS_Table);
    fprintf('For SNR = %.2f dB, using CQI = %d (lower), CQI = %d (current), and CQI = %d (upper)\n', SNRdB, lowerCQI.CQI, currentCQI.CQI, upperCQI.CQI);
    CQIs = [lowerCQI, currentCQI, upperCQI];
cqiRange = 1:3;
numCQI = numel(cqiRange);


    for cqiIdx = 1:numCQI
    currentCQI = CQI_MCS_Table([CQI_MCS_Table.CQI] == cqiRange(cqiIdx));
    fprintf('\nSimulating for CQI = %d\n', cqiRange(cqiIdx));
    simParameters.PDSCH.Modulation = currentCQI.Modulation;
    simParameters.PDSCHExtension.TargetCodeRate = currentCQI.CodeRate;
%    simParameters.PDSCH.Modulation = CQIs(cqiIdx).Modulation;
 %   simParameters.PDSCHExtension.TargetCodeRate = CQIs(cqiIdx).CodeRate;
    % Specify the fixed order in which we cycle through the HARQ process IDs
    harqSequence = 0:pdschextra.NHARQProcesses-1;

    % Initialize the state of all HARQ processes
    harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);
        
    % Reset the channel so that each CQI point will experience the same
    % channel realization
    reset(channel);
    
    % Total number of slots in the simulation period
    NSlots = simLocal.NFrames * carrier.SlotsPerFrame;
    
    % Obtain a precoding matrix (wtx) to be used in the transmission of the
    % first transport block
    estChannelGridAnts = getInitialChannelEstimate(carrier,simLocal.NTxAnts,channel,simLocal.DataType);
    initialChannelResponse = estChannelGridAnts;
    newWtx = hSVDPrecoders(carrier,pdsch,estChannelGridAnts,pdschextra.PRGBundleSize);
    
    % Timing offset, updated in every slot for perfect synchronization and
    % when the correlation is strong for practical synchronization
    offset = 0;

    % Initialize variables for BER calculation
    totalBits = 0;
    totalBitErrors = 0;
totSlotsBinaryVec = zeros(1,NSlots);

    % Loop over the entire waveform length
    for nslot = 0:NSlots-1
              % % Create carrier resource grid for one slot
              %    csirsSlotGrid = nrResourceGrid(carrier,csirsPorts);
              % % Update slot number in carrier configuration object
              %    carrier.NSlot = nslot;
              % % Generate CSI-RS indices and symbols
              %    csirsInd = nrCSIRSIndices(carrier,csirs);
              %    csirsSym = nrCSIRS(carrier,csirs);
              % 
              % % Map CSI-RS to slot grid
              %    csirsSlotGrid(csirsInd) = csirsSym;
                reset(channel);
        % Calculate the transport block sizes for the transmission in the slot
        [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschextra.TargetCodeRate,pdschextra.XOverhead);
        
        % HARQ processing
        for cwIdx = 1:pdsch.NumCodewords
            % If new data for current process and codeword then create a new DL-SCH transport block
            if harqEntity.NewData(cwIdx) 
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCHLocal,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end
                                
        % Encode the DL-SCH transport blocks
        codedTrBlocks = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers, ...
            pdschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
        % Get precoding matrix (wtx) calculated in previous slot
        wtx = newWtx;
        %wtx = eye(csirsPorts,4);

        % Create resource grid for a slot
        pdschGrid = nrResourceGrid(carrier,simLocal.NTxAnts,OutputDataType=simLocal.DataType);
        


        % PDSCH reserved REs for CSI-RS
        pdsch.ReservedRE = csirsInd-1; % 0-based indices

        % PDSCH modulation and precoding
        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlocks);
        [pdschAntSymbols,pdschAntIndices] = nrPDSCHPrecode(carrier,pdschSymbols,pdschIndices,wtx);
        
  
        % PDSCH mapping in grid associated with PDSCH transmission period
        pdschGrid(pdschAntIndices) = pdschAntSymbols;
        
        % PDSCH DM-RS precoding and mapping
        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
        [dmrsAntSymbols,dmrsAntIndices] = nrPDSCHPrecode(carrier,dmrsSymbols,dmrsIndices,wtx);
        pdschGrid(dmrsAntIndices) = dmrsAntSymbols;

        % PDSCH PT-RS precoding and mapping
        ptrsSymbols = nrPDSCHPTRS(carrier,pdsch);
        ptrsIndices = nrPDSCHPTRSIndices(carrier,pdsch);
        [ptrsAntSymbols,ptrsAntIndices] = nrPDSCHPrecode(carrier,ptrsSymbols,ptrsIndices,wtx);
        pdschGrid(ptrsAntIndices) = ptrsAntSymbols;
 


        %% Apply predistortion using the inverse channel matrix
%       if nslot>0
%         H_est = nextEstChannelGridPorts(:,:,1:4,1:4); % siso
% H_est = H_est / sqrt(612);
% 
%         % apply predistortion
%         pdschGrid = pdschGrid ./ H_est;
% 
% pdschGrid = pdschGrid /sqrt(612);
% 
%       end

%% Apply water-filling 
if nslot > 0
    % Get practical SINR per subband per codeword for the previous slot (linear scale)
    currentSNR = SINRPerSubbandPerCWPractical(2:14, :, nslot);   % [NumSubbands x Ncw]
    numSubbands = size(currentSNR, 1);

    % 1) Compute a representative SNR for each subband and convert to alpha = N0/g = 1/SNR
    % Use the nanmean across codewords as the subband SNR (for single codeword, it's just that value)
    snr_sb = nanmean(currentSNR, 2);                  % [NumSubbands x 1] representative SNR per subband
    % Handle invalid or non-positive SNR values by assigning a very large alpha (to minimize their power allocation)
    alpha = 1 ./ max(snr_sb, eps);                    % alpha_k = 1/SNR_k (treat eps for zero SNR to avoid Inf)
    alpha(~isfinite(alpha) | alpha <= 0) = 1e9;       % Extremely large alpha for invalid SNR (acts like a very bad channel)

    % 2) Water-filling algorithm to determine relative power p_k for each subband
    Ptot_rel = 1.0;                                   % Total relative power (normalize to 1)
    p_rel = zeros(numSubbands, 1);                    % Initialize power allocation vector
    % Determine set of subbands to include in water-filling
    include = true(numSubbands, 1);
    while true
        idx = find(include);
        n_inc = numel(idx);
        if n_inc == 0
            break; % no subbands to allocate (edge case)
        end
        % Compute water level (mu) for currently included subbands
        mu = (Ptot_rel + sum(alpha(idx))) / n_inc;
        % Identify subbands whose alpha is higher than the water level (these will get zero power)
        newInclude = include;
        newInclude(alpha > mu) = false;
        if all(newInclude == include)
            % No change in inclusion set, break out
            break;
        end
        include = newInclude;
    end
    % Compute power allocation for included subbands at final water level
    finalIdx = find(include);
    if ~isempty(finalIdx)
        mu = (Ptot_rel + sum(alpha(finalIdx))) / numel(finalIdx);
        p_rel(finalIdx) = max(mu - alpha(finalIdx), 0);
    end

    % 3) Ensure no subband is completely zero-powered (to avoid decoding issues with wideband coding)
    zeroIdx = (p_rel == 0);
    if any(zeroIdx)
        % Assign a tiny fraction of power to subbands that would otherwise get 0
        p_rel(zeroIdx) = 1e-3;  % small epsilon power
    end
    % Normalize the power allocation so that the total remains 1 (relative)
    if sum(p_rel) > 0
        p_rel = p_rel / sum(p_rel);
    end

    % 4) Scale the PDSCH grid per subband according to the power allocation
    % Divide the subcarriers into equal subbands
    nSC = size(pdschGrid, 1);
    subbandSize = floor(nSC / numSubbands);
    totalPowerBefore = sum(abs(pdschGrid(:)).^2);
    for sb = 1:numSubbands
        % Determine subcarrier indices for this subband
        scStart = (sb-1)*subbandSize + 1;
        if sb < numSubbands
            scEnd = sb * subbandSize;
        else
            scEnd = nSC;  % last subband takes any remaining subcarriers
        end
        subIdx = scStart:scEnd;
        % Current power in this subband
        P_before_sb = sum(abs(pdschGrid(subIdx, :, :)).^2, 'all');
        % Target power for this subband based on water-filling allocation
        P_target_sb = totalPowerBefore * p_rel(sb);
        % Compute scaling factor for this subband
        if P_before_sb > 0
            scaleFactor = sqrt(P_target_sb / P_before_sb);
        else
            scaleFactor = sqrt(P_target_sb / (length(subIdx) * simLocal.NTxAnts));  % approximate if no signal present
        end
        % Apply scaling to the PDSCH grid for this subband
        pdschGrid(subIdx, :, :) = pdschGrid(subIdx, :, :) * scaleFactor;
    end
    % Final sanity check: normalize total power (should be very close to original)
    totalPowerAfter = sum(abs(pdschGrid(:)).^2);
    if totalPowerAfter > 0
        pdschGrid = pdschGrid * sqrt(totalPowerBefore / totalPowerAfter);
    end
end




 % CSI-RS mapping to the slot resource grid
        [csirsInd,csirsInfo] = nrCSIRSIndices(carrier,csirs);
        csirsSym = nrCSIRS(carrier,csirs);
        pdschGrid(csirsInd) = csirsSym;
        csirsTransmission = ~isempty(csirsInd);
 



%% 
for txAnt = 1:size(pdschGrid, 3)
    pdschGrid(:, :, txAnt) = h(:, :, txAnt) .* pdschGrid(:, :, txAnt);
end        
        % OFDM modulation

        rxWaveform = nrOFDMModulate(carrier,pdschGrid);

        % Pass data through channel model. Append zeros at the end of the
        % transmitted waveform to flush channel content. These zeros take
        % into account any delay introduced in the channel. This is a mix
        % of multipath delay and implementation delay. This value may 
        % change depending on the sampling rate, delay profile and delay
        % spread
       % txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))]; %#ok<AGROW>
        % [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
        % % 
        %  release(channel); 
% numRxAnts = 1; 
% numTxAnts = size(txWaveform, 2); 
% numSamples = size(txWaveform, 1);
% rxWaveform = zeros(numSamples, numRxAnts);
% h_time = ifft(h, [], 1); %
% %h_time = reshape(h_time, [], numRxAnts, numTxAnts); 
% 
% for rxAnt = 1:numRxAnts
%     for txAnt = 1:numTxAnts
%         rxWaveform(:, rxAnt) = rxWaveform(:, rxAnt) + conv(txWaveform(:, txAnt), squeeze(h_time(:, rxAnt, txAnt)), 'same');
%     end
% end
        % Add AWGN to the received time domain waveform
        % Normalize noise power by the IFFT size used in OFDM modulation,
        % as the OFDM modulator applies this normalization to the
        % transmitted waveform. Also normalize by the number of receive
        % antennas, as the channel model applies this normalization to the
        % received waveform, by default
signalPower = mean(abs(rxWaveform).^2, 'all');
N0 = signalPower / (10^(SNRdB/10));
noise = sqrt(N0/2) * (randn(size(rxWaveform)) + 1i * randn(size(rxWaveform)));
rxWaveform = rxWaveform + noise;

        % if (simLocal.PerfectChannelEstimator)
        %     % Perfect synchronization. Use information provided by the
        %     % channel to find the strongest multipath component
        %     pathFilters = getPathFilters(channel);
        %     [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
        % else
        %     % Practical synchronization. Correlate the received waveform
        %     % with the PDSCH DM-RS to give timing offset estimate 't' and
        %     % correlation magnitude 'mag'. The function
        %     % hSkipWeakTimingOffset is used to update the receiver timing
        %     % offset. If the correlation peak in 'mag' is weak, the current
        %     % timing estimate 't' is ignored and the previous estimate
        %     % 'offset' is used
        %     [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols); 
        %     offset = hSkipWeakTimingOffset(offset,t,mag);
        %     % Display a warning if the estimated timing offset exceeds the
        %     % maximum channel delay
        %     if offset > maxChDelay
        %         warning(['Estimated timing offset (%d) is greater than the maximum channel delay (%d).' ...
        %             ' This will result in a decoding failure. This may be caused by low SNR,' ...
        %             ' or not enough DM-RS symbols to synchronize successfully.'],offset,maxChDelay);
        %     end
        % end
        % rxWaveform = rxWaveform(1+offset:end,:); 


        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
        [K,L,R] = size(rxGrid);
        if (L < carrier.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,carrier.SymbolsPerSlot-L,R));
        end

        if (simLocal.PerfectChannelEstimator)
            % Perfect channel estimation, using the value of the path gains
            % provided by the channel. This channel estimate does not
            % include the effect of transmitter precoding
            estChannelGridAnts = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);

            % Get perfect noise estimate (from the noise realization)
            noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));
            % Get PDSCH resource elements from the received grid and 
            % channel estimate
            [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,estChannelGridAnts);

            % Apply precoding to channel estimate
            pdschHest = nrPDSCHPrecode(carrier,pdschHest,pdschHestIndices,permute(wtx,[2 1 3]));
        else
            % Practical channel estimation between the received grid and
            % each transmission layer, using the PDSCH DM-RS for each
            % layer. This channel estimate includes the effect of
            % transmitter precoding
            [estChannelGridPorts,noiseEst] = hSubbandChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,pdschextra.PRGBundleSize,'CDMLengths',pdsch.DMRS.CDMLengths); 
nextEstChannelGridPorts=estChannelGridPorts;

            % Calculate the inverse of the estimated channel matrix
            H_est_inv = pagepinv(estChannelGridPorts);
            %H_est_inv = reshape(H_est_inv, size(pdschAntSymbols));

            % Average noise estimate across PRGs and layers
            noiseEst = mean(noiseEst,'all');

            % Get PDSCH resource elements from the received grid and
            % channel estimate
            [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGridPorts);

            % Remove precoding from estChannelGridPorts to get channel
            % estimate w.r.t. antennas
            estChannelGridAnts = precodeChannelEstimate(carrier,estChannelGridPorts,conj(wtx));

        end
%%  %% CSI-RS

    % Consider only the NZP-CSI-RS symbols and indices for channel estimation
    nzpCSIRSSym = csirsSym(csirsSym ~= 0);
    nzpCSIRSInd = csirsInd(csirsSym ~= 0);

    % Calculate practical channel estimate. Use a time-averaging window
    % that covers all the transmitted CSI-RS symbols.
    [PracticalHest,nVarPractical] = nrChannelEstimate(carrier,rxGrid, ...
        nzpCSIRSInd,nzpCSIRSSym,'CDMLengths',cdmLengths,'AveragingWindow',[0 5]);

    % Perform perfect channel estimation
    %PerfectHest = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offsetPerfect,sampleTimes);
    % Get perfect noise estimate value from noise realization
    if ~isempty(nzpCSIRSInd)
        % Set the totSlotsBinaryVec value corresponding to the slot
        % index where NZP-CSI-RS is present to 1
        totSlotsBinaryVec(nslot+1) = 1;

        % Calculate the RI value using practical channel estimate
        numLayersPractical = hRISelect(carrier,csirs,reportConfig,PracticalHest,nVarPractical,'MaxSE');

        % Calculate CQI and PMI values using practical channel estimate
        [cqiPractical,pmiPractical,cqiInfoPractical,pmiInfoPractical] = hCQISelect(carrier,csirs, ...
            reportConfig,numLayersPractical,PracticalHest,nVarPractical);
        numCodeWordsPr = size(cqiPractical,2);
        numSBs = size(cqiPractical,1);

        % Store CQI, PMI, RI, and subband SINR values of each slot for the
        % practical channel estimation scenario. Because the number of
        % codewords can vary based on the rank, append NaNs to the
        % CQI-related variables to account for the second codeword
        % information in the slots where only one codeword is present.
        riPracticalPerSlot(1,nslot+1) = numLayersPractical; %#ok<SAGROW> 
        cqiPracticalPerSlot(:,:,nslot+1) = [cqiPractical NaN(numSBs,2-numCodeWordsPr)]; %#ok<SAGROW> 
        pmiPracticalPerSlot(nslot+1) = pmiPractical;
        subbandCQIPractical(:,:,nslot+1) = [cqiInfoPractical.SubbandCQI NaN(numSBs,2-numCodeWordsPr)]; %#ok<SAGROW> 
        SINRPerSubbandPerCWPractical(:,:,nslot+1) = [cqiInfoPractical.SINRPerSubbandPerCW NaN(numSBs,2-numCodeWordsPr)]; %#ok<SAGROW> 

        % Calculate the RI value using perfect channel estimate
        %numLayersPerfect = hRISelect(carrier,csirs,reportConfig,PerfectHest,nVarPerfect,'MaxSE');

        % Calculate CQI and PMI values using perfect channel estimate
        %[cqiPerfect,pmiPerfect,cqiInfoPerfect,pmiInfoPerfect] = hCQISelect(carrier,csirs, ...
        %    reportConfig,numLayersPerfect,PerfectHest,nVarPerfect);
        %numCodeWordsPe = size(cqiPerfect,2);

        % Store CQI, PMI, RI, and subband SINR values of each slot for the
        % perfect channel estimation scenario. Because the number of
        % codewords can vary based on the rank, append NaNs to the
        % CQI-related variables to account for the second codeword
        % information in the slots where only one codeword is present.
        %riPerfectPerSlot(1,nslot+1) = numLayersPerfect; %#ok<SAGROW> 
        %cqiPerfectPerSlot(:,:,nslot+1) = [cqiPerfect NaN(numSBs,2-numCodeWordsPe)]; %#ok<SAGROW> 
        %subbandCQIPerfect(:,:,nslot+1) = [cqiInfoPerfect.SubbandCQI NaN(numSBs,2-numCodeWordsPe)]; %#ok<SAGROW> 
        %pmiPerfectPerSlot(nslot+1) = pmiPerfect;
        %SINRPerSubbandPerCWPerfect(:,:,nslot+1) = [cqiInfoPerfect.SINRPerSubbandPerCW NaN(numSBs,2-numCodeWordsPe)]; %#ok<SAGROW> 
    end
    % Get the active slot numbers (1-based) in which NZP-CSI-RS is present
    activeSlotNum = find(totSlotsBinaryVec);

    % Fill the CQI, PMI and RI variables with NaNs in the slots where NZP-CSI-RS is
    % absent according to codebook type
    [cqiPracticalPerSlot,subbandCQIPractical,pmiPracticalPerSlot,SINRPerSubbandPerCWPractical, ...
     cqiPerfectPerSlot,subbandCQIPerfect,pmiPerfectPerSlot,SINRPerSubbandPerCWPerfect,riPracticalPerSlot, ...
     riPerfectPerSlot] = fillInactiveSlots(cqiPracticalPerSlot,subbandCQIPractical, ...
    pmiPracticalPerSlot,SINRPerSubbandPerCWPractical,cqiPerfectPerSlot,subbandCQIPerfect,pmiPerfectPerSlot, ...
    SINRPerSubbandPerCWPerfect,riPracticalPerSlot,riPerfectPerSlot,reportConfig,totSlotsBinaryVec,activeSlotNum);
%% 


        % Equalization
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        % Common phase error (CPE) compensation
        if ~isempty(ptrsIndices)
            % Initialize temporary grid to store equalized symbols
            tempGrid = nrResourceGrid(carrier,pdsch.NumLayers);

            % Extract PT-RS symbols from received grid and estimated
            % channel grid
            [ptrsRx,ptrsHest,~,~,ptrsHestIndices,ptrsLayerIndices] = nrExtractResources(ptrsIndices,rxGrid,estChannelGridAnts,tempGrid);
            ptrsHest = nrPDSCHPrecode(carrier,ptrsHest,ptrsHestIndices,permute(wtx,[2 1 3]));

            % Equalize PT-RS symbols and map them to tempGrid
            ptrsEq = nrEqualizeMMSE(ptrsRx,ptrsHest,noiseEst);
            tempGrid(ptrsLayerIndices) = ptrsEq;

            % Estimate the residual channel at the PT-RS locations in
            % tempGrid
            cpe = nrChannelEstimate(tempGrid,ptrsIndices,ptrsSymbols);

            % Sum estimates across subcarriers, receive antennas, and
            % layers. Then, get the CPE by taking the angle of the
            % resultant sum
            cpe = angle(sum(cpe,[1 3 4]));

            % Map the equalized PDSCH symbols to tempGrid
            tempGrid(pdschIndices) = pdschEq;

            % Correct CPE in each OFDM symbol within the range of reference
            % PT-RS OFDM symbols
            symLoc = pdschIndicesInfo.PTRSSymbolSet(1)+1:pdschIndicesInfo.PTRSSymbolSet(end)+1;
            tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1i*cpe(symLoc));

            % Extract PDSCH symbols
            pdschEq = tempGrid(pdschIndices);
        end

        % Decode PDSCH physical channel
        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);
        
        % Display EVM per layer, per slot and per RB
        if (simLocal.DisplayDiagnostics)
            plotLayerEVM(NSlots,nslot,pdsch,size(pdschGrid),pdschIndices,pdschSymbols,pdschEq);
        end
        
        % Scale LLRs by CSI
        csi = nrLayerDemap(csi); % CSI layer demapping
        for cwIdx = 1:pdsch.NumCodewords
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % bits per symbol
            csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % expand by each bit per symbol
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale by CSI
        end
        
        % Decode the DL-SCH transport channel
        decodeDLSCHLocal.TransportBlockLength = trBlkSizes;
        [decbits,blkerr] = decodeDLSCHLocal(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        
        % Calculate bit errors for current slot
        for cwIdx = 1:pdsch.NumCodewords
            if harqEntity.NewData(cwIdx)
            trBlk = getTransportBlock(encodeDLSCH, cwIdx-1, harqEntity.HARQProcessID);
            totalBits = totalBits + length(trBlk);
            bitErrors = sum(decbits(:, cwIdx) ~= trBlk); % Correct indexing for array
            totalBitErrors = totalBitErrors + bitErrors;
            fprintf('Slot %d, CW %d: Total Bits = %d, Bit Errors = %d\n', nslot, cwIdx, totalBits, bitErrors);

            end
        end

        % Store values to calculate throughput
        simThroughput(snrIdx) = simThroughput(snrIdx) + sum(~blkerr .* trBlkSizes);
        maxThroughput(snrIdx) = maxThroughput(snrIdx) + sum(trBlkSizes);
%how cellular do error detection/correction of the block
%how many blocks in 2 frames
        % Update current process with CRC error and advance to next process
        procstatus = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschIndicesInfo.G);
        if (simLocal.DisplaySimulationInformation)
            fprintf('\n(%3.2f%%) NSlot=%d, %s',100*(nslot+1)/NSlots,nslot,procstatus);
        end

        % Get precoding matrix for next slot
        newWtx = hSVDPrecoders(carrier,pdsch,estChannelGridAnts,pdschextra.PRGBundleSize);

     end
    
    % Display the results dynamically in the command window
    if (simLocal.DisplaySimulationInformation)
        fprintf('\n');
    end
    fprintf('\nThroughput(Mbps) for %d frame(s) = %.4f\n',simLocal.NFrames,1e-6*simThroughput(snrIdx)/(simLocal.NFrames*10e-3));
    fprintf('Throughput(%%) for %d frame(s) = %.4f\n',simLocal.NFrames,simThroughput(snrIdx)*100/maxThroughput(snrIdx));
    throughputResults(snrIdx, cqiIdx) = 1e-6 * simThroughput(snrIdx) / (simLocal.NFrames * 10e-3);

    % Calculate and display BER
    ber = totalBitErrors / totalBits;
    fprintf('Bit Error Rate (BER) = %.6f\n', ber);
 
    end
end

%% Results & Save data
% Display the measured throughput. This is calculated as the percentage of
% the maximum possible throughput of the link given the available resources
% for data transmission.
%writetable(resultsTable, 'simulation_results.xlsx');

% slotDuration = 1e-3;  
% totalTime = NSlots * slotDuration; 
% 
% throughputMbps = (simThroughput / totalTime) / 1e6;  
% 
% figure;
% plot(simParameters.SNRIn, throughputMbps, 'o-.')  
% xlabel('SNR (dB)'); ylabel('Throughput (Mbps)'); grid on;
% title(sprintf('%s (%dx%d) / NRB=%d / SCS=%dkHz', ...
%               simParameters.DelayProfile,simParameters.NTxAnts,simParameters.NRxAnts, ...
%               simParameters.Carrier.NSizeGrid,simParameters.Carrier.SubcarrierSpacing));
% 
% simResults.simParameters = simParameters;
% simResults.simThroughput = throughputMbps; 
% simResults.maxThroughput = maxThroughput;
figure;
hold on;
for cqiIdx = 1:numCQI
    plot(snrValues, throughputResults(:, cqiIdx), 'o-', 'DisplayName', sprintf('CQI = %d', cqiRange(cqiIdx)));
end
xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
legend show;
grid on;
title('Throughput vs SNR for Different CQI Levels');
hold off;
save('results_withoutwater_CQI7-9.mat', 'throughputResults', 'snrValues');

%%
% The figure below shows throughput results obtained simulating 10000
% subframes (|NFrames = 1000|, |SNRIn = -18:2:16|).
%
% <<../longRunThroughput.png>>
%

%% Selected Bibliography
% # 3GPP TS 38.211. "NR; Physical channels and modulation."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.212. "NR; Multiplexing and channel coding." 
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.213. "NR; Physical layer procedures for control." 
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.214. "NR; Physical layer procedures for data."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TR 38.901. "Study on channel model for frequencies from 0.5 to 100
% GHz."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.

%% Local Functions

function validateNumLayers(simParameters)
% Validate the number of layers, relative to the antenna geometry

    numlayers = simParameters.PDSCH.NumLayers;
    ntxants = simParameters.NTxAnts;
    nrxants = simParameters.NRxAnts;
    antennaDescription = sprintf('min(NTxAnts,NRxAnts) = min(%d,%d) = %d',ntxants,nrxants,min(ntxants,nrxants));
    if numlayers > min(ntxants,nrxants)
        error('The number of layers (%d) must satisfy NumLayers <= %s', ...
            numlayers,antennaDescription);
    end
    
    % Display a warning if the maximum possible rank of the channel equals
    % the number of layers
    if (numlayers > 2) && (numlayers == min(ntxants,nrxants))
        warning(['The maximum possible rank of the channel, given by %s, is equal to NumLayers (%d).' ...
            ' This may result in a decoding failure under some channel conditions.' ...
            ' Try decreasing the number of layers or increasing the channel rank' ...
            ' (use more transmit or receive antennas).'],antennaDescription,numlayers); %#ok<SPWRN>
    end

end

function estChannelGrid = getInitialChannelEstimate(carrier,nTxAnts,propchannel,dataType)
% Obtain channel estimate before first transmission. This can be used to
% obtain a precoding matrix for the first slot.

    ofdmInfo = nrOFDMInfo(carrier);
    
    chInfo = info(propchannel);
    maxChDelay =chInfo.MaximumChannelDelay;
    
    % Temporary waveform (only needed for the sizes)
    tmpWaveform = zeros((ofdmInfo.SampleRate/1000/carrier.SlotsPerSubframe)+maxChDelay,nTxAnts,dataType);
    
    % Filter through channel    
    [~,pathGains,sampleTimes] = propchannel(tmpWaveform);
    
    % Perfect timing synch    
    pathFilters = getPathFilters(propchannel);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
    
end

function estChannelGrid = precodeChannelEstimate(carrier,estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate

    [K,L,R,P] = size(estChannelGrid);
    estChannelGrid = reshape(estChannelGrid,[K*L R P]);
    estChannelGrid = nrPDSCHPrecode(carrier,estChannelGrid,reshape(1:numel(estChannelGrid),[K*L R P]),W);
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end

function plotLayerEVM(NSlots,nslot,pdsch,siz,pdschIndices,pdschSymbols,pdschEq)
% Plot EVM information

    persistent slotEVM;
    persistent rbEVM
    persistent evmPerSlot;
    
    if (nslot==0)
        slotEVM = comm.EVM;
        rbEVM = comm.EVM;
        evmPerSlot = NaN(NSlots,pdsch.NumLayers);
        figure;
    end
    evmPerSlot(nslot+1,:) = slotEVM(pdschSymbols,pdschEq);
    subplot(2,1,1);
    plot(0:(NSlots-1),evmPerSlot,'o-');
    xlabel('Slot number');
    ylabel('EVM (%)');
    legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
    title('EVM per layer per slot');

    subplot(2,1,2);
    [k,~,p] = ind2sub(siz,pdschIndices);
    rbsubs = floor((k-1) / 12);
    NRB = siz(1) / 12;
    evmPerRB = NaN(NRB,pdsch.NumLayers);
    for nu = 1:pdsch.NumLayers
        for rb = unique(rbsubs).'
            this = (rbsubs==rb & p==nu);
            evmPerRB(rb+1,nu) = rbEVM(pdschSymbols(this),pdschEq(this));
        end
    end
    plot(0:(NRB-1),evmPerRB,'x-');
    xlabel('Resource block');
    ylabel('EVM (%)');
    legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
    title(['EVM per layer per resource block, slot #' num2str(nslot)]);
    
    drawnow;
    
end
%% *Local Functions*
% This example uses these local functions to validate the CSI-RS configuration 
% object and to plot the computed CQI, PMI, and RI values.

function validateCSIRSConfig(carrier,csirs,nTxAnts)
%   Validates the CSI-RS configuration, given the carrier specific
%   configuration object, CSI-RS configuration object, and the number of
%   transmit antennas.

    % Validate the number of CSI-RS ports
    if ~isscalar(unique(csirs.NumCSIRSPorts))
        error('nr5g:InvalidCSIRSPorts', ...
            'All the CSI-RS resources must be configured to have the same number of CSI-RS ports.');
    end

    % Validate the CDM lengths
    if ~iscell(csirs.CDMType)
        cdmType = {csirs.CDMType};
    else
        cdmType = csirs.CDMType;
    end
    if (~all(strcmpi(cdmType,cdmType{1})))
        error('nr5g:InvalidCSIRSCDMTypes', ...
            'All the CSI-RS resources must be configured to have the same CDM lengths.');
    end
    if nTxAnts ~= csirs.NumCSIRSPorts(1)
        error('nr5g:InvalidNumTxAnts',['Number of transmit antennas (' num2str(nTxAnts) ...
            ') must be equal to the number of CSI-RS ports (' num2str(csirs.NumCSIRSPorts(1)) ').']);
    end

    % Check for the overlap between the CSI-RS indices
    csirsInd = nrCSIRSIndices(carrier,csirs,"OutputResourceFormat",'cell');
    numRes = numel(csirsInd);
    csirsIndAll = cell(1,numRes);
    ratioVal = csirs.NumCSIRSPorts(1)/prod(getCDMLengths(csirs));
    for resIdx = 1:numRes
        if ~isempty(csirsInd{resIdx})
            grid = nrResourceGrid(carrier,csirs.NumCSIRSPorts(1));
            [~,tempInd] = nrExtractResources(csirsInd{resIdx},grid);
            if numel(tempInd)/numel(csirsInd{resIdx}) ~= ratioVal
                error('nr5g:OverlappedCSIRSREsSingleResource',['CSI-RS indices of resource ' ...
                    num2str(resIdx) ' must be unique. Try changing the symbol or subcarrier locations.']);
            end
            csirsIndAll{resIdx} = tempInd(:);
            for idx = 1:resIdx-1
                overlappedInd = ismember(csirsIndAll{idx},csirsIndAll{resIdx});
                if any(overlappedInd)
                    error('nr5g:OverlappedCSIRSREsMultipleResources',['The resource elements of the ' ...
                        'configured CSI-RS resources must not overlap. Try changing the symbol or ' ...
                        'subcarrier locations of CSI-RS resource ' num2str(idx) ' and resource ' num2str(resIdx) '.']);
                end
            end
        end
    end
end

function cdmLengths = getCDMLengths(csirs)
%   Returns the CDM lengths, given the CSI-RS configuration object.

    CDMType = csirs.CDMType;
    if ~iscell(csirs.CDMType)
        CDMType = {csirs.CDMType};
    end
    CDMTypeOpts = {'noCDM','fd-CDM2','CDM4','CDM8'};
    CDMLengthOpts = {[1 1],[2 1],[2 2],[2 4]};
    cdmLengths = CDMLengthOpts{strcmpi(CDMTypeOpts,CDMType{1})};
end

function [cqiPracticalPerSlot,subbandCQIPractical,pmiPracticalPerSlot,SINRPerSubbandPerCWPractical,cqiPerfectPerSlot, ...
    subbandCQIPerfect,pmiPerfectPerSlot,SINRPerSubbandPerCWPerfect,riPracticalPerSlot,riPerfectPerSlot] = fillInactiveSlots(cqiPracticalPerSlot, ...
    subbandCQIPractical,pmiPracticalPerSlot,SINRPerSubbandPerCWPractical,cqiPerfectPerSlot,subbandCQIPerfect,pmiPerfectPerSlot, ...
    SINRPerSubbandPerCWPerfect,riPracticalPerSlot,riPerfectPerSlot,reportConfig,totSlotsBinaryVec,activeSlots)
%   Returns the CQI, PMI, and RI related variables filled with NaNs in the
%   slots where NZP-CSI-RS is not present according to the codebook type from
%   the report configuration structure. Note that the CQI, PMI, and RI
%   variables are returned as empty if there are no NZP-CSI-RS resources,
%   that is, no active slots in the entire simulation duration.

    % Compute the indices of the slots and the number of slots in which
    % NZP-CSI-RS is not present
    inactiveSlotIdx = ~totSlotsBinaryVec;
    numInactiveSlots = nnz(inactiveSlotIdx);
    
    if ~isempty(activeSlots)
        numCQISBs = size(cqiPracticalPerSlot,1);
    
        % Get the codebook type
        codebookType = 'Type1SinglePanel';
        if isfield(reportConfig,'CodebookType')
            codebookType = validatestring(reportConfig.CodebookType,{'Type1SinglePanel','Type1MultiPanel','Type2','eType2'},'fillInactiveSlots','CodebookType field');
        end
    
        % Fill the CQI, PMI, and RI variables with NaNs in the slots where NZP-CSI-RS is
        % not present
        cqiPracticalPerSlot(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        subbandCQIPractical(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        SINRPerSubbandPerCWPractical(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        cqiPerfectPerSlot(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        subbandCQIPerfect(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        SINRPerSubbandPerCWPerfect(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        riPracticalPerSlot(inactiveSlotIdx) = NaN;
        riPerfectPerSlot(inactiveSlotIdx) = NaN;
    
        numi1Indices = 3;
        numi2Indices = 1;
        if strcmpi(codebookType,'Type1MultiPanel')
            numi1Indices = 6;
            numi2Indices = 3;
        end
        numPMISBs = size(pmiPerfectPerSlot(activeSlots(1)).i2,2);
        [pmiPerfectPerSlot(inactiveSlotIdx),pmiPracticalPerSlot(inactiveSlotIdx)] = deal(struct('i1',NaN(1,numi1Indices),'i2',NaN(numi2Indices,numPMISBs)));
    end
end

function plotWidebandCQIAndSINR(cqiPracticalPerSlot,cqiPerfectPerSlot,SINRPerSubbandPerCWPractical,SINRPerSubbandPerCWPerfect,activeSlotNum)
%   Plots the wideband SINR and wideband CQI values for each codeword
%   across all specified active slots (1-based) (in which the CQI is
%   reported as other than NaN) for practical and perfect channel
%   estimation cases.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No CQI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    cqiPracticalPerCW = permute(cqiPracticalPerSlot(1,:,:),[1 3 2]);
    cqiPerfectPerCW = permute(cqiPerfectPerSlot(1,:,:),[1 3 2]);
    SINRPerCWPractical = permute(SINRPerSubbandPerCWPractical(1,:,:),[1 3 2]);
    SINRPerCWPerfect = permute(SINRPerSubbandPerCWPerfect(1,:,:),[1 3 2]);

    % Extract wideband CQI indices for slots where NZP-CSI-RS is present
    cqiPracticalPerCWActiveSlots = cqiPracticalPerCW(1,activeSlotNum,:);
    cqiPerfectPerCWActiveSlots = cqiPerfectPerCW(1,activeSlotNum,:);
    widebandSINRPractical = 10*log10(SINRPerCWPractical(1,activeSlotNum,:));
    widebandSINRPerfect = 10*log10(SINRPerCWPerfect(1,activeSlotNum,:));

    if isempty(reshape(cqiPracticalPerCWActiveSlots(:,:,1),1,[]))
        disp('No CQI data to plot, because all CQI values are NaNs.');
        return;
    end

    figure();
    plotWBCQISINR(widebandSINRPerfect,widebandSINRPractical,211,activeSlotNum,'SINR');
    plotWBCQISINR(cqiPerfectPerCWActiveSlots,cqiPracticalPerCWActiveSlots,212,activeSlotNum,'CQI');
end

function plotWBCQISINR(perfectVals,practicalVals,subplotIdx,activeSlotNum,inpText)
%   Plots the wideband SINR and wideband CQI values for each codeword
%   across all specified active slots (1-based) (in which the CQI is
%   reported as other than NaN) for practical and perfect channel
%   estimation cases.

    subplot(subplotIdx)
    plot(perfectVals(:,:,1),'r-o');
    hold on;
    plot(practicalVals(:,:,1),'b-*');
    if ~all(isnan(perfectVals(:,:,2))) % Two codewords
        hold on;
        plot(perfectVals(:,:,2),'r:s');
        hold on;
        plot(practicalVals(:,:,2),'b:d');
        title(['Wideband ' inpText ' Values for Codeword 1&2']);
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.','Codeword 2:Perfect channel est.','Codeword 2:Practical channel est.'});
    else
        title(['Wideband ' inpText ' Values for Codeword 1']);
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.'});
    end
    xlabel('Slots');
    if strcmpi(inpText,'SINR')
        units = ' in dB';
    else
        units = '';
    end
    ylabel(['Wideband ' inpText ' Values' units]);
    xticks(1:size(perfectVals,2));
    xTickLables = num2cell(activeSlotNum(:)-1);
    xticklabels(xTickLables);
    [lowerBound,upperBound] = bounds([practicalVals(:);perfectVals(:)]);
    ylim([lowerBound-1 upperBound+3.5]);
end

function plotSubbandCQIAndSINR(subbandCQIPractical,subbandCQIPerfect,SINRPerCWPractical,SINRPerCWPerfect,activeSlotNum,nslot)
%   Plots the SINR and CQI values for each codeword across all the subbands
%   for practical and perfect channel estimation cases for the given slot
%   number (0-based) among all specified active slots (1-based). The
%   function does not plot the values if CQIMode is 'Wideband' or if the
%   CQI and SINR values are all NaNs in the given slot.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No CQI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    numSubbands = size(subbandCQIPractical,1);
    if numSubbands > 1 && ~any(nslot+1 == activeSlotNum) % Check if the CQI values are reported in the specified slot
        disp(['For the specified slot (' num2str(nslot) '), CQI values are not reported. Please choose another slot number.']);
        return;
    end

    % Plot subband CQI values
    if numSubbands > 1 % Subband mode
        subbandCQIPerCWPractical = subbandCQIPractical(2:end,:,nslot+1);
        subbandCQIPerCWPerfect = subbandCQIPerfect(2:end,:,nslot+1);
        subbandSINRPerCWPractical = 10*log10(SINRPerCWPractical(2:end,:,nslot+1));
        subbandSINRPerCWPerfect = 10*log10(SINRPerCWPerfect(2:end,:,nslot+1));
        figure();
        plotSBCQISINR(subbandSINRPerCWPerfect,subbandSINRPerCWPractical,numSubbands,211,nslot,'SINR')
        plotSBCQISINR(subbandCQIPerCWPerfect,subbandCQIPerCWPractical,numSubbands,212,nslot,'CQI');
    end
end

function plotSBCQISINR(perfectVals,practicalVals,numSubbands,subplotIdx,nslot,inpText)
%   Plots the SINR and CQI values for each codeword across all the subbands
%   for practical and perfect channel estimation cases for the given slot
%   number (0-based). The function does not plot the values if CQIMode is
%   'Wideband' or if the CQI and SINR values are all NaNs in the given
%   slot.

    subplot(subplotIdx)
    plot(perfectVals(:,1),'ro-');
    hold on;
    plot(practicalVals(:,1),'b*-');
    if ~all(isnan(perfectVals(:,2))) % Two codewords
        hold on;
        plot(perfectVals(:,2),'rs:');
        hold on;
        plot(practicalVals(:,2),'bd:');
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.','Codeword 2:Perfect channel est.','Codeword 2:Practical channel est.'});
        title(['Estimated Subband ' inpText ' Values for Codeword 1&2 in Slot ' num2str(nslot)]);
    else % Single codeword
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.'});
        title(['Estimated Subband ' inpText ' Values for Codeword 1 in Slot ' num2str(nslot)]);
    end

    if strcmpi(inpText,'SINR')
        units = ' in dB';
    else
        units = '';
    end
    xlabel('Subbands');
    ylabel(['Subband ' inpText ' Values' units]);
    xticks(1:numSubbands);
    xTickLables = num2cell(1:numSubbands);
    xticklabels(xTickLables);
    xlim([0 numSubbands+1]);
    [lowerBound,upperBound] = bounds([perfectVals(:);practicalVals(:)]);
    ylim([lowerBound-1 upperBound+3.5]);
end

function plotType1PMIAndRI(pmiPracticalPerSlot,pmiPerfectPerSlot,riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,nslot)
%   Plots the RI and PMI i1 indices across all specified active slots
%   (1-based), for practical and perfect channel estimation scenarios. The
%   function also plots the i2 indices of practical and perfect channel
%   estimation scenarios across all specified active slots when the PMI
%   mode is 'Wideband' or plots i2 indices across all the subbands for the
%   specified slot number (0-based) when the PMI mode is 'Subband'.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No PMI and RI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    
    numi1Indices = numel(pmiPracticalPerSlot(activeSlotNum(1)).i1);
    if numi1Indices == 6
        codebookType = 'Type1MultiPanel';
    else
        codebookType = 'Type1SinglePanel';
    end
    
    % Extract wideband PMI indices (i1 values) for slots where NZP-CSI-RS
    % is present
    i1PerfectValsActiveSlots = reshape([pmiPerfectPerSlot(activeSlotNum).i1],numi1Indices,[])';
    i1PracticalValsActiveSlots = reshape([pmiPracticalPerSlot(activeSlotNum).i1],numi1Indices,[])';
    
    if isempty(i1PerfectValsActiveSlots)
        disp('No PMI and RI data to plot, because all PMI and RI values are NaNs.');
        return;
    end
    
    figure;
    % Plot RI
    plotRI(riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,411);
    
    % Extract and plot i11 indices
    i11PerfectVals = i1PerfectValsActiveSlots(:,1);
    i11PracticalVals = i1PracticalValsActiveSlots(:,1);
    plotIxxIndices(i11PerfectVals,i11PracticalVals,activeSlotNum,412,'i11');

    % Extract and plot i12 indices
    i12PerfectVals = i1PerfectValsActiveSlots(:,2);
    i12PracticalVals = i1PracticalValsActiveSlots(:,2);
    plotIxxIndices(i12PerfectVals,i12PracticalVals,activeSlotNum,413,'i12');

    % Extract and plot i13 indices
    i13PerfectVals = i1PerfectValsActiveSlots(:,3);
    i13PracticalVals = i1PracticalValsActiveSlots(:,3);
    plotIxxIndices(i13PerfectVals,i13PracticalVals,activeSlotNum,414,'i13');
    
    % Plot the i141, i142 and i143 indices in type I multi-panel case
    if strcmpi(codebookType,'Type1MultiPanel')
        figure()
        % Extract and plot i141 indices
        i141PerfectVals = i1PerfectValsActiveSlots(:,4);
        i141PracticalVals = i1PracticalValsActiveSlots(:,4);
        plotIxxIndices(i141PerfectVals,i141PracticalVals,activeSlotNum,311,'i141');

        % Extract and plot i142 indices
        i142PerfectVals = i1PerfectValsActiveSlots(:,5);
        i142PracticalVals = i1PracticalValsActiveSlots(:,5);
        plotIxxIndices(i142PerfectVals,i142PracticalVals,activeSlotNum,312,'i142');
    
        % Extract and plot i143 indices
        i143PerfectVals = i1PerfectValsActiveSlots(:,6);
        i143PracticalVals = i1PracticalValsActiveSlots(:,6);
        plotIxxIndices(i143PerfectVals,i143PracticalVals,activeSlotNum,313,'i143');
    end

    % Get the number of subbands
    numSubbands = size(pmiPracticalPerSlot(activeSlotNum(1)).i2,2);
    % Get the number of i2 indices according to codebook type
    numi2Indices = 1;
    if strcmpi(codebookType,'Type1MultiPanel')
        numi2Indices = 3;
    end

    % Get number of active slots
    numActiveSlots = numel(activeSlotNum);
    % Extract i2 values
    i2PerfectVals = reshape([pmiPerfectPerSlot(activeSlotNum).i2],[numSubbands,numi2Indices,numActiveSlots]);     % Of size numActiveSlots-by-numi2Indices-numSubbands
    i2PracticalVals = reshape([pmiPracticalPerSlot(activeSlotNum).i2],[numSubbands,numi2Indices,numActiveSlots]); % Of size numActiveSlots-by-numi2Indices-numSubbands

    % Plot i2 values
    if numSubbands == 1 % Wideband mode
        figure;

        % In type I single-panel case, there is only one i2 index. The
        % first column of i2PerfectVals and i2PracticalVals corresponds to
        % i2 index. In type I multi-panel case, the i2 values are a set of
        % three indices i20, i21, and i22. Each column of i2PerfectVals and
        % i2PracticalVals correspond to i20, i21, and i22 indices. Extract
        % and plot the respective index values
        if strcmpi(codebookType,'Type1SinglePanel')
            % Extract and plot i2 values in each slot
            i2PerfectVals = reshape(i2PerfectVals(:,1,:),[],numActiveSlots).';
            i2PracticalVals = reshape(i2PracticalVals(:,1,:),[],numActiveSlots).';
            plotIxxIndices(i2PerfectVals,i2PracticalVals,activeSlotNum,111,'i2');
        else
            % Extract and plot i20 values in each slot
            i20PerfectVals = reshape(i2PerfectVals(:,1,:),[],numActiveSlots).';
            i20PracticalVals = reshape(i2PracticalVals(:,1,:),[],numActiveSlots).';
            plotIxxIndices(i20PerfectVals,i20PracticalVals,activeSlotNum,311,'i20');

            % Extract and plot i21 values in each slot
            i21PerfectVals = reshape(i2PerfectVals(:,2,:),[],numActiveSlots).';
            i21PracticalVals = reshape(i2PracticalVals(:,2,:),[],numActiveSlots).';
            plotIxxIndices(i21PerfectVals,i21PracticalVals,activeSlotNum,312,'i21');

            % Extract and plot i22 values in each slot
            i22PerfectVals = reshape(i2PerfectVals(:,3,:),[],numActiveSlots).';
            i22PracticalVals = reshape(i2PracticalVals(:,3,:),[],numActiveSlots).';
            plotIxxIndices(i22PerfectVals,i22PracticalVals,activeSlotNum,313,'i22');
        end
    else % Subband mode
        if any(nslot+1 == activeSlotNum)
    
            % In subband mode, plot the PMI i2 indices corresponding to the
            % specified slot number
            figure;

            if strcmpi(codebookType,'Type1SinglePanel')
                % Extract and plot i2 values
                pmiSBi2Perfect = pmiPerfectPerSlot(nslot+1).i2(1,:);
                pmiSBi2Practical = pmiPracticalPerSlot(nslot+1).i2(1,:);
                plotI2xIndices_SB(pmiSBi2Perfect,pmiSBi2Practical,numSubbands,nslot,111,'i2');
            else
                % Extract and plot i20 values
                pmiSBi20Perfect = pmiPerfectPerSlot(nslot+1).i2(1,:);
                pmiSBi20Practical = pmiPracticalPerSlot(nslot+1).i2(1,:);
                plotI2xIndices_SB(pmiSBi20Perfect,pmiSBi20Practical,numSubbands,nslot,311,'i20');
                
                % Extract and plot i21 values
                pmiSBi21Perfect = pmiPerfectPerSlot(nslot+1).i2(2,:);
                pmiSBi21Practical = pmiPracticalPerSlot(nslot+1).i2(2,:);
                plotI2xIndices_SB(pmiSBi21Perfect,pmiSBi21Practical,numSubbands,nslot,312,'i21');
    
                % Extract and plot i22 values
                pmiSBi22Perfect = pmiPerfectPerSlot(nslot+1).i2(3,:);
                pmiSBi22Practical = pmiPracticalPerSlot(nslot+1).i2(3,:);
                plotI2xIndices_SB(pmiSBi22Perfect,pmiSBi22Practical,numSubbands,nslot,313,'i22');
            end
        else
            disp(['For the specified slot (' num2str(nslot) '), PMI i2 indices are not reported. Please choose another slot number.'])
        end
    end
end

function plotType2PMIAndRI(pmiPracticalPerSlot,pmiPerfectPerSlot,riPracticalPerSlot,riPerfectPerSlot,panelDims,numBeams,activeSlotNum,nslot)
%   Plots the grid of beams by highlighting the beams that are used for the
%   precoding matrix generation for the specified slot number (0-based),
%   for practical and perfect channel estimation scenarios.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No PMI and RI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    plotRI(riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,111);
    if ~any(nslot+1 == activeSlotNum)
        disp(['For the specified slot (' num2str(nslot) '), PMI values are not reported. Please choose another slot number.']);
    else
        pmiPractical = pmiPracticalPerSlot(nslot+1);
        pmiPerfect = pmiPerfectPerSlot(nslot+1);
        figure();
        plotType2GridOfBeams(pmiPractical,panelDims,numBeams,'Practical Channel Estimation Scenario',1);
        hold on;
        plotType2GridOfBeams(pmiPerfect,panelDims,numBeams,'Perfect Channel Estimation Scenario',2);
    end
end

function plotRI(riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,subplotIndex)
%   Plots the RI values across all specified active slots (1-based), for
%   practical and perfect channel estimation scenarios.

    % Get number of active slots
    numActiveSlots = numel(activeSlotNum);

    % Extract RI values for slots where NZP-CSI-RS is present
    RIPerfectValsActiveSlots = riPerfectPerSlot(activeSlotNum)';
    RIPracticalValsActiveSlots = riPracticalPerSlot(activeSlotNum)';
    
    if isempty(RIPerfectValsActiveSlots)
        disp('No RI data to plot, because all RI values are NaNs.');
        return;
    end
    
    figure;
    subplot(subplotIndex);
    plot(RIPerfectValsActiveSlots,'r-o');
    hold on;
    plot(RIPracticalValsActiveSlots,'b-*');
    xlabel('Slots')
    ylabel('RI Values');
    xticks(1:numActiveSlots);
    xTickLables = num2cell(activeSlotNum(:)-1);
    xticklabels(xTickLables);
    [~,upperBound] = bounds([RIPerfectValsActiveSlots; RIPracticalValsActiveSlots]);
    xlim([0 numActiveSlots+8]);
    ylim([0 upperBound+1]);
    yticks(0:upperBound+1);
    title('RI Values')
    legend({'Perfect channel est.','Practical channel est.'});
end

function plotType2GridOfBeams(PMISet,panelDims,numBeams,chEstType,subplotNum)
%   Plots the grid of beams by highlighting the beams that are used for the
%   type II codebook based precoding matrix generation.

    N1 = panelDims(1);
    N2 = panelDims(2);    
    % Get the oversampling factors
    O1 = 4;
    O2 = 1 + 3*(N2 ~= 1);

    % Extract q1, q2 values
    qSet = PMISet.i1(1:2);
    q1 = qSet(1)-1;
    q2 = qSet(2)-1;

    % Extract i12 value
    i12 = PMISet.i1(3);
    s = 0;
    % Find the n1, n2 values for all the beams, as defined in TS 38.214
    % Section 5.2.2.2.3
    n1_i12 = zeros(1,numBeams);
    n2_i12 = zeros(1,numBeams);
    for beamIdxI = 0:numBeams-1
        i12minussVal = i12 - s;
        xValues = numBeams-1-beamIdxI:N1*N2-1-beamIdxI;
        CValues = zeros(numel(xValues),1);
        for xIdx = 1:numel(xValues)
            if xValues(xIdx) >= numBeams-beamIdxI
                CValues(xIdx) = nchoosek(xValues(xIdx),numBeams-beamIdxI);
            end
        end
        indices = i12minussVal >= CValues;
        maxIdx = find(indices,1,'last');
        xValue = xValues(maxIdx);
        ei = CValues(maxIdx);
        s = s+ei;
        ni = N1*N2 - 1 - xValue;
        n1_i12(beamIdxI+1) = mod(ni,N1);
        n2_i12(beamIdxI+1) = (ni-n1_i12(beamIdxI+1))/N1;
    end
    m1 = O1*(0:N1-1) + q1;
    m2 = O2*(0:N2-1) + q2;

    % Calculate the indices of orthogonal basis set which corresponds to
    % the reported i12 value
    m1_LBeams = O1*(n1_i12) + q1;
    m2_LBeams = O2*(n2_i12) + q2;
    OrthogonalBeams = [repmat(m1,1,length(m2));reshape(repmat(m2,length(m1),1),1,[])]';

    % Plot the grid of beams
    numCirlcesInRow = N1*O1;
    numCirlcesInCol = N2*O2;
    subplot(2,1,subplotNum);
    circleRadius = 1;
    for colIdx = 0:numCirlcesInCol-1
        for rowIdx = 0:numCirlcesInRow-1
            p = nsidedpoly(1000, 'Center', [2*rowIdx 2*colIdx], 'Radius', circleRadius);
            if any(prod(OrthogonalBeams == [rowIdx colIdx],2))
                h2 = plot(p, 'FaceColor', 'w','EdgeColor','r','LineWidth',2.5);
                hold on;
                if any(prod([m1_LBeams' m2_LBeams'] == [rowIdx colIdx],2))
                    h3 = plot(p, 'FaceColor', 'g','LineStyle','-.');                
                end
            else
                h1 = plot(p, 'FaceColor', 'w');
            end
            hold on;
        end
    end
    rowLength = 2*circleRadius*O1;
    colLength = 2*circleRadius*O2;
    for n2 = 0:N2-1
        for n1 = 0:N1-1
            x1 = -1*circleRadius + rowLength*n1;
            x2 = x1 + rowLength;
            y1 = -1*circleRadius + colLength*n2;
            y2 = y1 + colLength;
            x = [x1, x2, x2, x1, x1];
            y = [y1, y1, y2, y2, y1];
            plot(x, y, 'b-', 'LineWidth', 2);
            hold on;
        end
    end
    
    xlabel('N1O1 beams');
    ylabel('N2O2 beams');
    axis equal;
    set(gca,'xtick',[],'ytick',[]);
    legend([h1 h2 h3],{'Oversampled DFT beams',['Orthogonal basis set with [q1 q2] = [' num2str(q1) ' ' num2str(q2) ']'],'Selected beam group'},'Location','northeast');
    title(['Grid of Beams or DFT Vectors for ' chEstType]);
end

function plotIxxIndices(ixxPerfectVals,ixxPracticalVals,activeSlotNum,subplotInp,pmiIdxType)
%   Plots i11, i12, i13 indices in case of type I single-panel codebooks
%   and plots i141, i142, and i143 in case of type I multi-panel codebooks.

    % Plot ixx values
    subplot(subplotInp)
    plot(ixxPerfectVals,'r-o');
    hold on;
    plot(ixxPracticalVals,'b-*');
    xlabel('Slots')
    ylabel([pmiIdxType ' Indices']);
    % Get number of active slots
    numActiveSlots = numel(activeSlotNum);
    xticks(1:numActiveSlots);
    xTickLables = num2cell(activeSlotNum(:)-1);
    xticklabels(xTickLables);
    [lowerBound,upperBound] = bounds([ixxPerfectVals; ixxPracticalVals]);
    xlim([0 numActiveSlots+8]);
    ylim([lowerBound-2 upperBound+2]);
    title(['PMI: ' pmiIdxType ' Indices']);
    legend({'Perfect channel est.','Practical channel est.'});
end

function plotI2xIndices_SB(pmiSBi2Perfect,pmiSBi2Practical,numSubbands,nslot,subplotInp,pmiIdxType)
%   Plots i2 indices in case of type I single-panel codebooks and plots
%   i20, i21, and i22 in case of type I multi-panel codebooks.

    subplot(subplotInp)
    plot(pmiSBi2Perfect,'r-o');
    hold on;
    plot(pmiSBi2Practical,'b-*');
    title(['PMI: ' pmiIdxType ' Indices for All Subbands in Slot ' num2str(nslot)]);
    xlabel('Subbands')
    ylabel([pmiIdxType ' Indices']);
    xticks(1:numSubbands);
    xticklabels(num2cell(1:numSubbands));
    [lowerBound,upperBound] = bounds([pmiSBi2Perfect pmiSBi2Practical]);
    yticks(lowerBound:upperBound);
    yticklabels(num2cell(lowerBound:upperBound));
    xlim([0 numSubbands+1])
    ylim([lowerBound-1 upperBound+1]);
    legend({'Perfect channel est.','Practical channel est.'});
end

function [lowerCQI, currentCQI, upperCQI] = getAdjacentCQIs(snr, CQI_MCS_Table)
    snrValues = cell2mat({CQI_MCS_Table.SNRdB});
    
    [~, idx] = min(abs(snrValues - snr));
    
    if idx > 1 && idx < numel(snrValues)
        lowerCQI = CQI_MCS_Table(idx - 1);  % 下一个SNR点的CQI
        currentCQI = CQI_MCS_Table(idx);    % 当前SNR点的CQI
        upperCQI = CQI_MCS_Table(idx + 1);  % 上一个SNR点的CQI
    elseif idx == 1
        lowerCQI = CQI_MCS_Table(idx);      % 如果SNR是最小值，选择第一个CQI作为lower
        currentCQI = CQI_MCS_Table(idx);    % 当前SNR点的CQI
        upperCQI = CQI_MCS_Table(idx + 1);  % 选择下一个CQI作为upper
    elseif idx == numel(snrValues)
        lowerCQI = CQI_MCS_Table(idx - 1);  % 选择上一个CQI作为lower
        currentCQI = CQI_MCS_Table(idx);    % 当前SNR点的CQI
        upperCQI = CQI_MCS_Table(idx);      % 如果SNR是最大值，选择最后一个CQI作为upper
    end
end

function p = localWaterfill(alpha, Ptot)
% alpha: [Kx1], alpha_k = N0_k/g_k（数值越小频道越好）；可含 Inf
% Ptot : 标量总功率（相对值）
% 输出 p: [Kx1]，各子带功率（线性），若alpha=Inf则p=0
    K = numel(alpha);
    p = zeros(K,1);
    % 仅对有限alpha参与分配
    finiteIdx = find(isfinite(alpha));
    if isempty(finiteIdx) || Ptot<=0
        return;
    end
    a = alpha(finiteIdx);
    [aSort, order] = sort(a, 'ascend');      % 由好到差
    m = 0;
    while m < numel(aSort)
        m = m + 1;
        mu = (Ptot + sum(aSort(1:m))) / m;   % 理想水位
        if m < numel(aSort) && mu <= aSort(m+1)
            % 水位还没高到能"浇"到下一个更差子带，继续扩m
            continue;
        else
            % 计算功率
            pLoc = max(mu - aSort(1:m), 0);
            p_f = zeros(numel(aSort),1);
            p_f(1:m) = pLoc;
            % 回填到原索引
            tmp = zeros(numel(a),1);
            tmp(order) = p_f;
            p(finiteIdx) = tmp;
            break;
        end
    end
end