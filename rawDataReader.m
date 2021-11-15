
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      (C) Copyright 2018 Texas Instruments, Inc.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions
%  are met:
%
%    Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
%
%    Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the
%    distribution.
%
%    Neither the name of Texas Instruments Incorporated nor the names of
%    its contributors may be used to endorse or promote products derived
%    from this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overview:
% ===========
%   This script reads a JSON file generated from Radar Studio which contains
% data card captured radar data and its cpature and front end parameters.
%   It reads in data from binary file(s) name stored in JSON file and constructs
% raw data in frames, chirps and channels. It then performs range FFT to generate
% radar cube.
%   If file name is provided to export raw ADC data or radar cube, it
% exports data in Matlab mat file.
%   If debug plot option is enabled, it plots raw ADC data and radar cube
% data based on the selection of frame index, chirp index and RX channel.
%
%
% Matlab version:
% =================
%   2017a and later to support JSON files.
%
%
% Assumptions:
% ==============
%     1. Number of ADC bits is 16. (12 and 14 bits are not supported
%     2. Data logging mode
%        - Raw
%     3. LVDS packet format
%        - ADC
%     4. LVDS lanes
%        - xWR12xx/xWR14xx/xWR22xx : 4 lanes
%        - xWR16xx/xWR18xx/xWR68xx : 2 lanes
%     5. Data capture data format
%        - Complex only
%     6. Format
%        - Format 0
%     7. Number of RX antenna
%        - 1, 2, or 4
%     8. Frame configuration
%        - Single frame
%     9. Catpure device
%        - DCA1000 only
%
%
% Syntax:
% =========
%     rawDataReader('JSON setup File name','raw data file name', 'radarCube data file name', 'debug plot')
%          - 'JSON setup File name', JSON setup file
%          - 'raw data file name', file to save the raw ADC data
%          - 'radarCube data file name', file to save radar cube data
%          - 'debug plot', flag to enable plot of raw data and radar cube
%
%
function rawDataReader(setupJsonFileName, rawDataFileName, radarCubeDataFileName, debugPlot)
    close all;

    % Global parameters
    global Params;
    global ui;
    global dataSet;
    global EXIT_KEY_PRESSED
    EXIT_KEY_PRESSED = 0;

    % Read configuration and setup files
    setupJSON = jsondecode(fileread(setupJsonFileName));

    % Read mmwave JSON file
    jsonMmwaveFileName = setupJSON.configUsed;
    mmwaveJSON = jsondecode(fileread(jsonMmwaveFileName));

    % Print parsed current system parameter
    fprintf('mmwave Device:%s\n', setupJSON.mmWaveDevice);

    % Read bin file name
    numBinFiles = length(setupJSON.capturedFiles.files);
    if( numBinFiles < 1)
        error('Bin File is not available');
    end
    Params.numBinFiles = numBinFiles;

    for idx=1:numBinFiles
        binFileName{idx} = strcat(setupJSON.capturedFiles.files(idx).processedFileName);
    end

    % Generate ADC data parameters
    adcDataParams = dp_generateADCDataParams(mmwaveJSON);

    % Generate radar cube parameters
    radarCubeParams = dp_generateRadarCubeParams(mmwaveJSON);

    % Generate RF parameters
    Params.RFParams = dp_generateRFParams(mmwaveJSON, radarCubeParams, adcDataParams);
    Params.NSample = adcDataParams.numAdcSamples;
    Params.NChirp = adcDataParams.numChirpsPerFrame;
    Params.NChan = adcDataParams.numRxChan;
    Params.NTxAnt = radarCubeParams.numTxChan;
    Params.numRangeBins = radarCubeParams.numRangeBins;
    Params.numDopplerBins = radarCubeParams.numDopplerChirps;
    Params.rangeWinType = 0;

    % Validate configuration
    validConf = dp_validateDataCaptureConf(setupJSON, mmwaveJSON);
    if(validConf == false)
        error("Configuraion from JSON file is not Valid");
    end

    % Open raw data from file
    Params.NFrame = 0;
    for idx = 1:numBinFiles
        [Params.fid_rawData(idx), errmsg] = fopen(binFileName{idx}, 'r');
        if(Params.fid_rawData(idx) == -1)
            fprintf("Can not open Bin file %s, - %s\n",binFileName{idx}, errmsg);
            error('Quit with error');
        end

        % Calculate number of Frames in bin File
        try
            Params.NFramePerFile(idx) = dp_getNumberOfFrameFromBinFile(binFileName{idx});
            Params.NFrame = Params.NFrame + Params.NFramePerFile(idx);
        catch
            if(Params.NFramePerFile(idx) == 0)
                error("Not enough data in binary file");
            end
        end
    end

    % Export data
    dp_exportData(rawDataFileName, radarCubeDataFileName);

    % Start example UI and update time domain/range Profile plots
    if(debugPlot)
        % Start up processing display page
        ui.figHandle = initDisplayPage(setupJSON.mmWaveDevice);

        % load and Plot the first frame
        dp_updateFrameData(1);
        ui_updateFramePlot();

        % Wait for UI interactions
        while (~EXIT_KEY_PRESSED)
            pause(0.01);
        end

        close(ui.figHandle);
    end

    %close and delete handles before exiting
    for idx = 1: numBinFiles
        fclose(Params.fid_rawData(idx));
    end
    close all;
end


% ============================================================
% Configuration and Data File Parsing Functions
% ============================================================

%   -----------------------------------------------------------------------
%   Description:    This function loads one frame data and perform range
%   FFT
%   Input:          exportRawDataFile - file name to export raw data
%                   export1DFFTDataFile - file name to export 1D FFT data
%   Output:         mat files
%   -----------------------------------------------------------------------
function dp_exportData(rawDataFileName, radarCubeDataFileName)
    global dataSet
    global Params

    % Prepare data to be saved in mat-file
    if ((~strcmp(rawDataFileName, '')) || (~strcmp(rawDataFileName, '')))
        for frameIdx=1:Params.NFrame
            dp_updateFrameData(frameIdx);
            rawADCData{frameIdx} = dataSet.rawDataUint16;
            radarCubeData{frameIdx} = single(dataSet.radarCubeData);
        end
    end
    % Export raw ADC data
    if (~strcmp(rawDataFileName, ''))
        adcRawData.rfParams = Params.RFParams;
        adcRawData.data = rawADCData;
        adcRawData.dim.numFrames = Params.NFrame;
        adcRawData.dim.numChirpsPerFrame = Params.adcDataParams.numChirpsPerFrame;
        adcRawData.dim.numRxChan = Params.NChan;
        adcRawData.dim.numSamples = Params.NSample;

        % Save params and data to mat file
        save (rawDataFileName, 'adcRawData', '-v7.3');
    end

    % Export rangeFFT data
    if (~strcmp(rawDataFileName, ''))
        radarCubeParams = Params.radarCubeParams;
        radarCube.rfParams = Params.RFParams;
        radarCube.data = radarCubeData;
        radarCube.dim.numFrames = Params.NFrame;
        radarCube.dim.numChirps = radarCubeParams.numTxChan * radarCubeParams.numDopplerChirps;
        radarCube.dim.numRxChan = radarCubeParams.numRxChan;
        radarCube.dim.numRangeBins = radarCubeParams.numRangeBins;
        radarCube.dim.iqSwap = radarCubeParams.iqSwap;

        % Save params and data to mat file
        save (radarCubeDataFileName,'radarCube', '-v7.3');
    end
end
%   -----------------------------------------------------------------------
%   Description:    This function loads one frame data and perform range
%   FFT
%   Input:          frameIdx - frame index
%   Output:         dataSet.rawFrameData(complex)
%                   dataSet.radarCubeData(complex)
%   -----------------------------------------------------------------------
function dp_updateFrameData(frameIdx)
    global Params
    global dataSet

    % Find binFin index
    currFrameIdx = 0;
    fidIdx = 0;
    for idx = 1: Params.numBinFiles
        if frameIdx <= (Params.NFramePerFile(idx) + currFrameIdx)
            fidIdx = idx;
            break;
        else
            currFrameIdx = currFrameIdx + Params.NFramePerFile(idx);
        end
    end

    if(fidIdx <= Params.numBinFiles)
        % Load raw data from bin file
        rawDataComplex = dp_loadOneFrameData(Params.fid_rawData(fidIdx), Params.dataSizeOneFrame, frameIdx - currFrameIdx);

        % Read in raw data in uint16
        dataSet.rawDataUint16 = uint16(rawDataComplex);

        % time domain data y value adjustments
        timeDomainData = rawDataComplex - ( rawDataComplex >=2.^15).* 2.^16;

        % reshape data based on capture configurations
        dp_generateFrameData(timeDomainData);

        % Perform rangeFFT
        dataSet.radarCubeData = processingChain_rangeFFT(Params.rangeWinType);
    end
end


%   -----------------------------------------------------------------------
%   Description:    This function calcultes number of frames of data available
%                   in binary file
%   Input:          binFileName - binary file name
%   Output:         NFrame - number of Frames
%   -----------------------------------------------------------------------
function [NFrame] = dp_getNumberOfFrameFromBinFile(binFileName)
    global Params
    try
        binFile = dir(binFileName);
        fileSize = binFile.bytes;
    catch
        error('Reading Bin file failed');
    end
    NFrame = floor(fileSize/Params.dataSizeOneFrame);
end

%   -----------------------------------------------------------------------
%   Description:    This function load one frame data from binary file
%   Input:          fid_rawData - fid for binary file
%                   dataSizeOneFrame - size of one frame data
%                   frameIdx - frame index
%   Output:         rawData - one frame of raw ADC data
%   -----------------------------------------------------------------------
function [rawData] = dp_loadOneFrameData(fid_rawData, dataSizeOneFrame, frameIdx)
    % find the first byte of the frame
    fseek(fid_rawData, (frameIdx - 1)*dataSizeOneFrame, 'bof');

    try
        % Read in raw data in complex single
        rawData = fread(fid_rawData, dataSizeOneFrame/2, 'uint16=>single');
    catch
        error("error reading binary file");
    end
    if(dataSizeOneFrame ~= length(rawData)*2)
        fprintf("dp_loadOneFrameData, size = %d, expected = %d \n",length(rawData), dataSizeOneFrame);
        error("read data from bin file, have wrong length");
    end
end

%   -----------------------------------------------------------------------
%   Description:    This function validates configuration from JSON files
%   Input:          setupJson - setup JSON configuration structure
%                   mmwaveJSON - mmwave JSON configuration structure
%   Output:         confValid - true if the configuration is valid
%   -----------------------------------------------------------------------
function [confValid] = dp_validateDataCaptureConf(setupJson, mmwaveJSON)
    global Params

    mmWaveDevice = setupJson.mmWaveDevice;

    % Supported platform list
    supportedPlatform = {'awr1642',...
                         'iwr1642',...
                         'awr1243',...
                         'awr1443',...
                         'iwr1443',...
                         'awr1843',...
                         'iwr1843',...
                         'iwr6843',...
                         'awr2243'};

    confValid = true;

    % Validate if the device is supported
    index = find(contains(supportedPlatform,mmWaveDevice));
    if(index == 0)
        fprintf("Platform not supported : %s \n", mmWaveDevice);
        confValid = false;
    end

    % Validate the captureHardware
    if(setupJson.captureHardware ~= 'DCA1000')
        confValid = false;
        fprintf("Capture hardware is not supported : %s \n", setupJson.captureHardware);
    end

    % Validate ADC_ONLY capture
    if (mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevDataPathCfg_t.transferFmtPkt0 ~= '0x1')
        confValid = false;
        fprintf("Capture data format is not supported : %s \n", mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevDataPathCfg_t.transferFmtPkt0);
    end

    % Validate the dataLoggingMode
    if(setupJson.DCA1000Config.dataLoggingMode  ~= 'raw')
        confValid = false;
        fprintf("Capture data logging mode is not supported : %s \n", setupJson.DCA1000Config.dataLoggingMode);
    end

    % Validate the Capture configuration
    Params.numLane = dp_numberOfEnabledChan(sscanf(mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevLaneEnable_t.laneEn, '0x%x'));
    Params.chInterleave = mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevDataFmtCfg_t.chInterleave;
    if ((mmWaveDevice == 'awr1443') | (mmWaveDevice == 'iwr1443') | (mmWaveDevice == 'awr1243') | (mmWaveDevice == 'awr2243'))
        if(Params.numLane ~= 4)
            fprintf(" %d LVDS Lane is not supported for device : %s ", Params.numLane, mmWaveDevice);
            confValid = false;
        end

        if(Params.chInterleave ~= 0)
            fprintf(" Interleave mode %d is not supported for device : %s ", Params.chInterleave, mmWaveDevice);
            confValid = false;
        end
    else
        if(Params.numLane ~= 2)
            fprintf(" %d LVDS Lane is not supported for device : %s ", Params.numLane, mmWaveDevice);
            confValid = false;
        end
        if(Params.chInterleave ~= 1)
            fprintf(" Interleave mode %d is not supported for device : %s ", Params.chInterleave, mmWaveDevice);
            confValid = false;
        end
    end
end


%  -----------------------------------------------------------------------
%  Description:    This function reshape raw binary data based on capture
%                  configuration, generates data in cell of
%                  [number of chirps, number of RX channels, number of ADC samples]
%  Input:          rawData - raw ADC data
%  Output:         frameData - reshaped ADC data
%  -----------------------------------------------------------------------
function [frameData] = dp_generateFrameData(rawData)
    global Params
    global dataSet

    if(Params.numLane == 2)
        % Convert 2 lane LVDS data to one matrix
        frameData = dp_reshape2LaneLVDS(rawData);
    elseif (Params.numLane == 4)
        % Convert 4 lane LVDS data to one matrix
        frameData = dp_reshape4LaneLVDS(rawData);
    else
        fprintf("%d LVDS lane is not supported ", Params.numLane);
    end

    % checking iqSwap setting
    if(Params.adcDataParams.iqSwap == 1)
        % Data is in ReIm format, convert to ImRe format to be used in radarCube
        frameData(:,[1,2]) = frameData(:,[2,1]);
    end

    % Convert data to complex: column 1 - Imag, 2 - Real
    frameCplx = frameData(:,1) + 1i*frameData(:,2);

    % initialize frameComplex
    frameComplex = single(zeros(Params.NChirp, Params.NChan, Params.NSample));

    % Change Interleave data to non-interleave
    if(Params.chInterleave == 1)
        % non-interleave data
        temp = reshape(frameCplx, [Params.NSample * Params.NChan, Params.NChirp]).';
        for chirp=1:Params.NChirp
            frameComplex(chirp,:,:) = reshape(temp(chirp,:), [Params.NSample, Params.NChan]).';
        end
    else
        % interleave data
        temp = reshape(frameCplx, [Params.NSample * Params.NChan, Params.NChirp]).';
        for chirp=1:Params.NChirp
            frameComplex(chirp,:,:) = reshape(temp(chirp,:), [Params.NChan, Params.NSample]);
        end
    end

    % Save raw data
    dataSet.rawFrameData = frameComplex;
end

%  -----------------------------------------------------------------------
%  Description:    This function counts number of enabled channels from
%                  channel Mask.
%  Input:          chanMask
%  Output:         Number of channels
%  -----------------------------------------------------------------------
function [count] = dp_numberOfEnabledChan(chanMask)

    MAX_RXCHAN = 4;
    count = 0;
    for chan= 0:MAX_RXCHAN - 1
        bitVal = pow2(chan);
        if (bitand(chanMask,bitVal) == (bitVal))
            count = count + 1;
            chanMask = chanMask-bitVal;
            if(chanMask == 0)
                break;
            end
        end
    end
end

%  -----------------------------------------------------------------------
%  Description:    This function generates ADC raw data Parameters
%  Input:          mmwaveJSON
%  Output:         adcDataParams
%  -----------------------------------------------------------------------
function [adcDataParams] = dp_generateADCDataParams(mmwaveJSON)
    global Params
    frameCfg = mmwaveJSON.mmWaveDevices.rfConfig.rlFrameCfg_t;

    adcDataParams.dataFmt = mmwaveJSON.mmWaveDevices.rfConfig.rlAdcOutCfg_t.fmt.b2AdcOutFmt;
    adcDataParams.iqSwap = mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevDataFmtCfg_t.iqSwapSel;
    adcDataParams.chanInterleave = mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevDataFmtCfg_t.chInterleave;
    adcDataParams.numChirpsPerFrame = frameCfg.numLoops  * (frameCfg.chirpEndIdx - frameCfg.chirpStartIdx + 1);
    adcDataParams.adcBits = mmwaveJSON.mmWaveDevices.rfConfig.rlAdcOutCfg_t.fmt.b2AdcBits;
    rxChanMask = sscanf(mmwaveJSON.mmWaveDevices.rfConfig.rlChanCfg_t.rxChannelEn, '0x%x');

    adcDataParams.numRxChan = dp_numberOfEnabledChan(rxChanMask);
    adcDataParams.numAdcSamples = mmwaveJSON.mmWaveDevices.rfConfig.rlProfiles.rlProfileCfg_t.numAdcSamples;

    dp_printADCDataParams(adcDataParams);

    % Calculate ADC data size
    if(adcDataParams.adcBits == 2)
        if (adcDataParams.dataFmt == 0)
            % real data, one sample is 16bits=2bytes
            gAdcOneSampleSize = 2;
        elseif ((adcDataParams.dataFmt == 1) || (adcDataParams.dataFmt == 2))
            % complex data, one sample is 32bits = 4 bytes
            gAdcOneSampleSize = 4; %2 bytes
        else
            fprintf('Error: unsupported ADC dataFmt');
        end
    else
        fprintf('Error: unsupported ADC bits (%d)', adcDataParams.adcBits);
    end

    dataSizeOneChirp = gAdcOneSampleSize * adcDataParams.numAdcSamples * adcDataParams.numRxChan;
    Params.dataSizeOneFrame = dataSizeOneChirp * adcDataParams.numChirpsPerFrame;
    Params.dataSizeOneChirp = dataSizeOneChirp;

    Params.adcDataParams = adcDataParams;
end

%  -----------------------------------------------------------------------
%  Description:    This function prints ADC raw data Parameters
%  Input:          adcDataParams
%  Output:         None
%  -----------------------------------------------------------------------
function [] = dp_printADCDataParams(adcDataParams)
    fprintf('Input ADC data parameters:\n');
    fprintf('    dataFmt:%d\n',adcDataParams.dataFmt);
    fprintf('    iqSwap:%d\n',adcDataParams.iqSwap);
    fprintf('    chanInterleave:%d\n',adcDataParams.chanInterleave);
    fprintf('    numChirpsPerFrame:%d\n',adcDataParams.numChirpsPerFrame);
    fprintf('    adcBits:%d\n',adcDataParams.adcBits);
    fprintf('    numRxChan:%d\n',adcDataParams.numRxChan);
    fprintf('    numAdcSamples:%d\n',adcDataParams.numAdcSamples);
end

%  -----------------------------------------------------------------------
%  Description:    This function generates radar cube data Matrix Parameters
%  Input:          mmwaveJSON
%  Output:         radarCubeParams
%  -----------------------------------------------------------------------
function [radarCubeParams] = dp_generateRadarCubeParams(mmwaveJSON)
    global Params

    frameCfg = mmwaveJSON.mmWaveDevices.rfConfig.rlFrameCfg_t;

    radarCubeParams.iqSwap = mmwaveJSON.mmWaveDevices.rawDataCaptureConfig.rlDevDataFmtCfg_t.iqSwapSel;
    rxChanMask = sscanf(mmwaveJSON.mmWaveDevices.rfConfig.rlChanCfg_t.rxChannelEn, '0x%x');
    radarCubeParams.numRxChan = dp_numberOfEnabledChan(rxChanMask);
    radarCubeParams.numTxChan = frameCfg.chirpEndIdx - frameCfg.chirpStartIdx + 1;

    radarCubeParams.numRangeBins = pow2(nextpow2(mmwaveJSON.mmWaveDevices.rfConfig.rlProfiles.rlProfileCfg_t.numAdcSamples));
    radarCubeParams.numDopplerChirps = mmwaveJSON.mmWaveDevices.rfConfig.rlFrameCfg_t.numLoops;

    % 1D Range FFT output : cmplx16ImRe_t x[numChirps][numRX][numRangeBins]
    radarCubeParams.radarCubeFmt = 1; %RADAR_CUBE_FORMAT_1;

    dp_printRadarCubeParams(radarCubeParams);
    Params.radarCubeParams = radarCubeParams;
end

%  -----------------------------------------------------------------------
%  Description:    This function prints radar cube data Matrix Parameters
%  Input:          mmwaveJSON
%  Output:         radarCubeParams
%  -----------------------------------------------------------------------
function [] = dp_printRadarCubeParams(radarCubeParams)
    fprintf('Radarcube parameters:\n');
    fprintf('    iqSwap:%d\n',radarCubeParams.iqSwap);
    fprintf('    radarCubeFmt:%d\n',radarCubeParams.radarCubeFmt);
    fprintf('    numDopplerChirps:%d\n',radarCubeParams.numDopplerChirps);
    fprintf('    numRxChan:%d\n',radarCubeParams.numRxChan);
    fprintf('    numTxChan:%d\n',radarCubeParams.numTxChan);
    fprintf('    numRangeBins:%d\n',radarCubeParams.numRangeBins);
end

%  -----------------------------------------------------------------------
%  Description:    This function generates mmWave Sensor RF parameters
%  Input:          mmwaveJSON, radarCubeParams, adcDataParams
%  Output:         None
%  -----------------------------------------------------------------------
function [RFParams] = dp_generateRFParams(mmwaveJSON, radarCubeParams, adcDataParams)

    C = 3e8;
    profileCfg = mmwaveJSON.mmWaveDevices.rfConfig.rlProfiles.rlProfileCfg_t;


    RFParams.startFreq = profileCfg.startFreqConst_GHz;

    % Slope const (MHz/usec)
    RFParams.freqSlope = profileCfg.freqSlopeConst_MHz_usec;

    % ADC sampling rate in Msps
    RFParams.sampleRate = profileCfg.digOutSampleRate / 1e3;

    % Generate radarCube parameters
    RFParams.numRangeBins = pow2(nextpow2(adcDataParams.numAdcSamples));
    RFParams.numDopplerBins = radarCubeParams.numDopplerChirps;
    RFParams.bandwidth = abs(RFParams.freqSlope * profileCfg.numAdcSamples / profileCfg.digOutSampleRate);

    RFParams.rangeResolutionsInMeters = C * RFParams.sampleRate / (2 * RFParams.freqSlope * RFParams.numRangeBins * 1e6);
    RFParams.dopplerResolutionMps =  C  / (2*RFParams.startFreq * 1e9 *...
                                        (profileCfg.idleTimeConst_usec + profileCfg.rampEndTime_usec  ) *...
                                        1e-6 * radarCubeParams.numDopplerChirps * radarCubeParams.numTxChan);
    RFParams.framePeriodicity = mmwaveJSON.mmWaveDevices.rfConfig.rlFrameCfg_t.framePeriodicity_msec;

end

%  -----------------------------------------------------------------------
%  Description:    This function reshape raw data for 2 lane LVDS capture
%  Input:          rawData - raw ADC data from binary file
%  Output:         frameData
%  -----------------------------------------------------------------------
function [frameData] = dp_reshape2LaneLVDS(rawData)
    % Convert 2 lane LVDS data to one matrix
    rawData4 = reshape(rawData, [4, length(rawData)/4]);
    rawDataI = reshape(rawData4(1:2,:), [], 1);
    rawDataQ = reshape(rawData4(3:4,:), [], 1);

    frameData = [rawDataI, rawDataQ];
end

%  -----------------------------------------------------------------------
%  Description:    This function reshape raw data for 4 lane LVDS capture
%  Input:          rawData - raw ADC data from binary file
%  Output:         frameData
%  -----------------------------------------------------------------------
function [frameData] = dp_reshape4LaneLVDS(rawData)
    % Convert 4 lane LVDS data to one matrix
    rawData8 = reshape(rawData, [8, length(rawData)/8]);
    rawDataI = reshape(rawData8(1:4,:), [], 1);
    rawDataQ = reshape(rawData8(5:8,:), [], 1);

    frameData= [rawDataI, rawDataQ];
end


% ============================================================
% Processing Chain Functions
% ============================================================

%  -----------------------------------------------------------------------
%  Description:    This function is part of processing Chain to perform
%                  range FFT operation.
%  Input:          frameData in complex(time domain)
%  Output:         radarCube
%  -----------------------------------------------------------------------
function [radarCubeData] = processingChain_rangeFFT(rangeWinType)
    global Params
    global dataSet
    global ui

    NChirp = Params.NChirp;
    NChan = Params.NChan;
    NRangeBin = Params.numRangeBins;

    % generate windowing
    switch rangeWinType
        case 1 %hann
            win = hann(Params.NSample);
        case 2 %blackman
            win = blackman(Params.NSample);
        otherwise
            win = rectwin(Params.NSample);
    end

    radarCubeData = single(zeros(NChirp,NChan,NRangeBin));
    for chirpIdx=1:Params.NChirp
        for chIdx = 1: Params.NChan
            frameData(1,:) = dataSet.rawFrameData(chirpIdx,chIdx,:);
            frameData = fft(frameData .* win', NRangeBin);
            radarCubeData(chirpIdx,chIdx,:) = frameData(1,:);
        end
    end
end

% ============================================================
% ============================================================
% User Interface Functions
% ============================================================
% ============================================================
%  -----------------------------------------------------------------------
%  Description:    Callback function for linear plot checkbox
%  -----------------------------------------------------------------------
function ui_checkbox_LinRangeProf_Callback(hObject, ~, ~)
    global ui

    if (get(hObject,'Value') == get(hObject,'Max'))
      ui.rangePanel.linearPlot = 1;
    else
      ui.rangePanel.linearPlot = 0;
    end

    % Update plot
    ui_RangeProfile_updatePlot();
end

%  -----------------------------------------------------------------------
%  Description:    Callback function for key press for UI, if 'q' is
%                  pressed, exit from the script
%  -----------------------------------------------------------------------
function myKeyPressFcn(~, event)
    global EXIT_KEY_PRESSED

    if lower(event.Key) == 'q'
        EXIT_KEY_PRESSED  = 1;
    end
end

%  -----------------------------------------------------------------------
%  Description:    UI initialization
%  -----------------------------------------------------------------------
function [figHandle] = initDisplayPage(platformType)
    global ui

    %Setup the main figure
    figHandle = figure(1);
    clf(figHandle);

    set(figHandle,'Name',strcat('Texas Instruments mmWave device - ', upper(platformType), ' Post Processing Tool'),'NumberTitle','off');

    warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
    jframe=get(figHandle,'javaframe');
    jIcon=javax.swing.ImageIcon('texas_instruments.gif');
    jframe.setFigureIcon(jIcon);
    set(figHandle, 'Color', [0.8 0.8 0.8]);
    set(figHandle, 'KeyPressFcn', @myKeyPressFcn)

    % Preset panel positions in percentage
    timeDomainPanelPos = [0.0,0.4,0.5,0.6];
    rangePlotPanelPos = [0.5,0.4,0.5,0.6];

    % Adding panels
    ui.hTimedomainPanel = uipanel('Title','Time Domain Plot', 'FontSize',12, 'Units','norm', 'Pos',timeDomainPanelPos);
    ui.hRangePlotPanel = uipanel('Title','Range Profile Plot', 'FontSize',12, 'Units','norm', 'Pos',rangePlotPanelPos);

    % Adding Control panels
    ui.hControlPanel = uipanel('Title','Post Processing Control', 'FontSize',12, 'Units','norm', 'Pos',[0,0,.25,0.4]);
    % Adding Configuration panels
    ui.hChirpParamPanel = uipanel('Title','mmWave RF Configuration', 'FontSize',12, 'Units','norm', 'Pos',[0.25,0,.25,0.4]);
    ui.hAdcParamPanel = uipanel('Title','mmWave ADC data Configuration', 'FontSize',12, 'Units','norm', 'Pos',[0.5,0,.25,0.4]);
    ui.hRadarCubeParamPanel = uipanel('Title','mmWave Radar Cube Configuration', 'FontSize',12, 'Units','norm', 'Pos',[0.75,0,.25,0.4]);

    % Adding Configuration Table
    ui_displayConfParams();

    % Setup Control panel
    ui_setupControlPanel();

    % Setup Plot panels
    ui_setupTimeDomainPanel();
    ui_setupRangeProfilePanel();
    %ui_setupDopplerPlotPanel();

    pause(0.00001);
    set(jframe,'Maximized',1);
    pause(0.00001);
end

%  -----------------------------------------------------------------------
%  Description:    This function update current frame time domain and range profile plot
%                  with current chirp and channel index settings
%  -----------------------------------------------------------------------
function ui_updateFramePlot()
    global Params
    global dataSet
    global ui

    % Plot Time domain and range profile
    chirpFFTData = zeros(Params.NChan, Params.numRangeBins);
    for ch= 1: Params.NChan
        chirpFFTData(ch,:) = dataSet.radarCubeData(ui.rangePanel.chirpIdx, ch, :);

        % Plot 1D range Profile
        ui_plot1DRangeProfile(ch, ui.rangePanel.chirpIdx, ui.rangePanel.linearPlot);

        % Plot time domain signal
        ui_plotTimeDomainData(ui.timeDomainPanel.chanIdx, ui.timeDomainPanel.chirpIdx);
    end

end

%  -----------------------------------------------------------------------
%  Description:    Callback function for frame Edit window
%  -----------------------------------------------------------------------
function ui_ControlPanel_frameEdit_Callback(~, ~)
    global Params
    global ui

    if ((str2double(ui.ctrlPanel.frameNumber.String) <= Params.NFrame) && ...
         (str2double(ui.ctrlPanel.frameNumber.String) >= 1))
        % Keep the original setting
    else
        set(ui.ctrlPanel.frameNumber, 'String',num2str(ui.ctrlPanel.frameIdx));
    end
    ui.ctrlPanel.frameIdx = str2double(ui.ctrlPanel.frameNumber.String);

    % Update data and Plot
    dp_updateFrameData(ui.ctrlPanel.frameIdx);
    ui_updateFramePlot();
end


%  -----------------------------------------------------------------------
%  Description:    Function to setup control panel
%  -----------------------------------------------------------------------
function [] = ui_setupControlPanel()
    global Params
    global ui

    % Default frame index
    ui.ctrlPanel.frameIdx = 1;

    % Create frameSliderText
    uiCtrlPosition = [0.05 0.9 .35 .05];
    uicontrol(ui.hControlPanel, 'Style', 'text', 'String', ['Frame Index(1-' num2str(Params.NFrame) ')'],...
                'FontSize', 9,...
                'FontWeight','bold',...
                'Units','normalized',...
                'Position', uiCtrlPosition);

    uiCtrlPosition = [0.4 0.9 .15 .05];
    ui.ctrlPanel.frameNumber = uicontrol(ui.hControlPanel, 'Style', 'edit', 'Value', 1,...
                'Units','normalized',...
                'Position', uiCtrlPosition,...
                'Max', Params.NFrame,...
                'Min', 1,...
                'Callback', @ui_ControlPanel_frameEdit_Callback);

    set(ui.ctrlPanel.frameNumber, 'String',num2str(ui.ctrlPanel.frameIdx));

    % Adding Windowing radio button
    uiCtrlPosition = [.05 0.5  .4 .3];
    ui.ctrlPanel.windowSel = uibuttongroup(ui.hControlPanel,'Visible','on',...
                  'Title', 'Range Windowing',...
                  'FontSize', 9,...
                  'Position',uiCtrlPosition,...
                  'SelectionChangedFcn', @ui_radioButton_winType_Callback);

    uiCtrlPosition = [.05 0.1  0.9 .25];
    uicontrol(ui.ctrlPanel.windowSel, 'Style', 'radioButton', ...
        'String','Rect',...
        'FontSize', 9,...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position', uiCtrlPosition,...
        'HandleVisibility', 'off');
    uiCtrlPosition = [.05 0.7  0.9 .25];
    uicontrol(ui.ctrlPanel.windowSel, 'Style', 'radioButton', 'String', 'Hanning',...
        'FontSize', 9,...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position', uiCtrlPosition,...
        'HandleVisibility', 'off');

    uiCtrlPosition = [.05 0.4  0.9 .25];
    uicontrol(ui.ctrlPanel.windowSel, 'Style', 'radioButton', ...
        'String','Blackman',...
        'FontSize', 9,...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position', uiCtrlPosition,...
        'HandleVisibility', 'off');

    ui.ctrlPanel.rangeWinType = 0;
end

function ui_radioButton_winType_Callback(~,~)
    global ui
    global Params

    switch ui.ctrlPanel.windowSel.SelectedObject.String
        case 'Hanning'
            Params.rangeWinType = 1;
        case 'Blackman'
            Params.rangeWinType = 2;
        case 'Rect'
            Params.rangeWinType = 0;
        otherwise
            Params.rangeWinType = 0;
    end

    % Update plot
    dp_updateFrameData(ui.ctrlPanel.frameIdx);
    ui_RangeProfile_updatePlot();
end
%  -----------------------------------------------------------------------
%  Description:    Setup Time domain plot panel
%  -----------------------------------------------------------------------
function [] = ui_setupTimeDomainPanel()
    global Params
    global ui

    % Create UIAxes
    ui.timeDomainPanel.hTabGroup = uitabgroup(ui.hTimedomainPanel, 'Position', [0.1 0.20 0.8 0.8]);
    for chan = 1: Params.NChan
        ui.tsPanel(chan).htab = uitab(ui.timeDomainPanel.hTabGroup, 'Title', ['Chan' num2str(chan)]);

        ui.tsPanel(chan).hax = axes('Parent', ui.tsPanel(chan).htab );
        axis(ui.tsPanel(chan).hax, [1 Params.NSample -50000 50000]);

        title(ui.tsPanel(chan).hax , 'Time Domain Plot')
        set(ui.tsPanel(chan).hax ,'Color',[0 0 0.5]);
        ui.tsPanel(chan).hplot  = plot(ui.tsPanel(chan).hax , 0,0,'g.', 'Marker', '.','MarkerSize',5);

        ui.tsPanel(chan).hplot = plot(ui.tsPanel(chan).hax,...
                                       (1:Params.NSample),...
                                       zeros(length((1:Params.NSample)),2),...
                                       '-');

        hline = findobj(ui.tsPanel(chan).hplot, 'type', 'line');
        set(hline(1),'LineStyle','-', 'color',[0 0 1]);
        set(hline(2),'LineStyle',':', 'color',[1 0 0]);
        hold on;
        xlabel('Sample index');
        ylabel('ADC time domain output');
        title('Time Domain Output');
        grid on;
    end
    ui.timeDomainPanel.chanIdx = 1;
    ui.timeDomainPanel.hTabGroup.SelectedTab = ui.tsPanel(1).htab;
    set(ui.timeDomainPanel.hTabGroup, 'SelectionChangedFcn', @ui_TimeDoman_tabGroup_Callback);

    % Create ChirpSlider
    uiCtrlPosition = [0.1 0.1 .1 .04];
    uicontrol(ui.hTimedomainPanel, 'Style', 'text', 'String', 'Chirp Index',...
                'FontSize', 9,...
                'FontWeight','bold',...
                'Units','normalized',...
                'Position', uiCtrlPosition);

    uiCtrlPosition = [0.2 0.1 .05 .04];
    ui.timeDomainPanel.ChirpNumber = uicontrol(ui.hTimedomainPanel, 'Style', 'edit', 'Value', 1,...
                'Units','normalized',...
                'Position', uiCtrlPosition,...
                'Max', Params.NChirp,...
                'Min', 1,...
                'Callback', @ui_Timedomain_chirpEdit_Callback);

    uiCtrlPosition = [0.25 0.1  .15 .04];
    ui.timeDomainPanel.ChirpSlider = uicontrol(ui.hTimedomainPanel, 'Style', 'slider', 'String', 'Chirp Index',...
        'Units','normalized',...
        'Position', uiCtrlPosition,...
        'Value', 1,...
        'Callback', @ui_Timedomain_chirpIdx_Callback);
    ui.timeDomainPanel.ChirpSlider.Max = 1;
    ui.timeDomainPanel.ChirpSlider.Value = 0;
    ui.timeDomainPanel.chirpIdx = 1;

    set(ui.timeDomainPanel.ChirpNumber, 'String',num2str(ui.timeDomainPanel.chirpIdx));
end


%  -----------------------------------------------------------------------
%  Description:    Callback function to set time domain plot channel index from
%                  channels tab group
%  -----------------------------------------------------------------------
function ui_TimeDoman_tabGroup_Callback(~, ~)
    global ui

    % Set channel index from active tab title
    switch ui.timeDomainPanel.hTabGroup.SelectedTab.Title
        case 'Chan1'
            ui.timeDomainPanel.chanIdx = 1;
        case 'Chan2'
            ui.timeDomainPanel.chanIdx = 2;
        case 'Chan3'
            ui.timeDomainPanel.chanIdx = 3;
        case 'Chan4'
            ui.timeDomainPanel.chanIdx = 4;
        otherwise
    end

    % Upate time domain plot with updated channel index
    ui_plotTimeDomainData(ui.timeDomainPanel.chanIdx, ui.timeDomainPanel.chirpIdx);
end


%  -----------------------------------------------------------------------
%  Description:    Callback function to set time domain plot chirp index from
%                  chirp index edit window
%  -----------------------------------------------------------------------
function ui_Timedomain_chirpEdit_Callback(~, ~)
    global Params
    global ui

    if ((str2double(ui.timeDomainPanel.ChirpNumber.String) <= Params.NChirp) && ...
         (str2double(ui.timeDomainPanel.ChirpNumber.String) >= 1))
    else
        set(ui.timeDomainPanel.ChirpNumber, 'String',num2str(ui.timeDomainPanel.chirpIdx));
    end

    % Update chirp index
    ui.timeDomainPanel.chirpIdx = str2double(ui.timeDomainPanel.ChirpNumber.String);
    set(ui.timeDomainPanel.ChirpSlider, 'Value', ui.timeDomainPanel.chirpIdx/Params.NChirp);

    % Upate time domain plot
    ui_plotTimeDomainData(ui.timeDomainPanel.chanIdx, ui.timeDomainPanel.chirpIdx);
end

%  -----------------------------------------------------------------------
%  Description:    Callback function to set time domain plot chirp index from
%                  chirp index slider
%  -----------------------------------------------------------------------
function ui_Timedomain_chirpIdx_Callback(~, ~)
    global Params
    global ui

    ui.timeDomainPanel.chirpIdx = round(ui.timeDomainPanel.ChirpSlider.Value * Params.NChirp, 0);
    if ui.timeDomainPanel.chirpIdx == 0
        ui.timeDomainPanel.chirpIdx = 1;
    end

    % Update chirp index
    set(ui.timeDomainPanel.ChirpNumber, 'String',num2str(ui.timeDomainPanel.chirpIdx));

    % Update time domain plot
    ui_plotTimeDomainData(ui.timeDomainPanel.chanIdx, ui.timeDomainPanel.chirpIdx);
end


%  -----------------------------------------------------------------------
%  Description:    Setup range profile Panel
%  -----------------------------------------------------------------------
function [] = ui_setupRangeProfilePanel()
    global Params
    global ui

    % Create UIAxes
    ui.rangePanel.hTabGroup = uitabgroup(ui.hRangePlotPanel, 'Position', [0.1 0.20 0.8 0.8]);
    for chan = 1:Params.NChan

        ui.rpPanel(chan).htab = uitab(ui.rangePanel.hTabGroup, 'Title', ['Chan' num2str(chan)]);
        ui.rpPanel(chan).hax = axes('Parent', ui.rpPanel(chan).htab );
        axis(ui.rpPanel(chan).hax, [0 Params.numRangeBins * Params.RFParams.rangeResolutionsInMeters 0 Params.numRangeBins]);
        set(ui.rpPanel(chan).hax ,'Color',[0 0 0.5]);
        %ui.rpPanel(chan).hplot  = plot(ui.rpPanel(chan).hax , 0,0,'g.');
        rangeBin = linspace(0, Params.numRangeBins * Params.RFParams.rangeResolutionsInMeters, Params.numRangeBins);
        ydata = zeros(Params.numRangeBins, 1).';
        ui.rpPanel(chan).hplot = plot(ui.rpPanel(chan).hax,...
                                       rangeBin,...
                                       ydata,...
                                       '-');
        hline = findobj(ui.rpPanel(chan).hplot, 'type', 'line');
        set(hline(1),'LineStyle','-', 'color',[0 0 1]);
        hold on;
        xlabel('Range (m)');
        ylabel('Range FFT output (dB)');
        title('Range Profile');
        grid on;
    end
    ui.rangePanel.chanIdx = 1;
    ui.rangePanel.hTabGroup.SelectedTab = ui.rpPanel(1).htab;
    set(ui.rangePanel.hTabGroup, 'SelectionChangedFcn', @ui_RangeProfile_tabGroup_Callback);

    % Adding Linear control radio button
    uiCtrlPosition = [.5 0.1  .15 .03];
    uicontrol(ui.hRangePlotPanel, 'Style', 'checkbox', 'String', 'Linear scale',...
        'FontSize', 9,...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position', uiCtrlPosition,...
        'Value', 0,...
        'Callback', @ui_checkbox_LinRangeProf_Callback);

    ui.rangePanel.linearPlot = 0;

    % Create ChirpSlider
    uiCtrlPosition = [0.1 0.1 .1 .04];
    uicontrol(ui.hRangePlotPanel, 'Style', 'text', 'String', 'Chirp Index',...
                'FontSize', 9,...
                'FontWeight','bold',...
                'Units','normalized',...
                'Position', uiCtrlPosition);

    uiCtrlPosition = [0.2 0.1 .05 .04];
    ui.rangePanel.ChirpNumber = uicontrol(ui.hRangePlotPanel, 'Style', 'edit', 'Value', 1,...
                'Units','normalized',...
                'Position', uiCtrlPosition,...
                'Max', Params.NChirp,...
                'Min', 1,...
                'Callback', @ui_RangeProfile_chirpEdit_Callback);

    uiCtrlPosition = [0.25 0.1  .15 .04];
    ui.rangePanel.ChirpSlider = uicontrol(ui.hRangePlotPanel, 'Style', 'slider', 'String', 'Chirp Index',...
        'Units','normalized',...
        'Position', uiCtrlPosition,...
        'Value', 1,...
        'Callback', @ui_RangeProfile_chirpIdx_Callback);
    ui.rangePanel.ChirpSlider.Max = 1;
    ui.rangePanel.ChirpSlider.Value = 0;
    ui.rangePanel.chirpIdx = 1;

    set(ui.rangePanel.ChirpNumber, 'String',num2str(ui.rangePanel.chirpIdx));
end

%  -----------------------------------------------------------------------
%  Description:    Callback function to update chan index from range
%                  profile tab group
%  -----------------------------------------------------------------------
function ui_RangeProfile_chirpEdit_Callback(~, ~)
    global Params
    global ui

    if ((str2double(ui.rangePanel.ChirpNumber.String) <= Params.NChirp) && ...
         (str2double(ui.rangePanel.ChirpNumber.String) >= 1))
    else
        set(ui.rangePanel.ChirpNumber, 'String',num2str(ui.rangePanel.chirpIdx));
    end

    ui.rangePanel.chirpIdx = str2double(ui.rangePanel.ChirpNumber.String);
    set(ui.rangePanel.ChirpSlider, 'Value', ui.rangePanel.chirpIdx/Params.NChirp);
    ui_RangeProfile_updatePlot();
end

%  -----------------------------------------------------------------------
%  Description:    callback function to update chan index from debug GUI
%  -----------------------------------------------------------------------
function ui_RangeProfile_tabGroup_Callback(~, ~)
    global ui

    switch ui.rangePanel.hTabGroup.SelectedTab.Title
        case 'Chan1'
            ui.rangePanel.chanIdx = 1;
        case 'Chan2'
            ui.rangePanel.chanIdx = 2;
        case 'Chan3'
            ui.rangePanel.chanIdx = 3;
        case 'Chan4'
            ui.rangePanel.chanIdx = 4;
        otherwise
            return
    end

    ui_RangeProfile_updatePlot();
end

%  -----------------------------------------------------------------------
%  Description:    Callback function to update chirp index from debug UI chirp
%                  text box
%  -----------------------------------------------------------------------
function ui_RangeProfile_chirpIdx_Callback(~, ~)
    global Params
    global ui

    ui.rangePanel.chirpIdx = round(ui.rangePanel.ChirpSlider.Value * Params.NChirp, 0);
    if ui.rangePanel.chirpIdx == 0
        ui.rangePanel.chirpIdx = 1;
    end

    % Update chirp index from UI
    set(ui.rangePanel.ChirpNumber, 'String',num2str(ui.rangePanel.chirpIdx));

    % Update range profile plot
    ui_RangeProfile_updatePlot();
end

%  -----------------------------------------------------------------------
%  Description:    Update range Profile plot with updated chanIdx and chirpIdx
%                   from debug GUI
%  -----------------------------------------------------------------------
function ui_RangeProfile_updatePlot()
    global ui

    ui_plot1DRangeProfile(ui.rangePanel.chanIdx, ui.rangePanel.chirpIdx, ui.rangePanel.linearPlot);
end


%  ------------------------------------------------------------------------
%  Description:    Plot time domain data Imginary and real data.
%                   Two plots for each channel, one for Imginary data, one for
%                   real data.
%  Input:      chanIdx - RX channel index
%              chirpIdx - one chirp of range profile data for all
%                         RX channels.
%  ------------------------------------------------------------------------
function [] = ui_plotTimeDomainData(chanIdx, chirpIdx)
    global ui
    global dataSet

    currChDataQ = real(dataSet.rawFrameData(chirpIdx,chanIdx,:));
    currChDataI = imag(dataSet.rawFrameData(chirpIdx,chanIdx,:));

    % Time domain data plot
    if(ishandle(ui.tsPanel(chanIdx).hplot))
        set(ui.tsPanel(chanIdx).hplot, {'Ydata'}, num2cell([currChDataI(:) currChDataQ(:)].', 2));
    end
end

%  -----------------------------------------------------------------------
%  Description:    Plot range FFT results using given one chirp of range profile data.
%
%  Input:      chanIdx - RX channel index
%              chirpIdx - one chirp of range profile data for all
%                         RX channels.
%              linearMode - Two options of the plot: linear and log
%  -----------------------------------------------------------------------
function [] = ui_plot1DRangeProfile(chanIdx, chirpIdx, linearMode)
    global ui
    global dataSet

    rangeProfileData = dataSet.radarCubeData(chirpIdx, chanIdx , :);

    if (linearMode == 1)
        channelData = abs(rangeProfileData(:));
    else
        chFFT = rangeProfileData(:);
        channelData = 20 * log10 (abs(chFFT));
    end

    % Plot channel range plot data
    if(ishandle(ui.rpPanel(chanIdx).hplot))
        set(ui.rpPanel(chanIdx).hplot, {'Ydata'}, {channelData});
    end
end

%  -----------------------------------------------------------------------
%  Description: Display configuration parameters on display panel
%  -----------------------------------------------------------------------
function ui_displayConfParams()
global Params
global ui

    radarCubeParams = Params.radarCubeParams;
    adcDataParams =Params.adcDataParams;
    RFParams = Params.RFParams;

    dat =  {'Start Frequency (Ghz)', RFParams.startFreq;...
            'Slope (MHz/us)', RFParams.freqSlope;...
            'Samples per chirp', adcDataParams.numAdcSamples;...
            'Chirps per frame',  adcDataParams.numChirpsPerFrame;...
            'Sampling rate (Msps)', RFParams.sampleRate;...
            'Bandwidth (GHz)', RFParams.bandwidth;...
            'Range resolution (m)', RFParams.rangeResolutionsInMeters;...
            'Velocity resolution (m/s)', RFParams.dopplerResolutionMps;...
            'Number of Tx (MIMO)', radarCubeParams.numTxChan;...
            'Frame periodicity (msec)', RFParams.framePeriodicity;...
            };

    columnname =   {'             Parameter (Units)             ', '   Value   '};
    columnformat = {'char', 'numeric'};

    uitable(ui.hChirpParamPanel,'Units','normalized','Position',...
                [0.1 0.05 0.8 0.9], 'Data', dat,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnWidth', 'auto',...
                'RowName',[]);

     % Convert dataFmt to string
     if((adcDataParams.dataFmt == 0) || ...
        (adcDataParams.dataFmt == 3))
         dataFmtString = 'Real';
     else
          dataFmtString = 'Complex';
     end

      % Convert iqswap mode to string
     if(adcDataParams.iqSwap == 0)
         iqSwapString = 'IQ';
     else
         iqSwapString = 'QI';
     end

     % Convert interleave mode to string
     if(adcDataParams.chanInterleave == 0)
         chanInterleaveString = 'Interleave';
     else
         chanInterleaveString = 'non-Interleave';
     end

      % Convert adcBits to string
     if(adcDataParams.adcBits == 0)
         adcBitsString = '12 bits';
     elseif(adcDataParams.adcBits == 1)
         adcBitsString = '14 bits';
     elseif(adcDataParams.adcBits == 2)
         adcBitsString = '16 bits';
     else
         error 'adcBits (%d) is not supported';
     end

     adcDataTbl =  {'dataFmt', dataFmtString;...
        'iqSwap', iqSwapString;...
        'chanInterleave', chanInterleaveString;...
        'numChirpsPerFrame',  num2str(adcDataParams.numChirpsPerFrame);...
        'adcBits', adcBitsString;...
        'numRxChan', num2str(adcDataParams.numRxChan);...
        'numAdcSamples', num2str(adcDataParams.numAdcSamples);
        };

    columnname =   {'             Parameter (Units)             ', '    Value    '};
    columnformat = {'char', 'numeric'};

    uitable(ui.hAdcParamPanel,'Units','normalized','Position',...
                [0.1 0.05 0.8 0.9], 'Data', adcDataTbl,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnWidth', 'auto',...
                'RowName',[]);

    if(radarCubeParams.iqSwap == 0)
        iqSwapString = 'IQ';
    else
        iqSwapString = 'QI';
    end

    if(radarCubeParams.radarCubeFmt == 1)
        radarFmtString = 'FMT1';
    else
        radarFmtString = 'UnKnown';
    end
    radarCubeTbl =  {'iqSwap', iqSwapString;...
                    'numRxChan', num2str(radarCubeParams.numRxChan);...
                    'numTxChan', num2str(radarCubeParams.numTxChan);...
                    'numRangeBins', num2str(radarCubeParams.numRangeBins);...
                    'numDopplerChirps', num2str(radarCubeParams.numDopplerChirps);...
                    'radarCubeFmt', radarFmtString;...
                    };

    columnname =   {'             Parameter (Units)             ', '   Value   '};
    columnformat = {'char', 'numeric'};

    uitable(ui.hRadarCubeParamPanel,'Units','normalized','Position',...
                [0.1 0.05 0.8 0.9], 'Data', radarCubeTbl,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnWidth', 'auto',...
                'RowName',[]);
end
