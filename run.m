% load run_config.json
runConfig = jsondecode(fileread("run_config.json"));

% open setupJson
setupJson = jsondecode(fileread(runConfig.setupJsonFileName));

% change $(setupJsonFileName).configUsed to $signalConfigJsonFileName
setupJson.configUsed = runConfig.signalConfigJsonFileName;

% change $(setupJsonFileName).processedFileName to $binFilename and remove other fields
numBinFiles = length(runConfig.binFilenames);
setupJson.capturedFiles.files = [];
try
    setupJson.capturedFiles = rmfield(setupJson.capturedFiles, "fileBasePath");
catch
end
for file_i = 1:numBinFiles
    setupJson.capturedFiles.files{file_i} = struct("processedFileName",runConfig.binFilenames{file_i});
end

% save setupJson
fid = fopen(runConfig.setupJsonFileName, 'w');
fprintf(fid, '%s', jsonencode(setupJson, "PrettyPrint", true));

% Reset
clear;



% load run_config.json
runConfig = jsondecode(fileread("run_config.json"));

rawDataReader(runConfig.setupJsonFileName, runConfig.rawDataFileName, runConfig.radarCubeDataFileName, runConfig.createDebugPlot);

