%randomized rng
rng('shuffle');
%clear previous figures
close all;
clear all;

%% List of all device {deviceName variantName}
deviceNames{11} = {'ptm130' 'W20uL1u_60t'};
deviceNames{12} = {'ptm130' 'W20uL1u_50t'};

deviceNames{21} = {'d1n4148' '5'};
deviceNames{22} = {'d1n4148' '8'};
deviceNames{23} = {'d1n4148' '10'};

deviceNames{31} = {'d1n3858' '4'};
deviceNames{32} = {'d1n3858' '7'};
deviceNames{33} = {'d1n3858' '10'};

deviceNames{41} = {'ptm180' 'WL'};

%% List of all reference Model
referenceNames{11} = 'ptm130_W20uL1u_5000t';
referenceNames{12} = 'irff130_ref_57800_t';
referenceNames{21} = 'd1n4148_ref';
referenceNames{31} = 'd1n3858_ref';

referenceNames{40} = 'ptm180_NMOS_W18,5L4,5';
referenceNames{41} = 'ptm180_NMOS_W93,5L4,5';
referenceNames{42} = 'ptm180_NMOS_W130L0,24';
referenceNames{43} = 'ptm180_NMOS_W193,5L4,5';
referenceNames{44} = 'ptm180_NMOS_W100000L0,18';
referenceNames{45} = 'ptm180_PMOS_W93,5L4,5';
referenceNames{46} = 'ptm180_PMOS_W360L0,24';
referenceNames{47} = 'ptm180_PMOS_W393,5L4,5';
referenceNames{48} = 'ptm180_PMOS_W400000L0,18';


%% Settings %%
% Select device and its reference model
deviceName = deviceNames{41}{1};
variantName = deviceNames{41}{2};
ModelName = [deviceName '_' variantName];

referenceName = referenceNames{48};

%Switch for Simulated Annealing
StartSA = 0;

%% Initialization

settingsFilename = 'settings.m';

paths.mainPath = pwd;  
paths.dataPath = fullfile(paths.mainPath,'data');
paths.examplesPath = fullfile(paths.mainPath,'examples');
paths.devicePath = fullfile(paths.examplesPath,deviceName);
paths.IsDiode = fullfile(paths.devicePath,'IsDiode');
paths.IsTransistor = fullfile(paths.devicePath,'IsTransistor');

if exist(paths.IsDiode) == 2
    
    paths.referenceSimDataPath = fullfile(paths.devicePath,[referenceName '.txt']);
    paths.resultPath = fullfile(paths.devicePath,[ModelName '.xml']);
    
    disp(['INFO: Reference model: ' referenceName]);
    disp(['INFO: Initial model: ' ModelName]);
    
    generate_diode_xml(deviceName, variantName, paths);
    
    disp(['DONE: ' '"' deviceName '_' variantName '" generated.']);

elseif exist(paths.IsTransistor) == 2
    
    paths.variantPath = fullfile(paths.devicePath,ModelName);
    
    paths.referenceModelPath = fullfile(paths.devicePath,referenceName);
    
    paths.initialModelPath = fullfile(paths.variantPath,[ModelName '.xml']);
    paths.settingsPath = fullfile(paths.variantPath,settingsFilename);
    paths.initialSimDataPath = fullfile(paths.variantPath,[ModelName '.txt']);

    paths.resultPath = fullfile(paths.dataPath,ModelName);

    if exist(paths.resultPath) == 0
        mkdir(paths.resultPath);
    end

    %Check if PWLModel exist, if not try to create it from the simulation data
    if exist(paths.initialModelPath) ~= 2
        if (exist(paths.initialSimDataPath, 'file') ~= 2)
            error(['No simulation data or xml available for the device "' ModelName '"!']);
        else
            % Parse simulation data and create .xml
            simData2XML(paths.initialSimDataPath, paths.initialModelPath, deviceName, variantName)
        end
    end
    
     if exist([paths.referenceModelPath,'.xml']) ~= 2
        if (exist(paths.initialSimDataPath, 'file') ~= 2)
            error(['No reference data or xml available for the device "' ModelName '"!']);
        else
            % Parse reference data and create .xml
            simData2XML([paths.referenceModelPath,'.txt'], [paths.referenceModelPath,'.xml'], deviceName, variantName)
        end
    end

    %% start SA run for optimize PWL model
    if StartSA
        disp(['INFO: Reference model: ' referenceName]);
        disp(['INFO: Initial model: ' ModelName]);
        [ PWLModel, PWLCharacteristics, variantName ] = simulatedannealing(paths);

        %write XML file
        variantName=[char(variantName), '_optimized'];
        paths.outputFilePath= fullfile(paths.resultPath,[ModelName '_optimized.xml']);
        write_xml(PWLModel.Points, PWLModel.Triangles, paths.outputFilePath, deviceName, variantName);

        disp(['DONE: ' '"' deviceName '_' variantName '" generated.']);
    end

end

function simData2XML(simDataPath, outputPath, deviceName, variantName)
    P = dlmread(simDataPath,'	',1,1);
    P((P(:,3)==0),:)=[];
    [~,I]=sort(P(:,1));
    P=P(I,:);
    F = delaunay(P(:,1),P(:,2));

    write_xml(P, F, outputPath, deviceName, variantName);
end