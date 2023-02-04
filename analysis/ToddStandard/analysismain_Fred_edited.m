
%%
%*************************************************************************
% Initializations
%*************************************************************************


% define plot color sequence, axis fonts
PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 16)

colormap([0 0 0])
scrsz = get(0,'ScreenSize');

% jauimodel stuff
loader = edu.washington.rieke.Analysis.getEntityLoader(); 
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();


%Data and export folder paths
dataFolder = '/Users/toddappleby/Documents/Data/2019_0620';
exportFolder = '/Users/toddappleby/Documents/Exports/2019_0620/Bc1/';

import auimodel.*
import vuidocument.*

%%
params.FrequencyCutoff = 500;
params.Amp = 'Amp1';
params.Verbose = 1;
params.DecimatePts = 200;
params.SamplingInterval = 0.0001; 
params.FrameRate = 60.1830;
params.OvationFlag = 1;
params.SpatialFlag = 0;
params.SaveToIgor = 0;
params.SaveGraphs = 0;
params.rootDir = '/Users/toddappleby/Documents/Data/Figures/';

params.initialSTD = 5; %initialize STD for spike detection at some value
%%
%*************
%noise test loader
%*************
list = loader.loadEpochList([exportFolder 'MotionNoise.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(backgroundClass)'});

gui = epochTreeGUI(tree);
%%
%*************
%Single Spot loader
%*************
list = loader.loadEpochList([exportFolder 'SingleSpot.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulus:Amp1:offset)','protocolSettings(backgroundIntensity)','protocolSettings(spotDiameter)'});

gui = epochTreeGUI(tree);
%% for single spot analysis
params.CellAttached = 1;
params.plotOffset = 200;    

node = gui.getSelectedEpochTreeNodes;

result = SingleSpotAnalysis(node, params);
%%
%*******
% Contrast Response
%*******
list = loader.loadEpochList([exportFolder 'ContrastResp.mat'], dataFolder);
dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(contrast)'});

gui = epochTreeGUI(tree);
%%
% Contrast Response Analysis
node = gui.getSelectedEpochTreeNodes;
results = ContrastResp(node, params);

%%
%*********
%RigB LED Loader
%*********
list = loader.loadEpochList([exportFolder 'LEDPulse.mat'], dataFolder);
dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(background:Amp1:value)','protocolSettings(tailTime)'});

gui = epochTreeGUI(tree);

%%
%*****
%LED Analysis
%*****
params.CellAttached = 1;

node = gui.getSelectedEpochTreeNodes;

results = SimpleLED(node, params);

%%
%***********
%MTF loader
%***********
%                                                                                                                 
list = loader.loadEpochList([exportFolder 'sMTF.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulusClass)','protocolSettings(temporalFrequency)','protocolSettings(radius)'});

gui = epochTreeGUI(tree);
%% for MTF analysis
params.CellAttached = 1; 
params.plotOffset = 200;
params.color = 'b';


node = gui.getSelectedEpochTreeNodes;

[result,diag] = MTFAnalysis(node, params);
%%
%*************************************************************************
% chirp
%*************************************************************************
list = loader.loadEpochList([exportFolder 'Chirp.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'protocolSettings(source:parent:parent:description)','cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulus:Amp1:offset)','protocolSettings(background:Microdisplay_Stage@localhost:microdisplayBrightness)',keywordSplitter_java});

gui = epochTreeGUI(tree);

%%

params.SpatialFlag = 1;
params.CellAttached = 1;
params.fileName = 'chirp';

node = gui.getSelectedEpochTreeNodes;


results = AnalyzeChirpStimulus(node{1}, params);
%%
%*************************************************************************
% anything oriented 
%*************************************************************************

list = loader.loadEpochList([exportFolder 'GratingDSOS.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java =  riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulus:Amp1:offset)','protocolSettings(intensity)','protocolSettings(barWidth)','protocolSettings(apertureRadius)','protocolSettings(temporalFrequency)','protocolSettings(orientation)'});

gui = epochTreeGUI(tree);

%%

params.CellAttached = 1;  
params.plotOffset = 200;-
params.fileName = 'DS-WC';


clear Resp orientation AveResponse;
node = gui.getSelectedEpochTreeNodes;

results = AnalyzeOrientedStimuli(node, params);

%%
%*************************************************************************
% Doves
%*************************************************************************

list = loader.loadEpochList([exportFolder 'Doves.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulus:Amp1:offset)','protocolSettings(stimulusIndex)'});

gui = epochTreeGUI(tree);

%%
params.CellAttached = 1;
params.plotOffset = 300;
params.fileName = 'Doves';
params.plotColors = 'kkkkkkkkkkkkkkkkk'

node = gui.getSelectedEpochTreeNodes;

MeanSelectedNodes(node, params);

%%
%*************************************************************************
% OMS
%*************************************************************************

list = loader.loadEpochList([exportFolder 'OMSGrating.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulus:Amp1:offset)','protocolSettings(stimulusClass)'});

gui = epochTreeGUI(tree);

%%
%TEXTURE FIRST IF GRAPHING
params.CellAttached = 1;
params.plotOffset = 120;
params.fileName = 'OMS';
params.plotColors = 'bgkr';

node = gui.getSelectedEpochTreeNodes;

results = MeanSelectedNodes(node, params);

%%
%*************************************************************************
% Expanding spots
%*************************************************************************
list = loader.loadEpochList([exportFolder 'ExpandingSpots.mat'], dataFolder);
   

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulus:Amp1:offset)','protocolSettings(spotIntensity)','protocolSettings(currentSpotSize)'});

gui = epochTreeGUI(tree);


%%

params.CellAttached = 1;
params.fileName = 'ExpandingSpots';
params.plotOffset = 250;
params.plotColors = 'kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk';

node = gui.getSelectedEpochTreeNodes;

results = MeanSelectedNodes(node, params);  





