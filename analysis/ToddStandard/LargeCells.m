%%
%*************************************************************************
% Initializations
%*************************************************************************


% define plot color sequence, axis fonts
PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 16)

colormap([0 0 0])
scrsz = get(0, 'ScreenSize');

% jauimodel stuff
loader = edu.washington.rieke.Analysis.getEntityLoader(); 
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();

%Data and export folder paths
dataFolder = '/Users/fred/Dropbox/Data/';
exportFolder = '/Users/fred/Dropbox/Data/';

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
params.rootDir = '~/Dropbox/LargeCells/';

%%
%*************************************************************************
% chirp
%*************************************************************************
list = loader.loadEpochList([exportFolder 'ribeye.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'protocolSettings(source:parent:parent:description)','cell.label','protocolSettings(background:Microdisplay_Stage@localhost:microdisplayBrightness)',keywordSplitter_java});

gui = epochTreeGUI(tree);

%%

params.SpatialFlag = 0;
params.CellAttached = 0;
params.fileName = 'chirp';

node = gui.getSelectedEpochTreeNodes;

results = AnalyzeChirpStimulus(node{1}, params);

%%
%*************************************************************************
% moving bars
%*************************************************************************

list = loader.loadEpochList([exportFolder 'MovingBars.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(orientation)'});

gui = epochTreeGUI(tree);

%%

params.CellAttached = 1;
params.plotOffset = 80;
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

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:recordingTechnique)','protocolSettings(stimulusIndex)'});

gui = epochTreeGUI(tree);

%%
params.CellAttached = 1;
params.plotOffset = 200;
params.fileName = 'Doves';
params.plotColors = 'kkkkkkkkkkkkkkkkk'

node = gui.getSelectedEpochTreeNodes;

MeanSelectedNodes(node, params);

%%
%*************************************************************************
% OMS
%*************************************************************************

list = loader.loadEpochList([exportFolder 'OMS.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label',keywordSplitter_java,'protocolSettings(stimulusClass)'});

gui = epochTreeGUI(tree);

%%

params.CellAttached = 1;
params.plotOffset = 120;
params.fileName = 'OMS';
params.plotColors = 'bgkr'

node = gui.getSelectedEpochTreeNodes;

MeanSelectedNodes(node, params);

%%
%*************************************************************************
% Expanding spots
%*************************************************************************
list = loader.loadEpochList([exportFolder 'ExpandingSpots.mat'], dataFolder);

dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);

tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label',keywordSplitter_java,'protocolSettings(currentSpotSize)'});

gui = epochTreeGUI(tree);


%%

params.CellAttached = 1;
params.fileName = 'ExpandingSpots';
params.plotOffset = 200;
params.plotColors = 'kkkkkkkkkkkk'

node = gui.getSelectedEpochTreeNodes;

results = MeanSelectedNodes(node, params);





