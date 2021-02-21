setup_multiregion_allstdp;  


RecordingSettings.saveDir = '~/Documents/MATLAB/Vertex_Results/VERTEX_results_multiregion/mr4reg_diamond_nodelay_stdp_nostim_05startweight_3kconnects';

RecordingSettings.LFP = true;  
RecordingSettings.weights_preN_IDs = [1:966:3864,3865:171:4546,4547:50:5000]; 
[meaX, meaY, meaZ] = meshgrid(0:100:2000, 200, 195:-95:5); 
RecordingSettings.meaXpositions = meaX; 
RecordingSettings.meaYpositions = meaY;
RecordingSettings.meaZpositions = meaZ;
RecordingSettings.minDistToElectrodeTip = 20; 
RecordingSettings.v_m = 100:200:1200; 
RecordingSettings.maxRecTime = 500; 
RecordingSettings.sampleRate = 1000; 
RecordingSettings.stdpvars = [1:2:500];
 
SimulationSettings.simulationTime = 500;
SimulationSettings.timeStep = 0.03125; 
SimulationSettings.parallelSim = false; 


RecordingSettings.weights_arr = [1 (SimulationSettings.simulationTime/SimulationSettings.timeStep)-1]; 



NeuronParams(1).Input(2).inputType = 'i_step';
NeuronParams(1).Input(2).timeOn = 200;
NeuronParams(1).Input(2).timeOff = 250;
NeuronParams(1).Input(2).amplitude = 1000;



%% initialise the network

[params, connections, electrodes] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

 clear TissueParams;

          
%% reset the second region          
setup_multiregion_allstdp;
[params2, connections2, electrodes2] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);
       

 %% reset the third region
% setup_multiregion_gexp;
% setup_multiregion_withinboundconnection;
setup_multiregion_allstdp;
[TissueParams.StimulationField, TissueParams.model] = invitroSliceStim('farapartlectrodesbig.stl',10000);
TissueParams.StimulationOn = [100:20:400];
TissueParams.StimulationOff = [105:20:405];
[params3, connections3, electrodes3] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);
clear Tissueparams;
 %% reset the forth region
% setup_multiregion_gexp;
% setup_multiregion_withinboundconnection;
% setup_multiregion_allstdp;

[params4, connections4, electrodes4] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);          

%% reset the fifth region
% setup_multiregion_gexp;
% setup_multiregion_withinboundconnection;
% setup_multiregion_allstdp;

[params5, connections5, electrodes5] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

 %% Set the connectivity between regions

 % defining the between region connectivity here:
 regionConnect.map = [0,1,0,0,0;
                      0,0,1,0,0;
                      0,0,0,1,0;
                      0,0,0,0,1;
                      0,0,0,0,0];

 regionConnect.exportingNeuronPops{1} = 1; % 1 represents neuron population 1, the excitatory Pyramidal cells.
 regionConnect.exportingNeuronPops{2} = 1;
 regionConnect.exportingNeuronPops{3} = 1;
 regionConnect.exportingNeuronPops{4} = 1;
 regionConnect.exportingNeuronPops{5} = 1;
 
 % identify which neuron pops are designated as dummy neurons to just
 % recieve external signals. Use [] if no dummy neurons are present.
 regionConnect.dummyNeuronPops{1} = 3;
 regionConnect.dummyNeuronPops{2} = 3;
 regionConnect.dummyNeuronPops{3} = 3;
 regionConnect.dummyNeuronPops{4} = 3;
 regionConnect.dummyNeuronPops{5} = 3;
 
 
 %% Run the simulation
 
 %stack the parameters for params, connections and electrodes into cell
 % arrays which can then be fed into runSimulationMultiregional
 paramStacked = {params, params2, params3, params4, params5};
 connectionStacked = {connections,connections2, connections3, connections4, connections5};
 electrodeStacked = {electrodes,electrodes2, electrodes3, electrodes4, electrodes5};
 
tic
runSimulationMultiregional(paramStacked,connectionStacked,electrodeStacked,regionConnect);
toc


%% Plotting!

% need to use a Multiregions variant of loadResults to load the results for
% every region in one structure: 
oneStimMultiRegResults = loadResultsMultiregions(RecordingSettings.saveDir);

Region1A=mean(oneStimMultiRegResults(1).LFP);
Region2A=mean(oneStimMultiRegResults(2).LFP);
Region3A=mean(oneStimMultiRegResults(3).LFP);
Region4A=mean(oneStimMultiRegResults(4).LFP);
Region5A=mean(oneStimMultiRegResults(5).LFP);

% plot the spike raster for each region
 plotSpikeRaster(oneStimMultiRegResults(1))
 title('Region 1 Spike Raster One Stimulation')
 xlabel('time (ms)');
 ylabel('Neuron ID');
 
 plotSpikeRaster(oneStimMultiRegResults(2))
 title('Region 2 Spike Raster One Stimulation')
 xlabel('time (ms)');
 ylabel('Neuron ID');
 
 plotSpikeRaster(oneStimMultiRegResults(3))
 title('Region 3 Spike Raster One Stimulation')
 xlabel('time (ms)');
 ylabel('Neuron ID');
 
 plotSpikeRaster(oneStimMultiRegResults(4))
 title('Region 4 Spike Raster One Stimulation')
 xlabel('time (ms)');
 ylabel('Neuron ID');
 
 plotSpikeRaster(oneStimMultiRegResults(5))
 title('Region 5 Spike Raster One Stimulation')
 xlabel('time (ms)');
 ylabel('Neuron ID');
 
figure
plot(mean(oneStimMultiRegResults(1).LFP))
title('LFP of stimulated artecraft ')
set(gca,'XLim',[20 100]);
set(gca,'XTick',[20 30 40 45 50 60 70 80 100]);
% set(gca, 'YLim', [-0.1 0.1]);
% set(gca,'YTick',[-0.1 0 0.1]);
 
 
% % plot the mean LFP for both regions, and the difference between them
 figure
 subplot(311)
 plot(mean(oneStimMultiRegResults(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(312)
 plot(mean(oneStimMultiRegResults(2).LFP))
 title('Region 2 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(313)
 plot(mean(oneStimMultiRegResults(1).LFP) - mean(oneStimMultiRegResults(2).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');

 figure
 subplot(411)
 plot(mean(oneStimMultiRegResults(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(412)
 plot(mean(oneStimMultiRegResults(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(413)
 plot(mean(oneStimMultiRegResults(1).LFP) - mean(oneStimMultiRegResults(3).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');


 figure
 subplot(511)
 plot(mean(oneStimMultiRegResults(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(512)
 plot(mean(oneStimMultiRegResults(4).LFP))
 title('Region 4 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(513)
 plot(mean(oneStimMultiRegResults(1).LFP) - mean(oneStimMultiRegResults(4).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');

 figure
 subplot(611)
 plot(mean(oneStimMultiRegResults(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(612)
 plot(mean(oneStimMultiRegResults(5).LFP))
 title('Region 5 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(613)
 plot(mean(oneStimMultiRegResults(1).LFP) - mean(oneStimMultiRegResults(5).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');

 
 figure
 plot(mean(oneStimMultiRegResults(1).LFP),'LineWidth',2)
 title ('Five Region average LFP Comparison')
  hold on
 plot(mean(oneStimMultiRegResults(2).LFP),'LineWidth',2)
 plot(mean(oneStimMultiRegResults(3).LFP),'LineWidth',2)
 plot(mean(oneStimMultiRegResults(4).LFP),'LineWidth',2)
 plot(mean(oneStimMultiRegResults(5).LFP),'LineWidth',2)
 legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 hold off
  

% get the weights for the whole network at the first and last time
% snapshots in a plottable form:
 Ar1time1weights=getSparseConnectivityWeights(oneStimMultiRegResults(1).weights_arr{1},oneStimMultiRegResults(1).syn_arr,oneStimMultiRegResults(1).params.TissueParams.N);
 Ar1time2weights=getSparseConnectivityWeights(oneStimMultiRegResults(1).weights_arr{2},oneStimMultiRegResults(1).syn_arr,oneStimMultiRegResults(1).params.TissueParams.N);
 Ar2time1weights=getSparseConnectivityWeights(oneStimMultiRegResults(2).weights_arr{1},oneStimMultiRegResults(2).syn_arr,oneStimMultiRegResults(2).params.TissueParams.N);
 Ar2time2weights=getSparseConnectivityWeights(oneStimMultiRegResults(2).weights_arr{2},oneStimMultiRegResults(2).syn_arr,oneStimMultiRegResults(2).params.TissueParams.N);
 Ar3time1weights=getSparseConnectivityWeights(oneStimMultiRegResults(3).weights_arr{1},oneStimMultiRegResults(3).syn_arr,oneStimMultiRegResults(3).params.TissueParams.N);
 Ar3time2weights=getSparseConnectivityWeights(oneStimMultiRegResults(3).weights_arr{2},oneStimMultiRegResults(3).syn_arr,oneStimMultiRegResults(3).params.TissueParams.N);
 Ar4time1weights=getSparseConnectivityWeights(oneStimMultiRegResults(4).weights_arr{1},oneStimMultiRegResults(4).syn_arr,oneStimMultiRegResults(4).params.TissueParams.N);
 Ar4time2weights=getSparseConnectivityWeights(oneStimMultiRegResults(4).weights_arr{2},oneStimMultiRegResults(4).syn_arr,oneStimMultiRegResults(4).params.TissueParams.N);
 Ar5time1weights=getSparseConnectivityWeights(oneStimMultiRegResults(5).weights_arr{1},oneStimMultiRegResults(5).syn_arr,oneStimMultiRegResults(5).params.TissueParams.N);
 Ar5time2weights=getSparseConnectivityWeights(oneStimMultiRegResults(5).weights_arr{2},oneStimMultiRegResults(5).syn_arr,oneStimMultiRegResults(5).params.TissueParams.N);
% plot the weight differences between the start and end of the simulation
% for each network. Colours represent the weight changes.

figure
imagesc(log(abs(Ar1time2weights-Ar1time1weights)))
title('Weight changes between different neurons in Region1');
colorbar;
caxis([-2 2.5])


figure
imagesc(log(abs(Ar2time2weights-Ar2time1weights)))
title('Weight changes between different neurons in Region2');
colorbar;
caxis([-2 2.5])

figure
imagesc(log(abs(Ar3time2weights-Ar3time1weights)))
title('Weight changes between different neurons in Region3');
colorbar;
caxis([-2 2.5])


figure
imagesc(log(abs(Ar4time2weights-Ar4time1weights)))
title('Weight changes between different neurons in Region4');
colorbar;
caxis([-2 2.5])


figure
imagesc(log(abs(Ar5time2weights-Ar5time1weights)))
title('Weight changes between different neurons in Region5');
colorbar;
caxis([-2 2.5])


figure
plot(mean(oneStimMultiRegResults(1).weights{1,:})/SimulationSettings.simulationTime);xlim([0 SimulationSettings.simulationTime])
hold on
plot(mean(oneStimMultiRegResults(2).weights{1,:})/SimulationSettings.simulationTime);xlim([0 SimulationSettings.simulationTime])
plot(mean(oneStimMultiRegResults(3).weights{1,:})/SimulationSettings.simulationTime);xlim([0 SimulationSettings.simulationTime])
plot(mean(oneStimMultiRegResults(4).weights{1,:})/SimulationSettings.simulationTime);xlim([0 SimulationSettings.simulationTime])
plot(mean(oneStimMultiRegResults(5).weights{1,:})/SimulationSettings.simulationTime);xlim([0 SimulationSettings.simulationTime])
title ('One stimulation average synaptic weight');
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast');
xlabel('time (ms)');
ylabel('Mean synaptic weight (pA)');
set(gca,'XLim',[50 200]);
set(gca,'XTick',[0 50 100 150 200]);
hold off

figure
plot(groupRates(oneStimMultiRegResults(1), 30, 100));
hold on
plot(groupRates(oneStimMultiRegResults(2), 30, 100));
plot(groupRates(oneStimMultiRegResults(3), 30, 100));
plot(groupRates(oneStimMultiRegResults(4), 30, 100));
plot(groupRates(oneStimMultiRegResults(5), 30, 100));
title ('Group average firing rate in each region');
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast');
set(gca,'XLim',[1 3]);
set(gca,'XTick',[1 2 3]);
xlabel('Neuron Group Number');
ylabel('Mean firing rate');
hold off










