multiRegion_multilayer_stdp;  

RecordingSettings.saveDir = '~/Documents/MATLAB/Vertex_Results/VERTEX_results_multiregion/mr4reg_diamond_nodelay_stdp_nostim_05startweight_3kconnects';

RecordingSettings.LFP = true; 
RecordingSettings.weights_preN_IDs = [1:966:3864,3865:171:4546,4547:50:5000]; 
[meaX, meaY, meaZ] = meshgrid(0:500:2000, 200, 700:-100:0); 
RecordingSettings.meaXpositions = meaX; 
RecordingSettings.meaYpositions = meaY;
RecordingSettings.meaZpositions = meaZ;
RecordingSettings.minDistToElectrodeTip = 20; 
RecordingSettings.v_m = 250:250:4750; 
RecordingSettings.maxRecTime = 500; 
RecordingSettings.sampleRate = 1000; 
SimulationSettings.simulationTime = 500; 
SimulationSettings.timeStep = 0.03125; 
SimulationSettings.parallelSim = false; 
RecordingSettings.weights_arr = [1 (SimulationSettings.simulationTime/SimulationSettings.timeStep)-1]; % simulation steps




% NeuronParams(1).Input(2).inputType = 'i_step';
% NeuronParams(1).Input(2).timeOn = 200;
% NeuronParams(1).Input(2).timeOff = 250;
% NeuronParams(1).Input(2).amplitude = 1000; 


%% initialise the network

[params, connections, electrodes] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

%% reset the second region          
multiRegion_multilayer_stdp;


[params2, connections2, electrodes2] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);



 %% reset the third region
multiRegion_multilayer_stdp;

[params3, connections3, electrodes3] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

%% reset the forth region
multiRegion_multilayer_stdp;

[params4, connections4, electrodes4] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);
          

%% reset the fifth region
multiRegion_multilayer_stdp;

[params5, connections5, electrodes5] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

 %% Set the connectivity between regions

 % defining the between region connectivity here:
 regionConnect.map = [0,1,0,0,0;
                      0,0,1,0,0;
                      0,0,0,1,0;
                      0,0,0,0,1;

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
 connectionStacked = {connections, connections2, connections3, connections4, connections5};
 electrodeStacked = {electrodes,electrodes2,electrodes3,electrodes4,electrodes5};
 
tic
runSimulationMultiregional(paramStacked,connectionStacked,electrodeStacked,regionConnect);
toc

%% Plotting!

% need to use a Multiregions variant of loadResults to load the results for
% every region in one structure: 
noStimRegResults = loadResultsMultiregions(RecordingSettings.saveDir);

Region1B=mean(noStimRegResults(1).LFP);
Region2B=mean(noStimRegResults(2).LFP);
Region3B=mean(noStimRegResults(3).LFP);
Region4B=mean(noStimRegResults(4).LFP);
Region5B=mean(noStimRegResults(5).LFP);


 

% % get the weights for the whole network at the first and last time
% % snapshots in a plottable form:
 Br1time1weights=getSparseConnectivityWeights(noStimRegResults(1).weights_arr{1},noStimRegResults(1).syn_arr,noStimRegResults(1).params.TissueParams.N);
 Br1time2weights=getSparseConnectivityWeights(noStimRegResults(1).weights_arr{2},noStimRegResults(1).syn_arr,noStimRegResults(1).params.TissueParams.N);
 Br2time1weights=getSparseConnectivityWeights(noStimRegResults(2).weights_arr{1},noStimRegResults(2).syn_arr,noStimRegResults(2).params.TissueParams.N);
 Br2time2weights=getSparseConnectivityWeights(noStimRegResults(2).weights_arr{2},noStimRegResults(2).syn_arr,noStimRegResults(2).params.TissueParams.N);
 Br3time1weights=getSparseConnectivityWeights(noStimRegResults(3).weights_arr{1},noStimRegResults(3).syn_arr,noStimRegResults(3).params.TissueParams.N);
 Br3time2weights=getSparseConnectivityWeights(noStimRegResults(3).weights_arr{2},noStimRegResults(3).syn_arr,noStimRegResults(3).params.TissueParams.N);
 Br4time1weights=getSparseConnectivityWeights(noStimRegResults(4).weights_arr{1},noStimRegResults(4).syn_arr,noStimRegResults(4).params.TissueParams.N);
 Br4time2weights=getSparseConnectivityWeights(noStimRegResults(4).weights_arr{2},noStimRegResults(4).syn_arr,noStimRegResults(4).params.TissueParams.N);
 Br5time1weights=getSparseConnectivityWeights(noStimRegResults(5).weights_arr{1},noStimRegResults(5).syn_arr,noStimRegResults(5).params.TissueParams.N);
 Br5time2weights=getSparseConnectivityWeights(noStimRegResults(5).weights_arr{2},noStimRegResults(5).syn_arr,noStimRegResults(5).params.TissueParams.N);
%  
% % 
% % % plot the weight differences between the start and end of the simulation
% % % for each network. Colours represent the weight changes.
figure
imagesc(log(abs(Br1time2weights-Br1time1weights)))
title('Weight changes between different neurons in Region1');
colorbar;



figure
imagesc(log(abs(Br2time2weights-Br2time1weights)))
title('Weight changes between different neurons in Region2');
colorbar;


figure
imagesc(log(abs(Br3time2weights-Br3time1weights)))
title('Weight changes between different neurons in Region3');
colorbar;



figure
imagesc(log(abs(Br4time2weights-Br4time1weights)))
title('Weight changes between different neurons in Region4');
colorbar;



figure
imagesc(log(abs(Br5time2weights-Br5time1weights)))
title('Weight changes between different neurons in Region5');
colorbar;




rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 1 No Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 511;
plotSpikeRaster(noStimRegResults(1), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 2 No Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 512;
plotSpikeRaster(noStimRegResults(2), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 3 No Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 513;
plotSpikeRaster(noStimRegResults(3), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 4 No stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 514;
plotSpikeRaster(noStimRegResults(4), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster 5';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 515;
plotSpikeRaster(noStimRegResults(5), rasterParams);


figure
plot(mean(noStimRegResults(1).LFP(:,1:490)),'LineWidth',2)
hold on
plot(mean(noStimRegResults(2).LFP(:,1:490)),'LineWidth',2)
plot(mean(noStimRegResults(3).LFP(:,1:490)),'LineWidth',2)
plot(mean(noStimRegResults(4).LFP(:,1:490)),'LineWidth',2)
plot(mean(noStimRegResults(5).LFP(:,1:490)),'LineWidth',2)
hold off
title ('Five Region average LFP Comparison')
set(gcf,'color','w');
set(gca,'FontSize',16)
xlabel('Time (ms)', 'FontSize', 16)
ylabel('LFP (mV)', 'FontSize', 16)
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')






 
% lmax=mean(noStimRegResults(1).LFP(1,:));
% imax=1;
% for i=1:40
%     if mean(noStimRegResults(1).LFP(i,:)) > lmax
%         imax=i;
%         lmax=mean(noStimRegResults(1).LFP(i,:));
%     end
% end
%    
%         
% lmin=mean(noStimRegResults(1).LFP(1,:));
% imin=1;
% for i=1:40
%     if mean(noStimRegResults(1).LFP(i,:)) < lmin
%         imin=i;
%         lmin=mean(noStimRegResults(1).LFP(i,:));
%     end
% end


groupRates(noStimRegResults(1), 0, 500)
groupRates(noStimRegResults(2), 0, 500)
groupRates(noStimRegResults(3), 0, 500)
groupRates(noStimRegResults(4), 0, 500)
groupRates(noStimRegResults(5), 0, 500)



% ElecArray=[1,16,27,40];
% for j=1:ElecArray
%      for i=1:5
%      Electrode{j}=[noStimRegResults(i).LFP(1,150:450)];
%      var(Electrode{j})
%      end
% end


for t=1:500
   spikes = noStimRegResults(2).spikes(noStimRegResults(2).spikes(:,2)>=0 & ...
                        noStimRegResults(2).spikes(:,2)<=t, :); 
  FRregionTwo(t)=length(spikes)/(t-0)
end


for t=1:500
   spikes = noStimRegResults(1).spikes(noStimRegResults(1).spikes(:,2)>=0 & ...
                        noStimRegResults(1).spikes(:,2)<=t, :); 
   FRregionOne(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = noStimRegResults(3).spikes(noStimRegResults(3).spikes(:,2)>=0 & ...
                        noStimRegResults(3).spikes(:,2)<=t, :); 
   FRregionThree(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = noStimRegResults(4).spikes(noStimRegResults(4).spikes(:,2)>=0 & ...
                        noStimRegResults(4).spikes(:,2)<=t, :); 
   FRregionFour(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = noStimRegResults(5).spikes(noStimRegResults(5).spikes(:,2)>=0 & ...
                        noStimRegResults(5).spikes(:,2)<=t, :); 
   FRregionFive(t)=length(spikes)/(t-0)
end

plot(FRregionOne,'LineWidth',2)
hold on;
plot(FRregionTwo,'LineWidth',2);
plot(FRregionThree,'LineWidth',2);
plot(FRregionFour,'LineWidth',2);
plot(FRregionFive,'LineWidth',2);
hold off
title('Firing Rate Comparison in five regions without stimulation three layer model', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
ylabel('Firing Rate (Hz)', 'FontSize', 14)
legend({'Region 1','Region 2','Region 3','Region 4','Region 5'},'Location','northwest')

figure
 total = getGroupWeights(noStimRegResults(5).params,(Br5time2weights-Br5time1weights));
 gnames = {'L3EX', 'L3IN', 'L4EX', 'L4IN', 'L5EX', 'L5IN'};
 pcolor(total);
set(gca,'YTickLabel',gnames)
set(gca,'YTick',[1:6]+0.5)
set(gca,'YTickLabel',gnames)
set(gca,'XTick',[1:6]+0.5)
set(gca,'XTickLabel',gnames)
set(gca,'XAxisLocation','Top')
xtickangle(45)
xlabel('Presynaptic Group','FontSize', 14)
ylabel('Postsynaptic Group','FontSize', 14)
c = colorbar;
%colormap (mycmap)
lims = [min(min(total)), max(max(total))];
set(gca,'clim',lims);
ylabel(c,'Change in Synaptic Weight(nS)in Region 5','FontSize', 14)
axis image
axis ij

%plot out the electrode potential in different regions
figure
plot(noStimRegResults(1).LFP(17,:)','b', 'LineWidth', 2)
hold on;
plot(noStimRegResults(2).LFP(17,:)','r', 'LineWidth', 2)
plot(noStimRegResults(3).LFP(17,:)','c', 'LineWidth', 2)
plot(noStimRegResults(4).LFP(17,:)','m', 'LineWidth', 2)
plot(noStimRegResults(5).LFP(17,:)','g', 'LineWidth', 2)
set(gcf,'color','w');
set(gca,'FontSize',16)
title(' Top left electrode No17 LFP Comparison between five regions', 'FontSize', 16)
xlabel('Time (ms)', 'FontSize', 16)
ylabel('LFP (mV)', 'FontSize', 16)
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')