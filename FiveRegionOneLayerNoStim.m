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
          
%% reset the second region          
setup_multiregion_allstdp;
[params2, connections2, electrodes2] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);
       

 %% reset the third region
% setup_multiregion_gexp;
% setup_multiregion_withinboundconnection;

setup_multiregion_allstdp;

[params3, connections3, electrodes3] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

 %% reset the forth region
% setup_multiregion_gexp;
% setup_multiregion_withinboundconnection;
setup_multiregion_allstdp;

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
 % for example [1,1;0,1] there are two regions and there is only an
 % external connection from region 1 to region 2, it is not returned, and
 % while they do connect to themselves internally for the sake of incoming external
 % connections the diagonals are set to 0.
 % NB: Values represent delays, so a binary connectivity map means no delay is
 % present. 
 
 % Identify the neuron types (e.g. NP(1) in this instance are the
 % excitatory PY cells) which will export signals. Use [] if not exporting.
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
noStimMultiRegReults = loadResultsMultiregions(RecordingSettings.saveDir);

Region1C=mean(noStimMultiRegReults(1).LFP);
Region2C=mean(noStimMultiRegReults(2).LFP);
Region3C=mean(noStimMultiRegReults(3).LFP);
Region4C=mean(noStimMultiRegReults(4).LFP);
Region5C=mean(noStimMultiRegReults(5).LFP);

% plot the spike raster for each region
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 1 Spike Raster No Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(noStimMultiRegReults(1),rasterParams)
 
rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 2 Spike Raster No Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(noStimMultiRegReults(2),rasterParams)
 
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 3 Spike Raster No Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(noStimMultiRegReults(3),rasterParams)
 
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 4 Spike Raster No Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(noStimMultiRegReults(4),rasterParams)
 
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 5 Spike Raster No Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(noStimMultiRegReults(5),rasterParams)
 


 
% % plot the mean LFP for both regions, and the difference between them
 figure
 subplot(311)
 plot(mean(noStimMultiRegReults(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(312)
 plot(mean(noStimMultiRegReults(2).LFP))
 title('Region 2 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(313)
 plot(mean(noStimMultiRegReults(1).LFP) - mean(noStimMultiRegReults(2).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 

 figure
 subplot(411)
 plot(mean(noStimMultiRegReults(2).LFP))
 title('Region 2 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(412)
 plot(mean(noStimMultiRegReults(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(413)
 plot(mean(noStimMultiRegReults(2).LFP) - mean(noStimMultiRegReults(3).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 

 figure
 subplot(911)
 plot(mean(noStimMultiRegReults(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(912)
 plot(mean(noStimMultiRegReults(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(913)
 plot(mean(noStimMultiRegReults(1).LFP) - mean(noStimMultiRegReults(3).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 
 
 figure
 subplot(511)
 plot(mean(noStimMultiRegReults(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(512)
 plot(mean(noStimMultiRegReults(4).LFP))
 title('Region 4 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(513)
 plot(mean(noStimMultiRegReults(3).LFP) - mean(noStimMultiRegReults(4).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');


 figure
 subplot(611)
 plot(mean(noStimMultiRegReults(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(612)
 plot(mean(noStimMultiRegReults(5).LFP))
 title('Region 5 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(613)
 plot(mean(noStimMultiRegReults(3).LFP) - mean(noStimMultiRegReults(5).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');

 %plot out the mean LFP of each region in one plot
figure
 plot(mean(noStimMultiRegReults(1).LFP(:,1:490)),'LineWidth',2)
 title ('Five Region average LFP Comparison')
  hold on
 plot(mean(noStimMultiRegReults(2).LFP(:,1:490)),'LineWidth',2)
 plot(mean(noStimMultiRegReults(3).LFP(:,1:490)),'LineWidth',2)
 plot(mean(noStimMultiRegReults(4).LFP(:,1:490)),'LineWidth',2)
 plot(mean(noStimMultiRegReults(5).LFP(:,1:490)),'LineWidth',2)
 legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 hold off
  

% get the weights for the whole network at the first and last time
% snapshots in a plottable form:
 Cr1time1weights=getSparseConnectivityWeights(noStimMultiRegReults(1).weights_arr{1},noStimMultiRegReults(1).syn_arr,noStimMultiRegReults(1).params.TissueParams.N);
 Cr1time2weights=getSparseConnectivityWeights(noStimMultiRegReults(1).weights_arr{2},noStimMultiRegReults(1).syn_arr,noStimMultiRegReults(1).params.TissueParams.N);
 Cr2time1weights=getSparseConnectivityWeights(noStimMultiRegReults(2).weights_arr{1},noStimMultiRegReults(2).syn_arr,noStimMultiRegReults(2).params.TissueParams.N);
 Cr2time2weights=getSparseConnectivityWeights(noStimMultiRegReults(2).weights_arr{2},noStimMultiRegReults(2).syn_arr,noStimMultiRegReults(2).params.TissueParams.N);
 Cr3time1weights=getSparseConnectivityWeights(noStimMultiRegReults(3).weights_arr{1},noStimMultiRegReults(3).syn_arr,noStimMultiRegReults(3).params.TissueParams.N);
 Cr3time2weights=getSparseConnectivityWeights(noStimMultiRegReults(3).weights_arr{2},noStimMultiRegReults(3).syn_arr,noStimMultiRegReults(3).params.TissueParams.N);
 Cr4time1weights=getSparseConnectivityWeights(noStimMultiRegReults(4).weights_arr{1},noStimMultiRegReults(4).syn_arr,noStimMultiRegReults(4).params.TissueParams.N);
 Cr4time2weights=getSparseConnectivityWeights(noStimMultiRegReults(4).weights_arr{2},noStimMultiRegReults(4).syn_arr,noStimMultiRegReults(4).params.TissueParams.N);
 Cr5time1weights=getSparseConnectivityWeights(noStimMultiRegReults(5).weights_arr{1},noStimMultiRegReults(5).syn_arr,noStimMultiRegReults(5).params.TissueParams.N);
 Cr5time2weights=getSparseConnectivityWeights(noStimMultiRegReults(5).weights_arr{2},noStimMultiRegReults(5).syn_arr,noStimMultiRegReults(5).params.TissueParams.N);

 % plot the weight differences between the start and end of the simulation
% for each network. Colours represent the weight changes. Take the log
% scale to see the results clearly

figure
imagesc(log(abs(Cr1time2weights-Cr1time1weights)));
title('Weight changes between different neurons in Region1');
colorbar;
caxis([-12 6])


figure
imagesc(log(abs(Cr2time2weights-Cr2time1weights)));
title('Weight changes between different neurons in Region2');
colorbar;
caxis([-12 6])

figure
imagesc(log(abs(Cr3time2weights-Cr3time1weights)));
title('Weight changes between different neurons in Region3');
colorbar;
caxis([-12 6])


figure
imagesc(log(abs(Cr4time2weights-Cr4time1weights)));
title('Weight changes between different neurons in Region4');
colorbar;
caxis([-12 6])


figure
imagesc(log(abs(Cr5time2weights-Cr5time1weights)));
title('Weight changes between different neurons in Region5');
colorbar;
caxis([-12 6])



%calculate group rates of each region
groupRates(noStimMultiRegReults(1), 0, 500)
groupRates(noStimMultiRegReults(2), 0, 500)
groupRates(noStimMultiRegReults(3), 0, 500)
groupRates(noStimMultiRegReults(4), 0, 500)
groupRates(noStimMultiRegReults(5), 0, 500)

figure
plot(mean(noStimMultiRegReults(3).LFP(:,1:490)),'LineWidth',2);
title('local field potential in Region 3', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
ylabel('LFP (mV)', 'FontSize', 14)
legend({'No stimulation','Stimulation'},'Location','northwest')



%plot out the total firing rates of all neurons in each region
for t=1:500
   spikes = noStimMultiRegReults(2).spikes(noStimMultiRegReults(2).spikes(:,2)>=0 & ...
                        noStimMultiRegReults(2).spikes(:,2)<=t, :); 
  FRregionTwo(t)=length(spikes)/(t-0)
end


for t=1:500
   spikes = noStimMultiRegReults(1).spikes(noStimMultiRegReults(1).spikes(:,2)>=0 & ...
                        noStimMultiRegReults(1).spikes(:,2)<=t, :); 
   FRregionOne(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = noStimMultiRegReults(3).spikes(noStimMultiRegReults(3).spikes(:,2)>=0 & ...
                        noStimMultiRegReults(3).spikes(:,2)<=t, :); 
   FRregionThree(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = noStimMultiRegReults(4).spikes(noStimMultiRegReults(4).spikes(:,2)>=0 & ...
                        noStimMultiRegReults(4).spikes(:,2)<=t, :); 
   FRregionFour(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = noStimMultiRegReults(5).spikes(noStimMultiRegReults(5).spikes(:,2)>=0 & ...
                        noStimMultiRegReults(5).spikes(:,2)<=t, :); 
   FRregionFive(t)=length(spikes)/(t-0)
end

plot(FRregionOne,'LineWidth',2)
hold on;
plot(FRregionTwo,'LineWidth',2);
plot(FRregionThree,'LineWidth',2);
plot(FRregionFour,'LineWidth',2);
plot(FRregionFive,'LineWidth',2);
hold off
title('Firing Rate Comparison in five regions without stimulation', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
ylabel('Firing Rate (Hz)', 'FontSize', 14)
legend({'Region 1','Region 2','Region 3','Region 4','Region 5'},'Location','northwest')




%correlation
% for t=1:500
%    spikes = noStimMultiRegReults(5).spikes(noStimMultiRegReults(5).spikes(:,2)>=0 & ...
%                         noStimMultiRegReults(5).spikes(:,2)<=t, :); 
%     relationmatrix=corr(mean(noStimMultiRegReults(2).LFP), length(spikes));
% end
% imagesc(relationmatrix)


%plot out synaptic weight change
figure
 total = getGroupWeights(noStimMultiRegReults(3).params,(Dr3time2weights-Dr3time1weights));
 gnames = {'Group 1', 'Group 2', 'Group 3'};
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
ylabel(c,'Change in Synaptic Weight change (nS)','FontSize', 14)
axis image
axis ij




%Comparison within the regions
plot(Region1C(:,1:490),'LineWidth',2);
hold on;
plot(Region1D(:,1:490),'LineWidth',2);
legend({'No stim', 'One stim'},'Location','northeast');
title('Region 1 LFP Comparison');
hold off
