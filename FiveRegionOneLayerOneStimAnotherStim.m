setup_multiregion_allstdp;  



RecordingSettings.saveDir = '~/Documents/MATLAB/Vertex_Results/VERTEX_results_multiregion/mr4reg_diamond_nodelay_stdp_nostim_05startweight_3kconnects';
% we do want to record the local field potential
RecordingSettings.weights_preN_IDs = [1:966:3864,3865:171:4546,4547:50:5000]; % set connection weights to record
[meaX, meaY, meaZ] = meshgrid(0:100:2000, 200, 195:-95:5);  % make a grid of electrode positions
RecordingSettings.meaXpositions = meaX; % set the grid positions into the VERTEX readable structure
RecordingSettings.meaYpositions = meaY;
RecordingSettings.meaZpositions = meaZ;
RecordingSettings.minDistToElectrodeTip = 20; 
RecordingSettings.v_m = 100:200:1200; % specify neuron IDs to record from directly
RecordingSettings.maxRecTime = 500; 
RecordingSettings.sampleRate = 1000; 
RecordingSettings.stdpvars = [1:2:500];
% Recording a snapshot of the weights of the entire network at the
% specified timestep. 
%RecordingSettings.weights_arr = [1000:1000:6000];
SimulationSettings.simulationTime = 500; % simulation time! Modify as you like.
SimulationSettings.timeStep = 0.03125; 
SimulationSettings.parallelSim = false; % run in parallel if you like, not neccessarily worth it for short simulations.
%SimultationSetings.onTopsy = true; % if running on a HPC use this option.
%Otherwise ignore.

% specify time points to take snapshots of the entire network weightings
RecordingSettings.weights_arr = [1 (SimulationSettings.simulationTime/SimulationSettings.timeStep)-1]; % simulation steps


%%% optional - step current stimulation to neurons to see spread of activity
%%% through region to region connections

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
[TissueParams.StimulationField, TissueParams.model] = ...
    invitroSliceStim('farapartlectrodesbig.stl',10000);
TissueParams.StimulationOn = [100:20:400];
TissueParams.StimulationOff = [105:20:405];

[params3, connections3, electrodes3] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);
 clear TissueParams;

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
oneStimMultiRegResultsAno = loadResultsMultiregions(RecordingSettings.saveDir);

Region1D=mean(oneStimMultiRegResultsAno(1).LFP);
Region2D=mean(oneStimMultiRegResultsAno(2).LFP);
Region3D=mean(oneStimMultiRegResultsAno(3).LFP);
Region4D=mean(oneStimMultiRegResultsAno(4).LFP);
Region5D=mean(oneStimMultiRegResultsAno(5).LFP);

% plot the spike raster for each region
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 1 Spike Raster One Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(oneStimMultiRegResultsAno(1),rasterParams)
 
rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 2 Spike Raster One Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(oneStimMultiRegResultsAno(2),rasterParams)
 
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 3 Spike Raster One Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(oneStimMultiRegResultsAno(3),rasterParams)
 
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 4 Spike Raster One Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(oneStimMultiRegResultsAno(4),rasterParams)
 
 rasterParams.colors = {'k', 'm','k'};
 rasterParams.groupBoundaryLines = [0.7, 0.7, 0.7];
 rasterParams.title='Region 5 Spike Raster One Stimulation';
 rasterParams.xlabel='time (ms)';
 rasterParams.ylabel='Neuron ID';
 plotSpikeRaster(oneStimMultiRegResultsAno(5),rasterParams)
 

 
% % plot the mean LFP for both regions, and the difference between them
 figure
 subplot(311)
 plot(mean(oneStimMultiRegResultsAno(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(312)
 plot(mean(oneStimMultiRegResultsAno(2).LFP))
 title('Region 2 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(313)
 plot(mean(oneStimMultiRegResultsAno(1).LFP) - mean(oneStimMultiRegResultsAno(2).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 

 figure
 subplot(411)
 plot(mean(oneStimMultiRegResultsAno(2).LFP))
 title('Region 2 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(412)
 plot(mean(oneStimMultiRegResultsAno(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(413)
 plot(mean(oneStimMultiRegResultsAno(2).LFP) - mean(oneStimMultiRegResultsAno(3).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 

 figure
 subplot(911)
 plot(mean(oneStimMultiRegResultsAno(1).LFP))
 title('Region 1 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(912)
 plot(mean(oneStimMultiRegResultsAno(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(913)
 plot(mean(oneStimMultiRegResultsAno(1).LFP) - mean(oneStimMultiRegResultsAno(3).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 
 
 figure
 subplot(511)
 plot(mean(oneStimMultiRegResultsAno(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(512)
 plot(mean(oneStimMultiRegResultsAno(4).LFP))
 title('Region 4 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(513)
 plot(mean(oneStimMultiRegResultsAno(3).LFP) - mean(oneStimMultiRegResultsAno(4).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');


 figure
 subplot(611)
 plot(mean(oneStimMultiRegResultsAno(3).LFP))
 title('Region 3 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(612)
 plot(mean(oneStimMultiRegResultsAno(5).LFP))
 title('Region 5 averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 subplot(613)
 plot(mean(oneStimMultiRegResultsAno(3).LFP) - mean(oneStimMultiRegResultsAno(5).LFP))
 title('Difference in averaged LFP')
 xlabel('time (ms)');
 ylabel('LFP (mv)');

figure
 plot(mean(oneStimMultiRegResultsAno(1).LFP(:,1:490)),'LineWidth',2)
 title ('Five Region average LFP Comparison')
  hold on
 plot(mean(oneStimMultiRegResultsAno(2).LFP(:,1:490)),'LineWidth',2)
 plot(mean(oneStimMultiRegResultsAno(3).LFP(:,1:490)),'LineWidth',2)
 plot(mean(oneStimMultiRegResultsAno(4).LFP(:,1:490)),'LineWidth',2)
 plot(mean(oneStimMultiRegResultsAno(1).LFP(:,1:490)),'LineWidth',2)
 legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')
 xlabel('time (ms)');
 ylabel('LFP (mv)');
 hold off
  

% get the weights for the whole network at the first and last time
% snapshots in a plottable form:
 Dr1time1weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(1).weights_arr{1},oneStimMultiRegResultsAno(1).syn_arr,oneStimMultiRegResultsAno(1).params.TissueParams.N);
 Dr1time2weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(1).weights_arr{2},oneStimMultiRegResultsAno(1).syn_arr,oneStimMultiRegResultsAno(1).params.TissueParams.N);
 Dr2time1weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(2).weights_arr{1},oneStimMultiRegResultsAno(2).syn_arr,oneStimMultiRegResultsAno(2).params.TissueParams.N);
 Dr2time2weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(2).weights_arr{2},oneStimMultiRegResultsAno(2).syn_arr,oneStimMultiRegResultsAno(2).params.TissueParams.N);
 Dr3time1weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(3).weights_arr{1},oneStimMultiRegResultsAno(3).syn_arr,oneStimMultiRegResultsAno(3).params.TissueParams.N);
 Dr3time2weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(3).weights_arr{2},oneStimMultiRegResultsAno(3).syn_arr,oneStimMultiRegResultsAno(3).params.TissueParams.N);
 Dr4time1weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(4).weights_arr{1},oneStimMultiRegResultsAno(4).syn_arr,oneStimMultiRegResultsAno(4).params.TissueParams.N);
 Dr4time2weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(4).weights_arr{2},oneStimMultiRegResultsAno(4).syn_arr,oneStimMultiRegResultsAno(4).params.TissueParams.N);
 Dr5time1weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(5).weights_arr{1},oneStimMultiRegResultsAno(5).syn_arr,oneStimMultiRegResultsAno(5).params.TissueParams.N);
 Dr5time2weights=getSparseConnectivityWeights(oneStimMultiRegResultsAno(5).weights_arr{2},oneStimMultiRegResultsAno(5).syn_arr,oneStimMultiRegResultsAno(5).params.TissueParams.N);
% plot the weight differences between the start and end of the simulation
% for each network. Colours represent the weight changes.

figure
imagesc(log(abs(Dr1time2weights-Dr1time1weights)));
title('Weight changes between different neurons in Region1');
colorbar;
caxis([-10 6])


figure
imagesc(log(abs(Dr2time2weights-Dr2time1weights)));
title('Weight changes between different neurons in Region2');
colorbar;
caxis([-10 6])

figure
imagesc(log(abs(Dr3time2weights-Dr3time1weights)));
title('Weight changes between different neurons in Region3');
colorbar;
caxis([-10 6])


figure
imagesc(log(abs(Dr4time2weights-Dr4time1weights)));
title('Weight changes between different neurons in Region4');
colorbar;
caxis([-10 6])


figure
imagesc(log(abs(Dr5time2weights-Dr5time1weights)));
title('Weight changes between different neurons in Region5');
colorbar;
caxis([-10 6])



groupRates(oneStimMultiRegResultsAno(1), 0, 500)
groupRates(oneStimMultiRegResultsAno(2), 0, 500)
groupRates(oneStimMultiRegResultsAno(3), 0, 500)
groupRates(oneStimMultiRegResultsAno(4), 0, 500)
groupRates(oneStimMultiRegResultsAno(5), 0, 500)

figure
plot(mean(oneStimMultiRegResultsAno(3).LFP(:,1:490)),'LineWidth',2);
hold on; 
for i=100:20:400
 plot((i:i+5),Region3D(i:i+5),'Color','r','LineWidth',2)
end
hold off
title('local field potential in Region 3', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
ylabel('LFP (mV)', 'FontSize', 14)
legend({'No stimulation','Stimulation'},'Location','northwest')



for t=1:500
   spikes = oneStimMultiRegResultsAno(2).spikes(oneStimMultiRegResultsAno(2).spikes(:,2)>=0 & ...
                        oneStimMultiRegResultsAno(2).spikes(:,2)<=t, :); 
   FRregionTwo(t)=length(spikes)/(t-0)
end


for t=1:500
   spikes = oneStimMultiRegResultsAno(1).spikes(oneStimMultiRegResultsAno(1).spikes(:,2)>=0 & ...
                        oneStimMultiRegResultsAno(1).spikes(:,2)<=t, :); 
   FRregionOne(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = oneStimMultiRegResultsAno(3).spikes(oneStimMultiRegResultsAno(3).spikes(:,2)>=0 & ...
                        oneStimMultiRegResultsAno(3).spikes(:,2)<=t, :); 
   FRregionThree(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = oneStimMultiRegResultsAno(4).spikes(oneStimMultiRegResultsAno(4).spikes(:,2)>=0 & ...
                        oneStimMultiRegResultsAno(4).spikes(:,2)<=t, :); 
   FRregionFour(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = oneStimMultiRegResultsAno(5).spikes(oneStimMultiRegResultsAno(5).spikes(:,2)>=0 & ...
                        oneStimMultiRegResultsAno(5).spikes(:,2)<=t, :); 
   FRregionFive(t)=length(spikes)/(t-0)
end

plot(FRregionOne,'LineWidth',2)
hold on;
plot(FRregionTwo,'LineWidth',2);
plot(FRregionThree,'LineWidth',2);
plot(FRregionFour,'LineWidth',2);
plot(FRregionFive,'LineWidth',2);
hold off
title('Firing Rate Comparison in five regions', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
ylabel('Firing Rate (Hz)', 'FontSize', 14)
legend({'Region 1','Region x2','Region 3','Region 4','Region 5'},'Location','northwest')

%correlation
% for t=1:500
%    spikes = oneStimMultiRegResultsAno(5).spikes(oneStimMultiRegResultsAno(5).spikes(:,2)>=0 & ...
%                         oneStimMultiRegResultsAno(5).spikes(:,2)<=t, :); 
%     relationmatrix=corr(mean(oneStimMultiRegResultsAno(2).LFP), length(spikes));
% end
% imagesc(relationmatrix)


figure
 total = getGroupWeights(oneStimMultiRegResultsAno(3).params,(Dr3time2weights-Dr3time1weights));
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
ylabel(c,'Change in Synaptic Weight After TBS (nS)','FontSize', 14)
axis image
axis ij