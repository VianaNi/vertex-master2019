multiRegion_multilayer_stdp;  % call the setup script

% set a save directory for the simulation results. 
%%% NB!! Modify this to fit your computer's file system %%%
RecordingSettings.saveDir = '~/Documents/MATLAB/Vertex_Results/VERTEX_results_multiregion/mr4reg_diamond_nodelay_stdp_nostim_05startweight_3kconnects';

RecordingSettings.LFP = true;  % we do want to record the local field potential
RecordingSettings.weights_preN_IDs = [1:966:3864,3865:171:4546,4547:50:5000]; % set connection weights to record
[meaX, meaY, meaZ] = meshgrid(0:500:2000, 200, 700:-100:0);  % make a grid of electrode positions
RecordingSettings.meaXpositions = meaX; % set the grid positions into the VERTEX readable structure
RecordingSettings.meaYpositions = meaY;
RecordingSettings.meaZpositions = meaZ;
RecordingSettings.minDistToElectrodeTip = 20; 
RecordingSettings.v_m = 1000:1000:20000; % specify neuron IDs to record from directly
RecordingSettings.maxRecTime = 500; 
RecordingSettings.sampleRate = 1000; 
SimulationSettings.simulationTime = 500; % simulation time! Modify as you like.
SimulationSettings.timeStep = 0.03125; 
SimulationSettings.parallelSim = false; % run in parallel if you like, not neccessarily worth it for short simulations.
%SimultationSetings.onTopsy = true; % if running on a HPC use this option.
%Otherwise ignore.
% specify time points to take snapshots of the entire network weightings
RecordingSettings.weights_arr = [1 (SimulationSettings.simulationTime/SimulationSettings.timeStep)-1]; % simulation steps


%%% optional - step current stimulation to neurons to see spread of activity
%%% through region to region connections

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

[TissueParams.StimulationField, TissueParams.StimulationModel] = ...
invitroSliceStim('farapartlectrodes.stl',10000);
TissueParams.StimulationOn = [100:20:480]; 
TissueParams.StimulationOff = [105:20:485];

[params2, connections2, electrodes2] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

clear TissueParams;

 %% reset the third region
multiRegion_multilayer_stdp;

[params3, connections3, electrodes3] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);

%% reset the forth region
%multiRegion_multilayer_stdp;

[params4, connections4, electrodes4] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams,...
              RecordingSettings, SimulationSettings);
          

%% reset the fifth region
%multiRegion_multilayer_stdp;

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

 regionConnect.exportingNeuronPops{1} = 1; 
 regionConnect.exportingNeuronPops{2} = 1;
 regionConnect.exportingNeuronPops{3} = 1;
 regionConnect.exportingNeuronPops{4} = 1;
 regionConnect.exportingNeuronPops{5} = 1;
 

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
oneStimThreeLayerResults = loadResultsMultiregions(RecordingSettings.saveDir);

Region1A=mean(oneStimThreeLayerResults(1).LFP);
Region2A=mean(oneStimThreeLayerResults(2).LFP);
Region3A=mean(oneStimThreeLayerResults(3).LFP);
Region4A=mean(oneStimThreeLayerResults(4).LFP);
Region5A=mean(oneStimThreeLayerResults(5).LFP);

% % plot the mean LFP for both regions, and the difference between them
%  figure
%  subplot(311)
%  plot(mean(oneStimThreeLayerResults(1).LFP))
%  title('Region 1 averaged LFP')
%  subplot(312)
%  plot(mean(oneStimThreeLayerResults(2).LFP))
%  title('Region 2 averaged LFP')
%  subplot(313) 
%  plot(mean(oneStimThreeLayerResults(1).LFP) - mean(oneStimThreeLayerResults(2).LFP))
%  title('Difference in averaged LFP (region 1 and 2)')
%  
%  figure
%  subplot(411)
%  plot(mean(oneStimThreeLayerResults(2).LFP))
%  title('Region 2 averaged LFP')
%  subplot(412) 
%  plot(mean(oneStimThreeLayerResults(3).LFP))
%  title('Region 3 averaged LFP')
%  subplot(413) 
%  plot(mean(oneStimThreeLayerResults(2).LFP) - mean(oneStimThreeLayerResults(3).LFP))
%  title('Difference in averaged LFP (region 2 and 3)')
%  
% figure
%  subplot(611)
%  plot(mean(oneStimThreeLayerResults(1).LFP))
%  title('Region 1 averaged LFP')
%  subplot(612) 
%  plot(mean(oneStimThreeLayerResults(4).LFP))
%  title('Region 4 averaged LFP')
%  subplot(613) 
%  plot(mean(oneStimThreeLayerResults(1).LFP) - mean(oneStimThreeLayerResults(4).LFP))
%  title('Difference in averaged LFP (region 1 and 4)')
%  
% figure
%  subplot(711)
%  plot(mean(oneStimThreeLayerResults(1).LFP))
%  title('Region 1 averaged LFP')
%  subplot(712) 
%  plot(mean(oneStimThreeLayerResults(4).LFP))
%  title('Region 5 averaged LFP')
%  subplot(713) 
%  plot(mean(oneStimThreeLayerResults(1).LFP) - mean(oneStimThreeLayerResults(5).LFP))
%  title('Difference in averaged LFP (region 1 and 5)')
 

% % get the weights for the whole network at the first and last time
% % snapshots in a plottable form:
 r1time1weights=getSparseConnectivityWeights(oneStimThreeLayerResults(1).weights_arr{1},oneStimThreeLayerResults(1).syn_arr,oneStimThreeLayerResults(1).params.TissueParams.N);
 r1time2weights=getSparseConnectivityWeights(oneStimThreeLayerResults(1).weights_arr{2},oneStimThreeLayerResults(1).syn_arr,oneStimThreeLayerResults(1).params.TissueParams.N);
 r2time1weights=getSparseConnectivityWeights(oneStimThreeLayerResults(2).weights_arr{1},oneStimThreeLayerResults(2).syn_arr,oneStimThreeLayerResults(2).params.TissueParams.N);
 r2time2weights=getSparseConnectivityWeights(oneStimThreeLayerResults(2).weights_arr{2},oneStimThreeLayerResults(2).syn_arr,oneStimThreeLayerResults(2).params.TissueParams.N);
 r3time1weights=getSparseConnectivityWeights(oneStimThreeLayerResults(3).weights_arr{1},oneStimThreeLayerResults(3).syn_arr,oneStimThreeLayerResults(3).params.TissueParams.N);
 r3time2weights=getSparseConnectivityWeights(oneStimThreeLayerResults(3).weights_arr{2},oneStimThreeLayerResults(3).syn_arr,oneStimThreeLayerResults(3).params.TissueParams.N);
 r4time1weights=getSparseConnectivityWeights(oneStimThreeLayerResults(4).weights_arr{1},oneStimThreeLayerResults(4).syn_arr,oneStimThreeLayerResults(4).params.TissueParams.N);
 r4time2weights=getSparseConnectivityWeights(oneStimThreeLayerResults(4).weights_arr{2},oneStimThreeLayerResults(4).syn_arr,oneStimThreeLayerResults(4).params.TissueParams.N);
 r5time1weights=getSparseConnectivityWeights(oneStimThreeLayerResults(5).weights_arr{1},oneStimThreeLayerResults(5).syn_arr,oneStimThreeLayerResults(5).params.TissueParams.N);
 r5time2weights=getSparseConnectivityWeights(oneStimThreeLayerResults(5).weights_arr{2},oneStimThreeLayerResults(5).syn_arr,oneStimThreeLayerResults(5).params.TissueParams.N);
% % 
% % % plot the weight differences between the start and end of the simulation
% % % for each network. Colours represent the weight changes.
figure
imagesc(log(abs(r1time2weights-r1time1weights)));
title('Weight changes between different neurons in Region1');
colorbar;


figure
imagesc(log(abs(r2time2weights-r2time1weights)));
title('Weight changes between different neurons in Region2');
colorbar;


figure
imagesc(log(abs(r3time2weights-r3time1weights)));
title('Weight changes between different neurons in Region3');
colorbar;



figure
imagesc(log(abs((r4time2weights-r4time1weights))));
title('Weight changes between different neurons in Region4');
colorbar;


figure
imagesc(log(abs((r5time2weights-r5time1weights))));
title('Weight changes between different neurons in Region5');
colorbar;




rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 1 One Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 511;
plotSpikeRaster(oneStimThreeLayerResults(1), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 2 One Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 512;
plotSpikeRaster(oneStimThreeLayerResults(2), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 3 One Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 513;
plotSpikeRaster(oneStimThreeLayerResults(3), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster Region 4 One Stimulation';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 514;
plotSpikeRaster(oneStimThreeLayerResults(4), rasterParams);

rasterParams.colors = {'k','m','k','m','k','m'};
rasterParams.groupBoundaryLines = 'c';
rasterParams.title = 'Spike Raster 5';
rasterParams.xlabel = 'Time (ms)';
rasterParams.ylabel = 'Neuron ID';
rasterParams.figureID = 515;
plotSpikeRaster(oneStimThreeLayerResults(5), rasterParams);

% figure
plot(oneStimThreeLayerResults(1).LFP(1,:)','b', 'LineWidth', 2)
hold on;
plot(oneStimThreeLayerResults(2).LFP(1,:)','r', 'LineWidth', 2)
plot(oneStimThreeLayerResults(3).LFP(1,:)','c', 'LineWidth', 2)
plot(oneStimThreeLayerResults(4).LFP(1,:)','m', 'LineWidth', 2)
plot(oneStimThreeLayerResults(5).LFP(1,:)','g', 'LineWidth', 2)
set(gcf,'color','w');
set(gca,'FontSize',12)
title(' Top left electrode No1 LFP Comparison between five regions', 'FontSize', 12)
xlabel('Time (ms)', 'FontSize', 12)
ylabel('LFP (mV)', 'FontSize', 12)
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northwest')


figure
plot(oneStimThreeLayerResults(1).LFP(17,:)','b', 'LineWidth', 2)
hold on;
plot(oneStimThreeLayerResults(2).LFP(17,:)','r', 'LineWidth', 2)
plot(oneStimThreeLayerResults(3).LFP(17,:)','c', 'LineWidth', 2)
plot(oneStimThreeLayerResults(4).LFP(17,:)','m', 'LineWidth', 2)
plot(oneStimThreeLayerResults(5).LFP(17,:)','g', 'LineWidth', 2)
set(gcf,'color','w');
set(gca,'FontSize',16)
title(' Top left electrode No17 LFP Comparison between five regions', 'FontSize', 16)
xlabel('Time (ms)', 'FontSize', 16)
ylabel('LFP (mV)', 'FontSize', 16)
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')


figure
plot(oneStimThreeLayerResults(1).LFP(24,:)','b', 'LineWidth', 2)
hold on;
plot(oneStimThreeLayerResults(2).LFP(24,:)','r', 'LineWidth', 2)
plot(oneStimThreeLayerResults(3).LFP(24,:)','c', 'LineWidth', 2)
plot(oneStimThreeLayerResults(4).LFP(24,:)','m', 'LineWidth', 2)
plot(oneStimThreeLayerResults(5).LFP(24,:)','g', 'LineWidth', 2)
set(gcf,'color','w');
set(gca,'FontSize',16)
title(' Top left electrode No24 LFP Comparison between five regions', 'FontSize', 16)
xlabel('Time (ms)', 'FontSize', 16)
ylabel('LFP (mV)', 'FontSize', 16)
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')


figure
plot(oneStimThreeLayerResults(1).LFP(40,:)','b', 'LineWidth', 2)
hold on;
plot(oneStimThreeLayerResults(2).LFP(40,:)','r', 'LineWidth', 2)
plot(oneStimThreeLayerResults(3).LFP(40,:)','c', 'LineWidth', 2)
plot(oneStimThreeLayerResults(4).LFP(40,:)','m', 'LineWidth', 2)
plot(oneStimThreeLayerResults(5).LFP(40,:)','g', 'LineWidth', 2)
set(gcf,'color','w');
set(gca,'FontSize',16)
title(' Top left electrode No40 LFP Comparison between five regions', 'FontSize', 16)
xlabel('Time (ms)', 'FontSize', 16)
ylabel('LFP (mV)', 'FontSize', 16)
legend({'Region1','Region2','Region3','Region4','Region5'},'Location','northeast')



 figure
plot(mean(oneStimThreeLayerResults(1).LFP(:,1:490)),'Linewidth',2)
 title ('Five Region average LFP Comparison One Stimulation')
  hold on
plot(mean(oneStimThreeLayerResults(2).LFP(:,1:490)),'Linewidth',2)
plot(mean(oneStimThreeLayerResults(3).LFP(:,1:490)),'Linewidth',2)
plot(mean(oneStimThreeLayerResults(4).LFP(:,1:490)),'Linewidth',2)
plot(mean(oneStimThreeLayerResults(5).LFP(:,1:490)),'Linewidth',2)
 legend({'Region1','Region2','Region3','Region4','Region5'},'Location','southwest')
 hold off
 
lmax=mean(oneStimThreeLayerResults(1).LFP(1,:));
imax=1;
for i=1:40
    if mean(oneStimThreeLayerResults(1).LFP(i,:)) > lmax
        imax=i;
        lmax=mean(oneStimThreeLayerResults(1).LFP(i,:));
    end
end
   
        
lmin=mean(oneStimThreeLayerResults(1).LFP(1,:));
imin=1;
for i=1:40
    if mean(oneStimThreeLayerResults(1).LFP(i,:)) < lmin
        imin=i;
        lmin=mean(oneStimThreeLayerResults(1).LFP(i,:));
    end
end




%% Plot weight change by group
 
 figure
 total = getGroupWeights(oneStimThreeLayerResults(3).params,(r3time2weights-r3time1weights));
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
ylabel(c,'Change in Synaptic Weight(nS)in Region 3','FontSize', 14)
axis image
axis ij

figure
 total = getGroupWeights(oneStimThreeLayerResults(2).params,(r2time2weights-r2time1weights));
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
lims = [min(min(total)), max(max(total))];
set(gca,'clim',lims);
ylabel(c,'Change in Synaptic Weight (nS) in Region 2','FontSize', 14)
axis image
axis ij

figure
 total = getGroupWeights(oneStimThreeLayerResults(1).params,(r1time2weights-r1time1weights));
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
ylabel(c,'Change in Synaptic Weight in region 1','FontSize', 14)
axis image
axis ij
 
%% Plot weight change spatially


pars.toPlot = 1:1:oneStimThreeLayerResults(4).params.TissueParams.N;
plotSomaPositionsMembranePotential(oneStimThreeLayerResults(4).params.TissueParams,pars,sum(r4time2weights-r4time1weights));
colorbar;
title('Region 4 spatial weight change view');


pars.toPlot = 1:1:oneStimThreeLayerResults(3).params.TissueParams.N;
plotSomaPositionsMembranePotential(oneStimThreeLayerResults(3).params.TissueParams,pars,sum(r3time2weights-r3time1weights));
colorbar;
title('Region 3 spatial weight change view');


pars.toPlot = 1:1:oneStimThreeLayerResults(2).params.TissueParams.N;
plotSomaPositionsMembranePotential(oneStimThreeLayerResults(2).params.TissueParams,pars,sum(r2time2weights-r2time1weights));
colorbar;
title('Region 2 spatial weight change view');


pars.toPlot = 1:1:oneStimThreeLayerResults(5).params.TissueParams.N;
plotSomaPositionsMembranePotential(oneStimThreeLayerResults(5).params.TissueParams,pars,sum(r1time2weights-r1time1weights));
colorbar;
title('Region 5 spatial weight change view');


% ElecArray=[1,16,27,40];
% for j=1:ElecArray
%      for i=1:5
%      Electrode{j}=[oneStimThreeLayerResults(i).LFP(1,150:450)];
%      var(Electrode{j})
%      end
% end
% 

% figure;
% pdeplot3D(oneStimThreeLayerResults(2).params.TissueParams.StimulationModel, 'ColorMap', oneStimThreeLayerResults(2).params.TissueParams.StimulationField.NodalSolution, 'FaceAlpha', 0.8)
% hold on;
% 
% pars.colors = rasterParams.colors;
% 
% plotSomaPositions(oneStimThreeLayerResults(2).params.TissueParams,pars);



figure
plot(mean(oneStimThreeLayerResults(2).LFP(:,1:490)),'Linewidth',2);
hold on; 
for i=100:20:480
 plot((i:i+5),Region2A(i:i+5),'Color','r','Linewidth',2)
end
hold off
title('local field potential in Region 2', 'FontSize', 12)
xlabel('Time (ms)', 'FontSize', 12)
ylabel('LFP (mV)', 'FontSize', 12)
legend({'No stimulation','Stimulation'},'Location','southwest')

for t=1:500
   spikes = oneStimThreeLayerResults(2).spikes(oneStimThreeLayerResults(2).spikes(:,2)>=0 & ...
                        oneStimThreeLayerResults(2).spikes(:,2)<=t, :); 
  FRregionTwo(t)=length(spikes)/(t-0)
end


for t=1:500
   spikes = oneStimThreeLayerResults(1).spikes(oneStimThreeLayerResults(1).spikes(:,2)>=0 & ...
                        oneStimThreeLayerResults(1).spikes(:,2)<=t, :); 
   FRregionOne(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = oneStimThreeLayerResults(3).spikes(oneStimThreeLayerResults(3).spikes(:,2)>=0 & ...
                        oneStimThreeLayerResults(3).spikes(:,2)<=t, :); 
   FRregionThree(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = oneStimThreeLayerResults(4).spikes(oneStimThreeLayerResults(4).spikes(:,2)>=0 & ...
                        oneStimThreeLayerResults(4).spikes(:,2)<=t, :); 
   FRregionFour(t)=length(spikes)/(t-0)
end

for t=1:500
   spikes = oneStimThreeLayerResults(5).spikes(oneStimThreeLayerResults(5).spikes(:,2)>=0 & ...
                        oneStimThreeLayerResults(5).spikes(:,2)<=t, :); 
   FRregionFive(t)=length(spikes)/(t-0)
end




plot(FRregionOne,'LineWidth',2)
hold on;
plot(FRregionTwo,'LineWidth',2);
plot(FRregionThree,'LineWidth',2);
plot(FRregionFour,'LineWidth',2);
plot(FRregionFive,'LineWidth',2);
hold off
title('Firing Rate Comparison in five regions with one stimulation three layer model', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
ylabel('Firing Rate (Hz)', 'FontSize', 14)
legend({'Region 1','Region 2','Region 3','Region 4','Region 5'},'Location','northwest')

groupRates(oneStimThreeLayerResults(1), 0, 500)
groupRates(oneStimThreeLayerResults(2), 0, 500)
groupRates(oneStimThreeLayerResults(3), 0, 500)
groupRates(oneStimThreeLayerResults(4), 0, 500)
groupRates(oneStimThreeLayerResults(5), 0, 500)
