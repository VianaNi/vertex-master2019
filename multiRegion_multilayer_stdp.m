% Set up script for three layer regions

%% Tissue parameters
% First we specify the same tissue parameters as in tutorial 1:

stdprate=0.1;

TissueParams.X = 2200;
TissueParams.Y = 400;
TissueParams.Z = 1240;
TissueParams.neuronDensity = 20000;

TissueParams.numLayers = 3;
TissueParams.layerBoundaryArr = [1240, 450, 150, 0];
TissueParams.numStrips = 50;
TissueParams.tissueConductivity = 0.3;
TissueParams.maxZOverlap = [-1 , -1];

%% Neuron parameters
% Next we will specify the parameters for our two neuron groups. We will
% use the neuron models described in (Tomsett et al. 2014) for layer 2/3
% pyramidal neurons and basket interneurons. We are going to set 85% of the
% neurons to be pyramidal cells, in neuron group 1.clear al

NeuronParams(1).modelProportion = 0.4;
NeuronParams(1).somaLayer = 1;

%%
% We are going to use the adaptive exponential (AdEx) model to generate
% spiking dynamics (Brette & Gerstner 2005). Our models will be different
% from the original AdEx model as they will also contain dendritic
% compartments so that they can generate an extracellular potential. The
% AdEx dynamics are in the soma compartment, while the dendrites are
% passive. The AdEx model requires us to specify some extra parameters that
% control the model dynamics:

NeuronParams(1).neuronModel = 'adex';
NeuronParams(1).V_t = -50;
NeuronParams(1).delta_t = 2;
NeuronParams(1).a = 2.6;
NeuronParams(1).tau_w = 65;
NeuronParams(1).b = 220;
NeuronParams(1).v_reset = -60;
NeuronParams(1).v_cutoff = -45;

%%
% |V_t| is the spike generation threshold (in mV), |delta_t| is the spike
% steepness parameter (in mV), |a| is the scale factor of the spike
% after-hyperpolarisation (AHP) current (in nanoSiemens), |tau_w| is the
% AHP current time constant (in ms), |b| is the instantaneous change in the
% AHP current after a spike (in pA), |v_reset| is the membrane potential
% that the soma compartment is reset to after firing a spike (in mV), and
% |v_cutoff| is the membrane potential at which a spike is detected (in
% mV). We recommend that this parameter is set to |V_t + 5|; if it is set 
% much higher, the exponential term in the AdEx equations can lead to the 
% membrane potentialexploding to a not-a-number (NaN) value, which breaks 
% things.
%
% The remaining parameters defining the structure and passive properties
% are the same as in Tutorial 1:
NeuronParams(1).numCompartments = 8;
NeuronParams(1).compartmentParentArr = [0, 1, 2, 2, 4, 1, 6, 6];
NeuronParams(1).compartmentLengthArr = [13 48 124 145 137 40 143 143];
NeuronParams(1).compartmentDiameterArr = ...
  [29.8, 3.75, 1.91, 2.81, 2.69, 2.62, 1.69, 1.69];
NeuronParams(1).compartmentXPositionMat = ...
[   0,    0;
    0,    0;
    0,  124;
    0,    0;
    0,    0;
    0,    0;
    0, -139;
    0,  139];
NeuronParams(1).compartmentYPositionMat = ...
[   0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0];
NeuronParams(1).compartmentZPositionMat = ...
[ -13,    0;
    0,   48;
   48,   48;
   48,  193;
  193,  330;
  -13,  -53;
  -53, -139;
  -53, -139];
NeuronParams(1).axisAligned = 'z';
NeuronParams(1).C = 1.0*2.96;
NeuronParams(1).R_M = 20000/2.96;
NeuronParams(1).R_A = 150;
NeuronParams(1).E_leak = -70;
NeuronParams(1).somaID = 1;
NeuronParams(1).basalID = [6, 7, 8];
NeuronParams(1).apicalID = [2 3 4 5];
NeuronParams(1).labelNames = {'somaID', 'basalID','apicalID'};
NeuronParams(1).minCompartmentSize = 0.8;
%%
% In order to generate spikes, we need to provide the neurons with some
% input. We set the inputs to our neuron group in another structure array,
% called Inptut, that we create as a field in our NeuronParams(1)
% structure. Input is a structure array rather than just a structure
% so that we can specify multiple different inputs to the neurons in
% multiple Input array elements. For the moment, we're just going to use
% one input type: random currents with a mean value of 330 pA, a standard
% deviation of 90 pA, and a time constant of 2 ms. The type of random
% current we use is an Ornstein Uhlenbeck process, so the input type is set
% to |'i_ou'| (we could also apply our input as a conductance |'g_ou'|, in
% which case we would also need to set an |E_reversal| parameter to set the
% reversal potential).

NeuronParams(1).Input(1).inputType = 'i_ou';
NeuronParams(1).Input(1).meanInput = 210;
NeuronParams(1).Input(1).stdInput = 80;
NeuronParams(1).Input(1).tau = 2;

%%
% Next we set the parameters for the 2nd neuron group, which represent
% basket interneurons. These cells' dendrites are not aligned along a
% particular axis, so we set the |axisAligned| parameters to be empty. We
% set the parameters to give this neuron type fast-spiking behaviour.

NeuronParams(2).modelProportion = 0.1;
NeuronParams(2).somaLayer = 1;
NeuronParams(2).axisAligned = '';
NeuronParams(2).neuronModel = 'adex';
NeuronParams(2).V_t = -50;
NeuronParams(2).delta_t = 2;
NeuronParams(2).a = 0.04;
NeuronParams(2).tau_w = 10;
NeuronParams(2).b = 40;
NeuronParams(2).v_reset = -65;
NeuronParams(2).v_cutoff = -45;
NeuronParams(2).numCompartments = 7;
NeuronParams(2).compartmentParentArr = [0 1 2 2 1 5 5];
NeuronParams(2).compartmentLengthArr = [10 56 151 151 56 151 151];
NeuronParams(2).compartmentDiameterArr = ...
  [24 1.93 1.95 1.95 1.93 1.95 1.95];
NeuronParams(2).compartmentXPositionMat = ...
[   0,    0;
    0,    0;
    0,  107;
    0, -107; 
    0,    0; 
    0, -107;
    0,  107];
NeuronParams(2).compartmentYPositionMat = ...
[   0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0];
NeuronParams(2).compartmentZPositionMat = ...
[ -10,    0;
    0,   56;
   56,  163;
   56,  163; 
  -10,  -66;
  -66, -173;
  -66, -173];
NeuronParams(2).C = 1.0*2.93;
NeuronParams(2).R_M = 15000/2.93;
NeuronParams(2).R_A = 150;
NeuronParams(2).E_leak = -70;
NeuronParams(2).dendritesID = [2 3 4 5 6 7];
NeuronParams(2).labelNames = {'dendritesID'};
NeuronParams(2).minCompartmentSize = 0.8;

NeuronParams(2).Input(1).inputType = 'i_ou';
NeuronParams(2).Input(1).meanInput = 200;
NeuronParams(2).Input(1).tau = 1;
NeuronParams(2).Input(1).stdInput = 20;

%% 
%Now setting up passive neurons which are dummy representatives of the inputs from other areas

NeuronParams(3) = NeuronParams(2); 
NeuronParams(3).somaLayer = 2;
NeuronParams(3).modelProportion = 0.3;
NeuronParams(3).V_t = -50;         % and different AdEx parameters
NeuronParams(3).delta_t = 2.2;
NeuronParams(3).a = 0.35;
NeuronParams(3).tau_w = 150;
NeuronParams(3).b = 40;
NeuronParams(3).v_reset = -70;
NeuronParams(3).v_cutoff = -45;

% they are going to be excitatory pyramidal dummy cells in this
% instance. Modify this to have them be whatever the incoming signals are
% generated by.
NeuronParams(3).neuronModel = 'passive'; %but with passive dynamics
NeuronParams(3).modelProportion = 0.1; 

NeuronParams(3).Input(1).inputType = 'i_ou'; %need to give it zero input.
NeuronParams(3).Input(1).meanInput = 190;
NeuronParams(3).Input(1).tau = 2;
NeuronParams(3).Input(1).stdInput = 60;

NeuronParams(4) = NeuronParams(2); % basket cells same in every layer
NeuronParams(4).somaLayer = 2;     % these are in layer 4
NeuronParams(4).modelProportion = 0.08;
NeuronParams(4).Input(1).inputType = 'i_ou';
NeuronParams(4).Input(1).meanInput = 200;
NeuronParams(4).Input(1).stdInput = 20;
NeuronParams(4).Input(1).tau = 1;

NeuronParams(5).somaLayer = 3; % Pyramidal cells in layer 5
NeuronParams(5).modelProportion = 0.1;
NeuronParams(5).axisAligned = 'z';
NeuronParams(5).neuronModel = 'adex';
NeuronParams(5).V_t = -52;
NeuronParams(5).delta_t = 2;
NeuronParams(5).a = 10;
NeuronParams(5).tau_w = 75;
NeuronParams(5).b = 345;
NeuronParams(5).v_reset = -60;
NeuronParams(5).v_cutoff = -47;
NeuronParams(5).numCompartments = 9;
NeuronParams(5).compartmentParentArr = [0 1 2 2 4 5 1 7 7];
NeuronParams(5).compartmentLengthArr = [35 65 152 398 402 252 52 186 186];
NeuronParams(5).compartmentDiameterArr = ...
  [25 4.36 2.65 4.10 2.25 2.4 5.94 3.45 3.45];
NeuronParams(5).compartmentXPositionMat = ...
[   0,    0;
    0,    0;
    0,  152;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0, -193;
    0,  193];
NeuronParams(5).compartmentYPositionMat = ...
[   0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0;
    0,    0];
NeuronParams(5).compartmentZPositionMat = ...
[ -35,    0;
    0,   65;
   65,   65;
   65,  463;
  463,  865;
  865, 1117;
  -35,  -87;
  -87, -193;
  -87, -193];
NeuronParams(5).C = 1.0*2.95;
NeuronParams(5).R_M = 20000/2.95;
NeuronParams(5).R_A = 150;
NeuronParams(5).E_leak = -70;
NeuronParams(5).dendritesID = [2 3 4 5 6 7 8 9];
NeuronParams(5).Input(1).inputType = 'i_ou';
NeuronParams(5).Input(1).meanInput = 650;
NeuronParams(5).Input(1).stdInput = 160;
NeuronParams(5).Input(1).tau = 2;


NeuronParams(6) = NeuronParams(2); % Basket cells in layer 5
NeuronParams(6).somaLayer = 3;
NeuronParams(6).modelProportion = 0.02;
NeuronParams(6).Input(1).inputType = 'i_ou';
NeuronParams(6).Input(1).meanInput = 200;
NeuronParams(6).Input(1).stdInput = 20;
NeuronParams(6).Input(1).tau = 1;
%% Connectivity
% Connectivity parameteres are specified as before, except that now we have
% several layers, the numbers can be specified per layer. Parameters that
% can specified on a per-layer basis are |axonArborRadius|,
% |axonArborLimit|, and |numConnectionsToAllFromOne|: 

ConnectionParams(1).axonArborRadius = [300, 200, 100];
ConnectionParams(1).axonArborLimit = [600, 400, 200];

%% Connectivity parameters
% We set the connectivity parameters in the same way as in tutorial 1, but
% this time we need to specify the parameters for connections between the
% two groups. First we set the parameters for connections from group 1 (the
% pyramidal cells) to itself:
ConnectionParams(1).numConnectionsToAllFromOne{1} = [1500,    0,    0];
ConnectionParams(1).numConnectionsToAllFromOne{2} = [ 250,    0,    0];
ConnectionParams(1).numConnectionsToAllFromOne{3} = [   0,   50,    0];
ConnectionParams(1).numConnectionsToAllFromOne{4} = [   0,   20,    0];
ConnectionParams(1).numConnectionsToAllFromOne{5} = [  25,    0,  175];
ConnectionParams(1).numConnectionsToAllFromOne{6} = [   0,    0,   25];

%%
% Most neurons only reside in their soma layer, so many of the values are
% zero. However, layer 5 neurons are large and so span all three layers.
% Layer 3 pyramidal neurons can make connections with layer 5 pyramidal
% neurons in layer 3 and layer 5. VERTEX automatically calculates which
% compartments of each neuron type are in each layer.
%
% The other connection parameters are specified as before:
% i_exp_stdp gives us synapses with spike timing depedent plasticity. We
% also now are required to specify the preRate (maximum change in the
% synaptic weight when the presynaptic neuron fires) and postRate (maximum
% change in synaptic weight when post synaptic neuron fires), and the time
% constants for the decay in weight change occuring when pre or post
% synaptic neuron fires (tPre and tPost). We also specify a wmin and wmax
% to apply upper and lower boundaries on the weight.
ConnectionParams(1).synapseType = ...
  {'i_exp_stdp', 'i_exp', 'i_exp_stdp', 'i_exp', 'i_exp_stdp', 'i_exp'};
ConnectionParams(1).targetCompartments = ...
  {NeuronParams(1).dendritesID, NeuronParams(2).dendritesID, ...
   NeuronParams(3).dendritesID, NeuronParams(4).dendritesID, ...
   NeuronParams(5).dendritesID, NeuronParams(6).dendritesID};
ConnectionParams(1).weights = {2, 30, 1, 15, 1, 15};
ConnectionParams(1).tau = {2, 1, 2, 1, 2, 1};
ConnectionParams(1).preRate = {-0.01,-0.01,-0.01,-0.01,-0.01,-0.01} ;% 0.01;
ConnectionParams(1).postRate = {0.01,0.01,0.01,0.01,0.01,0.01};% 0.01;
ConnectionParams(1).tPre = {2,2,2,2,2,2};
ConnectionParams(1).tPost = {6,6,6,6,6,6};
ConnectionParams(1).wmin = {0.01,0.01,0.01,0.01,0.01,0.01};
ConnectionParams(1).wmax = {100,100,100,100,100,100};
ConnectionParams(1).axonArborSpatialModel = 'gaussian';
ConnectionParams(1).sliceSynapses = true;
ConnectionParams(1).axonConductionSpeed = 0.3;
ConnectionParams(1).synapseReleaseDelay = 0.5;

%%
% And now we set the connectivity parameters for the other neuron groups:

ConnectionParams(2).axonArborRadius = [150, 0, 0];
ConnectionParams(2).axonArborLimit = [300, 0, 0];
ConnectionParams(2).numConnectionsToAllFromOne = ...
  {[2000, 0, 0], [200, 0, 0], [0, 20, 0], [0, 50, 0], [15, 0, 50], [0, 0, 10]};
ConnectionParams(2).synapseType = ...
  {'i_exp', 'i_exp', 'i_exp', 'i_exp', 'i_exp', 'i_exp'};
ConnectionParams(2).targetCompartments = ...
  {1, NeuronParams(2).dendritesID, ...
   1, NeuronParams(4).dendritesID, ...
   1, NeuronParams(6).dendritesID};
ConnectionParams(2).weights = {-4, -4, -4, -4, -4, -4};
ConnectionParams(2).tau = {6, 3, 6, 3, 6, 3};
ConnectionParams(2).axonArborSpatialModel = 'gaussian';
ConnectionParams(2).sliceSynapses = true;
ConnectionParams(2).axonConductionSpeed = 0.3;
ConnectionParams(2).synapseReleaseDelay = 0.5;

ConnectionParams(3).axonArborRadius = [200, 300, 200];
ConnectionParams(3).axonArborLimit = [400, 600, 400];
ConnectionParams(3).numConnectionsToAllFromOne = ...
  {[50, 0, 0], [5, 0, 0], [0, 500, 0], [0, 100, 0], [0, 10, 40], [0, 0, 5]};
ConnectionParams(3).synapseType = ...
  {'i_exp_stdp', 'i_exp', 'i_exp_stdp', 'i_exp', 'i_exp_stdp', 'i_exp'};
ConnectionParams(3).targetCompartments = ...
  {NeuronParams(1).dendritesID, NeuronParams(2).dendritesID, ...
   NeuronParams(3).dendritesID, NeuronParams(4).dendritesID, ...
   NeuronParams(5).dendritesID, NeuronParams(6).dendritesID};
ConnectionParams(3).weights = {1, 15, 2, 30, 2, 30};
ConnectionParams(3).tau = {2, 1, 2, 1, 2, 1};
ConnectionParams(3).preRate = {-0.01,-0.01,-0.01,-0.01,-0.01,-0.01} ;% 0.01;
ConnectionParams(3).postRate = {0.01,0.01,0.01,0.01,0.01,0.01};% 0.01;
ConnectionParams(3).tPre = {2,2,2,2,2,2};
ConnectionParams(3).tPost = {6,6,6,6,6,6};
ConnectionParams(3).wmin = {0.01,0.01,0.01,0.01,0.01,0.01};
ConnectionParams(3).wmax = {100,100,100,100,100,100};
ConnectionParams(3).axonArborSpatialModel = 'gaussian';
ConnectionParams(3).sliceSynapses = true;
ConnectionParams(3).axonConductionSpeed = 0.3;
ConnectionParams(3).synapseReleaseDelay = 0.5;

ConnectionParams(4).axonArborRadius = [0, 150, 0];
ConnectionParams(4).axonArborLimit = [0, 300, 0];
ConnectionParams(4).numConnectionsToAllFromOne = ...
  {[100, 0, 0], [10, 0, 0], [0, 450, 0], [0, 150, 0], [0, 10, 15], [0, 0, 5]};
ConnectionParams(4).synapseType = ...
  {'i_exp', 'i_exp', 'i_exp', 'i_exp', 'i_exp', 'i_exp'};
ConnectionParams(4).targetCompartments = ...
  {1, NeuronParams(2).dendritesID, ...
   1, NeuronParams(4).dendritesID, ...
   1, NeuronParams(6).dendritesID};
ConnectionParams(4).weights = {-4, -4, -4, -4, -4, -4};
ConnectionParams(4).tau = {6, 3, 6, 3, 6, 3};
ConnectionParams(4).axonArborSpatialModel = 'gaussian';
ConnectionParams(4).sliceSynapses = true;
ConnectionParams(4).axonConductionSpeed = 0.3;
ConnectionParams(4).synapseReleaseDelay = 0.5;

ConnectionParams(5).axonArborRadius = [100, 200, 300];
ConnectionParams(5).axonArborLimit = [200, 400, 600];
ConnectionParams(5).numConnectionsToAllFromOne = ...
  {[250, 0, 0], [30, 0, 0], [0, 50, 0], [0, 20, 0], [15, 0, 200], [0, 0, 100]};
ConnectionParams(5).synapseType = ...
  {'i_exp_stdp', 'i_exp', 'i_exp_stdp', 'i_exp', 'i_exp_stdp', 'i_exp'};
ConnectionParams(5).targetCompartments = ...
  {NeuronParams(1).dendritesID, NeuronParams(2).dendritesID, ...
   NeuronParams(3).dendritesID, NeuronParams(4).dendritesID, ...
   NeuronParams(5).dendritesID, NeuronParams(6).dendritesID};
ConnectionParams(5).weights = {1, 15, 1, 15, 2, 30};
ConnectionParams(5).tau = {2, 1, 2, 1, 2, 1};
ConnectionParams(5).preRate = {-0.01,-0.01,-0.01,-0.01,-0.01,-0.01} ;% 0.01;
ConnectionParams(5).postRate = {0.01,0.01,0.01,0.01,0.01,0.01};% 0.01;
ConnectionParams(5).tPre = {2,2,2,2,2,2};
ConnectionParams(5).tPost = {6,6,6,6,6,6};
ConnectionParams(5).wmin = {0.01,0.01,0.01,0.01,0.01,0.01};
ConnectionParams(5).wmax = {100,100,100,100,100,100};
ConnectionParams(5).axonArborSpatialModel = 'gaussian';
ConnectionParams(5).sliceSynapses = true;
ConnectionParams(5).axonConductionSpeed = 0.3;
ConnectionParams(5).synapseReleaseDelay = 0.5;

ConnectionParams(6).axonArborRadius = [0, 0, 150];
ConnectionParams(6).axonArborLimit = [0, 0, 300];
ConnectionParams(6).numConnectionsToAllFromOne = ...
  {[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 400], [0, 0, 40]};
ConnectionParams(6).synapseType = {[], [], [], [], 'i_exp', 'i_exp'};
ConnectionParams(6).targetCompartments = ...
  {[], [], [], [], NeuronParams(5).dendritesID, NeuronParams(6).dendritesID};
ConnectionParams(6).weights = {[], [], [], [], -3, -3};
ConnectionParams(6).tau = {[], [], [], [], 6, 3};
ConnectionParams(6).axonArborSpatialModel = 'gaussian';
ConnectionParams(6).sliceSynapses = true;
ConnectionParams(6).axonConductionSpeed = 0.3;
ConnectionParams(6).synapseReleaseDelay = 0.5;