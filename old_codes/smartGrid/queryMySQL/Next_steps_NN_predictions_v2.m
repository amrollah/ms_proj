

%Hi Jack,
%When using narxnet, the network performs only a one-step ahead prediction after it has been trained. Therefore, you need to use closeloop to perform a multi-step-ahead prediction and turn the network into parallel configuration.
%Take a look at this example for a multi-step-ahead prediction, N steps. This uses the dataset magdata.mat which is available in the Neural Network Toolbox. Also, some of the inputs will be used for performing the multi-step-ahead prediction, and results validated with the original data. I hope the comments help to understand.

%% 1. Importing data
S = load('magdata');
X = con2seq(S.u);
T = con2seq(S.y);
%% 2. Data preparation

N = 300; % Multi-step ahead prediction
% Input and target series are divided in two groups of data:
% 1st group: used to train the network
inputSeries  = X(1:end-N);
targetSeries = T(1:end-N);

% 2nd group: this is the new data used for simulation. inputSeriesVal will 
% be used for predicting new targets. targetSeriesVal will be used for
% network validation after prediction
inputSeriesVal  = X(end-N+1:end);
targetSeriesVal = T(end-N+1:end); % This is generally not available
%% 3. Network Architecture

delay = 2;
neuronsHiddenLayer = 50;
% Network Creation
net = narxnet(1:delay,1:delay,neuronsHiddenLayer);

%% 4. Training the network

[Xs,Xi,Ai,Ts] = preparets(net,inputSeries,{},targetSeries); 
net = train(net,Xs,Ts,Xi,Ai);
view(net)
Y = net(Xs,Xi,Ai); 

% Performance for the series-parallel implementation, only 
% one-step-ahead prediction
perf = perform(net,Ts,Y);

%% 5. Multi-step ahead prediction

inputSeriesPred  = [inputSeries(end-delay+1:end),inputSeriesVal];
targetSeriesPred = [targetSeries(end-delay+1:end), con2seq(nan(1,N))];
netc = closeloop(net);
view(netc)
[Xs,Xi,Ai,Ts] = preparets(netc,inputSeriesPred,{},targetSeriesPred);
yPred = netc(Xs,Xi,Ai);
perf = perform(net,yPred,targetSeriesVal);


figure;

plot([cell2mat(targetSeries),nan(1,N);
      nan(1,length(targetSeries)),cell2mat(yPred);
      nan(1,length(targetSeries)),cell2mat(targetSeriesVal)]')
  
legend('Original Targets','Network Predictions','Expected Outputs')