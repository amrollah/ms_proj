

%Here is an example that may help. 
%A NARX network is trained on series inputs X and targets T, 
%then the simulation is picked up at the end of X using continuation input
%data X2 with a closed loop network. The final states after open loop simulation
%with X are used as the initial states for closed loop simulation with X2.

clear all
clc
load simplenarx_dataset


% DESIGN NETWORK
[x,t] = simplenarx_dataset;
net = narxnet;

[X,Xi,Ai,T] = preparets(net,x,{},t);
net = train(net,X,T,Xi,Ai);
view(net)

% SIMULATE NETWORK FOR ORIGINAL SERIES
[Y,Xf,Af] = sim(net,X,Xi,Ai);

plot(cell2mat(T)), hold on, plot(cell2mat(Y), 'r')


% CONTINUE SIMULATION FROM FINAL STATES XF & AF WITH ADDITIONAL
% INPUT DATA USING CLOSED LOOP NETWORK.
% Closed Loop Network
netc = closeloop(net);
view(netc)

% 10 More Steps for the first (now only) input
N=10;
X2 = num2cell(rand(1,N));

% Initial input states for closed loop continuation will be the
% first input's final states.
Xi2 = Xf(1,:);

% Initial 2nd layer states for closed loop contination will be the
% processed second input's final states.  Initial 1st layer states
% will be zeros, as they have no delays associated with them.
Ai2 = cell2mat(Xf(2,:));
for i=1:length(net.inputs{1}.processFcns)
  fcn = net.inputs{i}.processFcns{i};
  settings = net.inputs{i}.processSettings{i};
  Ai2 = feval(fcn,'apply',Ai2,settings);
end
Ai2 = mat2cell([zeros(10,2); Ai2],[10 1],ones(1,2));
% Closed loop simulation on X2 continues from open loop state after X.
Y2 = sim(netc,X2,Xi2,Ai2);


figure(2)
plot([cell2mat(t),nan(1,N);
       nan(1,length(t)),cell2mat(Y2)]');
legend('Original Targets','Network Predictions')
%      nan(1,length(t)),cell2mat(Y2);
%      nan(1,length(t)),cell2mat(Y2)]')
%legend('Original Targets','Network Predictions','Expected Outputs')
