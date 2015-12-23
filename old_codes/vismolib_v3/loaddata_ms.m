VMLCONF.measured_data = 'completeDataset2015.mat';

x1=load(conf.measured_data);
fn = fieldnames(x1);
x1 = x1.(fn{1});

data_Power = [x1.timestamp(:) x1.power(:)];
data_Irr = [x1.timestamp(:) x1.GHI(:)];
data_Temp = [x1.timestamp(:) x1.temperature(:)];
