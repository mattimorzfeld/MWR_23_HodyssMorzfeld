clearvars
close all
clc

%% Nice Colors for plots
% .........................................................................
Color(1,:) = [213,29,38]/255;
Color(2,:) = [251,173,104]/255;
Color(3,:) = [49,124,180]/255;
Color(4,:) = [94,63,151]/255;
Color(5,:) = [17,139,59]/255;
Color(6,:) = [165,97,36]/255;
Color(7,:) = [27,158,119]/255;
Color(8,:) = [217,95,2]/255;
Color(9,:) = [117,112,179]/255;
% .........................................................................


%% Simulation parameters
%% ------------------------------------------
T = 6;
dt = 0.05;
Gap = 2;
t = 0:dt:T;
Steps = length(t);
%% ------------------------------------------

%% Model parameters
%% ------------------------------------------
F = 8;
n = 500;
%% ------------------------------------------

%% Observations
%% ------------------------------------------
var_y = 1;
skip = 1;
H = getH(skip,n);
R = var_y*eye(size(H,1));
%% ------------------------------------------

%% Truth and observations
%% ------------------------------------------
FileName = strcat('SpinUpObs_n',num2str(n),'.txt');
z=load(FileName);
FileName = strcat('SpinUpTruth_n',num2str(n),'.txt');
yAll = load(FileName);
%% ------------------------------------------


%% EnKF
%% ------------------------------------------
load(strcat('LongSim_n',num2str(n),'.mat'));

NeNoLoc = 1e5;
XoNoLoc = X(:,randi([1 length(X)],1,NeNoLoc));

[x,traceP,X,Xf] = myEnKF_nl(1,NeNoLoc,XoNoLoc,z,R,H,F,Gap,Steps,dt);
NoLocRMSE = sqrt(mean( (x(:,Gap+1:Gap:end)-yAll(:,Gap+1:Gap:end)).^2 ));
NoLocSpread = sqrt(traceP'/n);

figure, hold on
plot(NoLocRMSE(1:end),'Color',Color(1,:),'Linewidth',2)
plot(NoLocSpread(1:end),'--','Color',Color(1,:),'Linewidth',2)
set(gcf,'Color','w')
set(gca,'FontSize',20)
box off
xlabel('Cycle number')
ylabel('RMSE and spread')

FileName = strcat('SpinUpResults_n',num2str(n),'.mat');
save(FileName,'traceP','Xf','NoLocRMSE','NoLocSpread','dt','Gap','var_y')


