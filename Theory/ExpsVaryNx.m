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


%% Setup
nos = 1e6;

Ne = 10;
nx = 11:10:101;
no = 10;
nr = 1;
ny = nr*no;
r = 5;

biasEnKF = zeros(length(nx),1);
biasEnKFnn = zeros(length(nx),1);
biasLocEnKF = zeros(length(nx),1);
biasLocEnKFnn = zeros(length(nx),1);
mseKF = zeros(length(nx),1);
mseEnKF = zeros(length(nx),1);
mseLocEnKF = zeros(length(nx),1);
biasPa11 = zeros(length(nx),1);
biasPa11Loc = zeros(length(nx),1);
biasPann = zeros(length(nx),1);
biasPannLoc = zeros(length(nx),1);

for kk = 1:length(nx)
    disp(kk)
    mseEnKF_tmp = zeros(nos,1);
    mseLocEnKF_tmp = zeros(nos,1);
    mseKF_tmp = zeros(nos,1);
    biasEnKF_tmp = zeros(nos,1);
    biasEnKFnn_tmp = zeros(nos,1);
    biasLocEnKF_tmp = zeros(nos,1);
    biasLocEnKFnn_tmp = zeros(nos,1);
    biasPa11_tmp = zeros(nos,1);
    biasPa11Loc_tmp = zeros(nos,1);
    biasPann_tmp = zeros(nos,1);
    biasPannLoc_tmp = zeros(nos,1);


    for jj=1:nos
        H = makeH(no,nr,nx(kk));
        Kt = H'/(H*H'+r*eye(ny));
        Pat = (eye(nx(kk))-Kt*H);

        %% Truth
        xt = randn(nx(kk),1);
        y = H*xt + sqrt(r)*randn(ny,1);

        %% KF
        xa = Kt*y;
        mseKF_tmp(jj) = mean((xa-xt).^2);

        %% EnKF
        X = randn(nx(kk),Ne);
        P = cov(X');
        K = (P*H')/(H*P*H'+r*eye(ny)); 
        xa = K*y;
        Pa = (eye(nx(kk))-K*H)*P; 
        mseEnKF_tmp(jj) = mean((xa-xt).^2);
        biasEnKF_tmp(jj) = K(1,1)-Kt(1,1);
        biasEnKFnn_tmp(jj) = K(end,end)-Kt(end,end);
        biasPa11_tmp(jj) = (xa(1)-xt(1))^2 - Pat(1,1);
        biasPann_tmp(jj) = (xa(nx(kk))-xt(nx(kk)))^2 - Pat(nx(kk),nx(kk));

        %% Loc. EnKF
        Ploc = eye(nx(kk)).*P;
        KLoc = (Ploc*H')/(H*Ploc*H'+r*eye(ny)); 
        xa = KLoc*y;
        PaLoc = (eye(nx(kk))-KLoc*H)*Ploc; 
        mseLocEnKF_tmp(jj) = mean((xa-xt).^2);
        biasLocEnKF_tmp(jj) = KLoc(1,1)-Kt(1,1);
        biasLocEnKFnn_tmp(jj) = KLoc(end,end)-Kt(end,end);
        biasPa11Loc_tmp(jj) = (xa(1)-xt(1))^2 - Pat(1,1);
        biasPannLoc_tmp(jj) = (xa(nx(kk))-xt(nx(kk)))^2 - Pat(nx(kk),nx(kk));

    end
    mseEnKF(kk) = mean(mseEnKF_tmp);
    mseLocEnKF(kk) = mean(mseLocEnKF_tmp);
    mseKF(kk) = mean(mseKF_tmp);

    biasEnKF(kk) = mean(biasEnKF_tmp);
    biasEnKFnn(kk) = mean(biasEnKFnn_tmp);
    biasLocEnKF(kk) = mean(biasLocEnKF_tmp);
    biasLocEnKFnn(kk) = mean(biasLocEnKFnn_tmp);
    biasPa11(kk) = mean(biasPa11_tmp);
    biasPa11Loc(kk) = mean(biasPa11Loc_tmp);
    biasPann(kk) = mean(biasPann_tmp);
    biasPannLoc(kk) = mean(biasPannLoc_tmp);
end


%% RMSE
thry = (1 - (nr./(r+nr)).*(no./nx))';

figure, hold on
plot(nx,mseEnKF-thry,'.-','Color',Color(9,:),'LineWidth',2,'MarkerSize',35)
plot(nx,mseLocEnKF-thry,'.-','Color',Color(7,:),'LineWidth',2,'MarkerSize',35)
% plot(nx,mseKF,'.-','Color',Color(8,:),'LineWidth',2,'MarkerSize',35)
set(gca,'XScale','log')
set(gca,'FontSize',18)
set(gcf,'Color','w')

thryNoLoc = (nr./(r+nr)).*((no+1)./(Ne-1)).*(1-(no./nx).*((nr.*(2*r+nr))./(r+nr).^2));
thryLoc = (1/(Ne-1)) .* (2*no./nx) .* ((r^2*nr)./(r+nr).^3); 
hold on, plot(nx,thryNoLoc,'k--','LineWidth',2)
hold on, plot(nx,thryLoc,'k-.','LineWidth',2)
% hold on, plot(nx,thry,'k-','LineWidth',2)

xlabel('Number of state variables')
ylabel('Excess post. MSE')
legend('EnKF (no loc.)','EnKF (loc.)',...
    'Theory (no loc.)', 'Theory (loc.)',...
    'Location','SouthEast')
ylim([-0.1 0.3])
box on
grid on
saveas(gcf, 'fig2_mse.tif')

%% Bias in the gain (observed)
figure, hold on
plot(nx,biasEnKF,'.-','Color',Color(9,:),'LineWidth',2,'MarkerSize',35)
plot(nx,biasLocEnKF,'.-','Color',Color(7,:),'LineWidth',2,'MarkerSize',35)

thry = -(r*nr./(r+nr).^3)*(no+1)/(Ne-1);
thryLoc = -(r*nr./(r+nr).^3)*2/(Ne-1);
plot(nx,thry*ones(size(nx)),'k--','LineWidth',2)
plot(nx,thryLoc*ones(size(nx)),'k-.','LineWidth',2)

set(gca,'XScale','log')
set(gca,'FontSize',18)
set(gcf,'Color','w')

xlabel('Number of state variables')
ylabel('Bias in the gain (observed)')
legend('EnKF (no loc.)','EnKF (loc.)','Theory (no loc.)', 'Theory (loc.)',...
    'Location','East')
% ylim([-.15 .02])
box on
grid on
saveas(gcf, 'fig2_biasGainObs.tif')

%% Bias in the gain (unobserved)
% figure, hold on
% plot(nx,biasEnKFnn,'.-','Color',Color(9,:),'LineWidth',2,'MarkerSize',35)
% plot(nx,biasLocEnKFnn,'.-','Color',Color(7,:),'LineWidth',2,'MarkerSize',35)
% 
% % thry = -(r*nr./(r+nr).^3)*(no+1)/(Ne-1);
% % thryLoc = -(r*nr./(r+nr).^3)*2/(Ne-1);
% % plot(ny,thry,'k--','LineWidth',2)
% % plot(ny,thryLoc*ones(size(ny)),'k.-','LineWidth',2)
% 
% set(gca,'XScale','log')
% set(gca,'FontSize',18)
% set(gcf,'Color','w')
% 
% xlabel('Number of observations')
% ylabel('Bias')
% legend('EnKF (no loc.)','EnKF (loc.)',...
%     'Location','SouthWest')

%% Bias in the posterior covariance (observed)
figure, hold on
plot(nx,biasPa11,'.-','Color',Color(9,:),'LineWidth',2,'MarkerSize',35)
plot(nx,biasPa11Loc,'.-','Color',Color(7,:),'LineWidth',2,'MarkerSize',35)

thry = (r^2*nr./(r+nr).^3)*(no+1)/(Ne-1);
thryLoc = (r^2*nr./(r+nr).^3)*2/(Ne-1);
plot(nx,thry*ones(size(nx)),'k--','LineWidth',2)
plot(nx,thryLoc*ones(size(nx)),'k-.','LineWidth',2)

set(gca,'XScale','log')
set(gca,'FontSize',18)
set(gcf,'Color','w')

xlabel('Number of state variables')
ylabel('Excess post. var. (observed)')
legend('EnKF (no loc.)','EnKF (loc.)','Theory (no loc.)', 'Theory (loc.)',...
    'Location','East')
ylim([-0.1 0.3])
box on
grid on
saveas(gcf, 'fig2_excessPostVarObs.tif')

%% Bias in the posterior covariance (unobserved)
figure, hold on
plot(nx,biasPann,'.-','Color',Color(9,:),'LineWidth',2,'MarkerSize',35)
plot(nx,biasPannLoc,'.-','Color',Color(7,:),'LineWidth',2,'MarkerSize',35)

thry = (nr./(r+nr)).*(no+1)/(Ne-1);
thryLoc = 0;
plot(nx,thry*ones(size(nx)),'k--','LineWidth',2)
plot(nx,thryLoc*ones(size(nx)),'k-.','LineWidth',2)

set(gca,'XScale','log')
set(gca,'FontSize',18)
set(gcf,'Color','w')

xlabel('Number of state variables')
ylabel('Excess post. var. (unobserved)')
legend('EnKF (no loc.)','EnKF (loc.)','Theory (no loc.)', 'Theory (loc.)',...
    'Location','East')
ylim([-0.1 0.3])
box on
grid on
saveas(gcf, 'fig2_excessPostVarUnobs.tif')
