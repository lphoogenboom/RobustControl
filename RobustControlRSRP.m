clear; clc; clear path;
load("Assignment_Data_SC42145.mat")

ops = sigmaoptions;
ops.Title.String = "Sigular values of W_i";
ops.Title.FontSize = 15;
ops.FreqUnits = 'Hz';
ops.InputLabels.FontSize = 12;
ops.OutputLabels.FontSize = 12;
ops.Grid = 'on';

%% 1.1

    % Block Diagram in Latex

%% 1.2

    % Equation in Latex

%% 1.3

% determine requires TFs
G = FWT(:,1:end-1); % MIMO TF
Gd = FWT(:,end); % TF for disturbance V

s = tf([1 0],1); % define s so we can easily setup TFs

Wp = [ (0.5556*s+2.0897)/(s+2.0897E-4) 0 ; 0 .2 ];
Wu = [ .01 0 ; 0  (5E-3*s^2+7E-4*s+5E-5)/(s^2+14E-4*s+10E-6)];

Wi = [((1/(16*pi))*s+0.3)/((1/(64*pi))*s+1) 0 ; 0 ((1/(16*pi))*s+0.3)/((1/(64*pi))*s+1)];
Wo = [(0.05*s+0.2)/(0.01*s+1) 0 ; 0 (0.05*s+0.2)/(0.01*s+1)];


% calc and plot singular values
[svi,freqi] = sigma(Wi); % singular values and frequencies for input uncertainty
[svo,freqo] = sigma(Wo); % singular values and frequencies for output uncertainty

% figure(1); clf; hold on; % plot singular values for uncertainties
% h1 = sigmaplot(Wi, ops, "b--o");
% h2 = sigmaplot(Wo, ops, "r--o");
% title("Singular values of W_i & W_o")
% legend("Wi","Wo");

H1 = [Wi(1,1)*ultidyn('D_1',[1,1]), 0; 0, Wi(2,2)*ultidyn('D_2',[1,1])];
H2 = [Wo(1,1)*ultidyn('D_3',[1,1]), 0; 0, Wo(2,2)*ultidyn('D_4',[1,1])];
delta = [H1, zeros(2); zeros(2), H2];
nom1 = H1(1,1).NominalValue;
figure(2); clf; hold on;
bode(H1, "b", .5*Wi, "r")
title("Frequency response of uncertainty W_i")
grid on;
% 
% figure(3); clf;
% bode(H2)
% title("Frequency response of uncertainty W_o")
% grid on;

