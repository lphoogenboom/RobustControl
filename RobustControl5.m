clear all

%% multivariable mixed-sensitivity
% load data
load('Assignment_Data_SC42145.mat');

% state space: A, B, C, D
% system: FWT
% wind data: Wind_Data (timeseries)

%% Generalized plant
G = FWT(:,1:end-1);
Gd = FWT(:,end);

% create performance weights:
s = tf('s');

% Wp:
Wp11 = (s/1.8 + 2.0897)/(s + 2.0897e-04);
Wp = [Wp11, 0;
      0   , 0.2];

% Wu:
Wu22 = ( 5*10^(-3)*s^2 + 7*10^(-4)*s + 5*10^(-5) ) / ...
       (s^2 + 14*10^(-4)*s +10^(-6));
Wu = [0.01, 0;
      0   , Wu22];

% Wi
Wid = ((1/(16*pi))*s + 0.3)/((1/(64*pi))*s + 1);
Wi = [Wid, 0;
      0   , Wid];

% Wo
Wod = (0.05*s + 0.2)/(0.01*s + 1);
Wo = [Wod, 0;
      0   , Wod];

% generalized plant
systemnames = 'G Gd Wp Wu Wi Wo';
inputvar = '[ud(2); kd(2); V; r; u(2)]';
input_to_G = '[ud + u]';
input_to_Gd = '[V]';
input_to_Wp = '[G(1)+Gd(1)+kd(1)-r; G(2)+Gd(2)+kd(2)]';
input_to_Wu = '[u]';
input_to_Wi = '[u]';
input_to_Wo = '[G + Gd]';
outputvar = '[Wi; Wo; Wp; Wu; -G(1)-Gd(1)+r-kd(1); -G(2)-Gd(2)-kd(2)]'; % Output generalized plant
sysoutname = 'P';
sysic;

H_1 = [ultidyn('D_1',[1,1]), 0; 0, ultidyn('D_2',[1,1])];
H_2 = [ultidyn('D_3',[1,1]), 0; 0, ultidyn('D_4',[1,1])];
delta = [H_1, zeros(2); zeros(2), H_2];

% iterate
Pcl = lft(delta,P);
[sscp,mu,logd0] = musyn(delta,2,2); 

%% Legacy:

% legacy: generalized plant
% systemnames = 'G Gd Wp Wu Wi Wo';
% inputvar = '[ud(2); kd(2); V(2); r(2); u(2)]';
% input_to_G = '[ud + u]';
% input_to_Gd = '[V(1)]';
% input_to_Wp = '[G(1)+Gd(1)-r(1); G(2)+Gd(2)]';
% input_to_Wu = '[u]';
% input_to_Wi = '[u]';
% input_to_Wo = '[G]'; % '[G + Gd]';
% outputvar = '[Wi; Wo; Wp; Wu; -G(1)-Gd(1)+r(1)-kd(1); -G(2)-Gd(2)-kd(2)]'; % Output generalized plant
% sysoutname = 'P';
% sysic;

% settings
n = 100;
omega = logspace(-4,4,n); % frequency space
blk = [1 1; 1 1; 1 1; 1 1; 2 2; 2 2];
d0 = 1;
iterations = 4;

% set placeholders
muRP = zeros(1,iterations);

% set initial condition
D = append(d0, d0, d0, d0, ...
        tf(eye(2)), tf(eye(2)), tf(eye(2)));

for i = 1:iterations
    % step 1
    [K, Nsc, gamma, info] = hinfsyn(D*P*inv(D),2,2, ...
                            'method','lmi','Tolgam',1e-3);
    Nf = frd(lft(P,K),omega);
    
    % step 2
    [mubnds,muinfo] = mussv(Nf, blk, 'c');
    % bodemag(mubnds(1,1), omega);
    muRP(i) = norm(mubnds(1,1),inf,1e-6)
    
    % step 3 
    [dsysl, dsysr] = mussvunwrap(muinfo);
    dsysl = dsysl/dsysl(3,3);
    dsys = fitfrd(genphase(dsysl(1,1)), 4);

    %D = append(tf(dsys), tf(dsys), tf(dsys), tf(dsys), ...
    %   tf(eye(2)), tf(eye(2)), tf(eye(2)))
end

%% Legacy: dksyn
% settings
n = 1000;
omega = logspace(-5,5,n); % frequency space

% Blockify incertainty
H_1 = [ultidyn('D_1',[1,1]), 0; 0, ultidyn('D_2',[1,1])];
H_2 = [ultidyn('D_3',[1,1]), 0; 0, ultidyn('D_4',[1,1])];
delta = [H_1, zeros(2); zeros(2), H_2];

% iterate
Punc = lft(delta,P);
opt = dkitopt('FrequencyVector',omega);
% [K, clp, bnd, dkinfo] = dksyn(Punc, 2, 2, opt);
