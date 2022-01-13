clear; clc;

%% multivariable mixed-sensitivity
% load data
load('Assignment_Data_SC42145.mat');

% state space: A, B, C, D
% system: FWT
% wind data: Wind_Data (timeseries)

%% part 1)
G = FWT(:,1:end-1);

% at w = 0:
G1 = evalfr(G, 0);
l1 = 1/(1 - (G1(1,2)*G1(2,1)/(G1(1,1)*G1(2,2)))  );
RGA1 = [l1 1-l1; 1-l1 l1];

% at w = 0.4*2*pi:
G2 = evalfr(G, 0.4*2*pi);
l1 = 1/(1 - (G2(1,2)*G2(2,1)/(G2(1,1)*G2(2,2)))  );
RGA2 = [l1 1-l1; 1-l1 l1];
clear l1;

%% part 2)
syms s real
% [y1; y2] = (1/d)*[M]*[u1; u2]
M = [-0.079878*(s^2 - 0.007693*s + 0.04)*(s^2 + 0.0492*s + 10.82), ...
     -0.0095639*(s^2 + 0.01085*s + 0.04)*(s^2 + 0.1651*s + 10.82);
     -0.048738*(s+0.6986)*(s^2 + 0.01562*s + 1.527), ...
     -0.0016139*(s^2 + 0.01562*s + 1.527)];
d = (s+0.4104)*(s^2 + 0.02113*s + 0.04101)*(s^2 + 0.1664*s + 10.85);

G = M/d;

sigma = svd(M) 

p = vpasolve(d==0,s)
z1 = vpasolve(sigma(1)==0,s)
z2 = vpasolve(sigma(2)==0,s)

% hence we hve 5 poles and 1 zero at -0.4105

%% part 3)
clear all
syms s wB real

wc = 0.4j*2*pi % cut-off frequency in rad/s
A = 1e-4; % disturbance attenuation
Hinf = 1.8; % bound on Hinf norm

Wp11 = (s/Hinf+wB)/(s+wB*A)

% find wB for the desired cut off frequency:
wBval = vpasolve(abs(subs(Wp11,s,wc))==1,wB) 
% wBval must be positive to ensure stable pole

Wp11 = subs(Wp11,wB,wBval)

% plot Wp11
clear syms

s = tf('s');
Wp11 = (0.5556*s - double(wBval))/(s - double(wBval)*10^(-4))
[Gm,Pm,W180,Wc] = margin(Wp11)



%% part 4 & 5:
clear G

% load data
load('Assignment_Data_SC42145.mat');
G = FWT; %(:,1:end-1);

% create performance weights:
s = tf('s');

% Wp:
Wp11 = (0.5556*s + 2.0897)/(s + 2.0897e-04);
Wp = [Wp11, 0;
      0   , 0.2];

% Wu:
Wu22 = (5*10^(-3)*s^2 + 7*10^(-4)*s + 5*10^-5) / ...
       (s^2 + 14*10^(-4)*s +10^(-5));
Wu = [0.01, 0;
      0   , Wu22];

Wp.InputName = {'r - w_r (rad/s)'; 'z (m)'};
Wp.OutputName = {'z11'; 'z12'};

Wu.InputName = {'B (deg)'; 't_e (Nm)'};
Wu.OutputName = {'z21'; 'z22'};


size(ss(Wp))
size(ss(Wu))

figure(1)
bode(1/Wp)

figure(2)
bode(Wu)

% systemnames = 'G Wp Wu'; % Define systems
% inputvar = '[w1; w2; u1; u2]'; % Input generalized plant
% input_to_G = '[u1; u2]';
% input_to_Wu = '[u1; u2]';
% input_to_Wp = '[w1 + G(1); w2 + G(2)]';
% outputvar = '[Wp; Wu; w1 + G(1); w2 + G(2)]'; % Output generalized plant
% sysoutname = 'P';
% sysic;

systemnames = 'G Wp Wu'; % Define systems
inputvar = '[r; V; B; T_e]'; % Input generalized plant
input_to_G = '[B; T_e; V]';
input_to_Wu = '[B; T_e]';
input_to_Wp = '[r - G(1); G(2)]';
outputvar = '[Wp; Wu; r - G(1); G(2)]'; % Output generalized plant
sysoutname = 'P';
sysic;

P.InputName = {'r (rad/s)'; 'V (m/s)'; 'B (deg)'; 'T_e (Nm'};
P.OutputName = {'z11'; 'z12'; 'z21'; 'z22'; 'err (rad/s)'; 'z (m)'};
P = minreal(P);
size(P)

[K2,CL2,GAM2,INFO2] = hinfsyn(P,2,2); % Hinf design

%% 2.6)
Pint = P;

figure(20); clf;
Pint(1,1) = P(1,1);
bode(Pint(1:2,1:2));

figure(21); clf;
bode(P(3:4,3:4))
grid on;

%% 2.7) 
K2 = minreal(K2);
CL2 = minreal(CL2);

% check stability
% P = eig(CL2)
% all eigenvalues have negavtive real part

% figure(4)
% nyquist(CL2(:,1))

% check in/ouputs and number of states
size(K2)
size(CL2)


% 2.8) using sysic

% for negative feedback:
% K2 = K2;

% % define signal names of K2
% K2.InputName = {'r-w_r'; 'z (m)'};
% K2.OutputName = {'B (deg)'; 't_e (Nm)'};
% 
% systemnames = 'FWT K2'; % Define systems
% inputvar = '[r; V]'; % Input generalized plant
% input_to_K2 = '[r - FWT(1); FWT(2)]';
% input_to_FWT = '[K2(1); K2(2); V]';
% outputvar = '[FWT(1); FWT(2); K2(1); K2(2)]'; % Output generalized plant
% sysoutname = 'CLsys';
% sysic;

% define signal names of K2


K2.InputName = {'r-w_r (rad/s)'; 'z (m)'};
K2.OutputName = {'B (deg)'; 't_e (Nm)'};

systemnames = 'P K2'; % Define systems
inputvar = '[r; V]'; % Input generalized plant
input_to_K2 = '[P(5); P(6)]';
input_to_P = '[r; V; K2(1); K2(2)]';
outputvar = '[P(5); P(6); K2(1); K2(2)]'; % Output generalized plant
sysoutname = 'CLsys';
sysic;

figure(7)
step(CLsys)

figure(8)
step(CL2)








