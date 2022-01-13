% clear workspace
clear; clc;

% load data
load('Assignment_Data_SC42145.mat');

% state space: A, B, C, D
% system: FWT
% wind data: Wind_Data (timeseries)

%% SISO design:

%% Question 1
[b,a] = ss2tf(A,B,C,D,1);

% formulate transfer function:
G1 = tf(b(1,:),a);
G = zpk(G1);

figure(1); clf;
title("Uncontrolled Open Loop system")
subplot(2,2,1);
bode(G1); grid on;
subplot(2,2,2);
pzmap(G1);  grid on;
subplot(2,2,[3,4])
step(G1); grid on;
title("Step Respone")

%% Question 2

%solve analytically
syms kp ki kd real
syms s real

G = -0.079878/(s+0.4104);

C = kp + ki/s + kd*s;

G1 = simplifyFraction(G*C);
G2 = simplifyFraction(G1/(1+G1));
%G2 = collect(G2);

clear syms

%% formulate transfer function:
s = zpk('s');

% formulate transfer function:
G1 = tf(b(1,:),a)
G = zpk(G1)

% PD:
kp = 5;
kd = 0;
ki = 2.65;
C = kp + kd*s + ki/s;

% zero cancelation - double notch
LP1 = (s^2 + 0.02113*s + 0.04101)/(s^2 - 0.007693*s + 0.04);
LP2 = (s^2 + 0.1664*s + 10.85)/(s^2 + 0.0492*s + 10.82);
LP3 = -1;

L = minreal(C*G*LP1*LP2*LP3,0.1);
T = minreal(L/(1+L),0.1);
S = minreal(1/(1+L),0.1);

[Gm,Pm,W180,Wc] = margin(L);
GMl = 1/abs(evalfr(L,W180));
[Ms, ~] = norm(S,Inf);
[Mt, ~] = norm(T,Inf);

P = pole(L);
Z = zero(L);

P = pole(T);
Z = zero(T);

figure(2); clf;
subplot(2,2,1)
bode(S); grid on;
title("Bode Diagram of Sensitivity Function")
% 
subplot(2,2,2)
bode(T); grid on;
title("Bode Diagram of Complementary Sensitivity Function")

% % check step response of tuned system:
subplot(2,2,[3,4])
step(T); grid on;
title("Step response of Complementary Sensitivity Function")

Evaluate = stepinfo(T);

figure(3); clf;
subplot(1,2,1)
bode(minreal(G*LP1*LP2,0.1)); grid on;
title("Bode Diagram of Filtered Plant")
subplot(1,2,2)
step(minreal(G*LP1*LP2,0.1)); grid on;
title("Step Response of Filtered Plant")

%% Q3

u = ones(1,1000);
t = 0:.1:100-.1;
figure(4); clf; hold on; grid on;
plot([-10, 0, 0 t], [0, 0, 1 u], 'r--');
y = lsim(G1,u,t);
plot(t,-y, "g")
y = lsim(T,u,t);
plot(t,y, "b")


title("Time-Domain Simulations of system")
xlabel("time [S]")
ylabel("\omega_r [rad/s]")
xlim([-10 100])
ylim([-0.1 1.02])
legend("Input","Uncontrolled","Notch+PI")

figure(5); clf; hold on; grid on;
plot(t,y, "b")
plot([4.1387,4.1387], [-.1,1.2], "b--")
plot([6.1784,6.1784], [-.1,1.2], "r--")
plot([0, 20], [1.0099, 1.0099], "g--")
title("Filtered and Controlled System")
xlim([0 20])
ylim([0.85, 1.02])
xticks([0 4.1387 5 6.1784 10 15 20])
legend("Controlled","Rise Time","Settling Time","Overshoot")

u = ones(1,100000);
t = 0:0.1:10000-0.1;
ssError = lsim(T,u,t);
ssError = 1-ssError(end);

%% Q4

G3 = tf(FWT(1,3));
L3 = minreal(C*G3*LP1*LP2*LP3,0.1);
T3 = minreal(L3/(1+L3),0.1);

u = ones(1,1000);
t = 0:.1:100-.1;
figure(6); clf; hold on; grid on;
plot([-10, 0, 0 t], [0, 0, 1 u], 'r--');
y = lsim(G3,u,t);
plot(t,y, "g")
y = lsim(T3,u,t);
plot(t,y, "b")

title("Time-Domain Step-Response Simulations of system")
xlabel("time [S]")
ylabel("\omega_r [rad/s]")
% xlim([-10 100])
% ylim([-0.1 1.02])
legend("Input","Uncontrolled","Notch+PI")


