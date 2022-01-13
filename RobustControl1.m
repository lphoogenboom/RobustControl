% clear workspace
clear all

% load data
load('Assignment_Data_SC42145.mat');

% state space: A, B, C, D
% system: FWT
% wind data: Wind_Data (timeseries)

%% SISO design:

%% Question 1
% figure(11)
% bode(FWT)

%% Question 2
[b,a] = ss2tf(A,B,C,D,1)
G1 = tf(b(1,:),a)

figure(1)
bode(G1)
grid on;

figure(2)
step(G1)
grid on;

figure(3)
pzmap(G1)
grid on;

%% solve analytically
syms kp ki kd real
syms s real

G = -0.079878/(s+0.4104);

C = kp + ki/s + kd*s;

G1 = simplifyFraction(G*C);
G2 = simplifyFraction(G1/(1+G1))
G2 = collect(G2);

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
% LP2 = 1; % (s^2 + 0.1664*s + 10.85)/(s^2 + 0.0492*s + 10.82);
% LP3 = -1;

% L = minreal(C*G*LP1*LP2*LP3,0.1)
L = minreal(-C*G*LP1, 0.1)
T = minreal(L/(1+L),0.1);
S = minreal(1/(1+L),0.1);

[Gm,Pm,W180,Wc] = margin(L)
GMl = 1/abs(evalfr(L,W180))
[Ms, ~] = norm(S,Inf)
[Mt, ~] = norm(T,Inf)

P = pole(L)
Z = zero(L)

P = pole(T);
Z = zero(T);

figure(4)
margin(S)
title('Bode plot of S')
grid on;

figure(5)
margin(T)
grid on;

figure(10)
pzmap(T)
grid on
title('Closed loop poles/zeros')
grid on;

Evaluate = stepinfo(T);
bandWidth = bandwidth(T);

% check step response of tuned system:

figure(6)
step(T)
grid on;

figure(7)
bode(minreal(G*LP1,0.1))
title('Bode plot of G1')
grid on;

figure(8)
step(minreal(G*LP1,0.1))
title('Step plot of G1')
grid on;

figure(9)
pzmap(minreal(L,0.1))
grid on
title('Open loop poles/zeros')
grid on;

% disturbance rejection
Gd = FWT(1,3)

figure(13)
step(S*Gd)
grid on;

%% sim comp

dt = 0.1;
t = (0:100-dt);
u = ones(1,length(t));
y = lsim(T,u,t);
yUncon = -lsim(G1,u,t);

filtered = minreal(G*LP1/(1+G*LP1));

figure(100); clf; hold on;
plot([-10 0 t],[0 0 u], 'r--', "lineWidth", 2)
plot(t,y, "b", 'LineWidth',2);

plot(t,yUncon, "g",  'LineWidth',2);
grid on;
xlim([-10 100])
ylim([-0.1 1.1])
legend("Unit Step Input","Uncontrolled", "Filtered" ,"Filtered+Control")
xlabel("Time [s]")
ylabel("\omega_r [rad/s]");
title("Unit Step Responses of Complementary Sensitivity Functions")




