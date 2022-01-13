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
inputvar = '[ud(2); kd(2); V(2); r(2); u(2)]';
input_to_G = '[ud + u]';
input_to_Gd = '[V(1)]';
input_to_Wp = '[G(1)+Gd(1)-r(1)+kd(1); G(2)+Gd(2)+kd(2)]';
input_to_Wu = '[u]';
input_to_Wi = '[u]';
input_to_Wo = '[G + Gd]'; % '[G + Gd]';
outputvar = '[Wi; Wo; Wp; Wu; -G(1)-Gd(1)+r(1)-kd(1); -G(2)-Gd(2)-kd(2)]'; % Output generalized plant
sysoutname = 'P';
sysic;

% mixed sensitivity controller (from A1)
systemnames = 'G Gd Wp Wu';
inputvar = '[V; r; u(2)]';
input_to_G = '[u]';
input_to_Gd = '[V]';
input_to_Wp = '[G(1)+Gd(1)-r; G(2)+Gd(2)]';
input_to_Wu = '[u]';
outputvar = '[Wp; Wu; -G(1)-Gd(1)+r; -G(2)-Gd(2)]'; % Output generalized plant
sysoutname = 'Pn';
sysic;

[Kn,CLn,GAMn,INFOn] = hinfsyn(Pn,2,2); % Hinf design


%% Question 5)
L = G*Kn;
S = inv(eye(size(L,1)) + L);
T = S*L;

N = lft(P,Kn);
p = pole(N);
z = tzero(N);

% P11 = P(1:8,1:8);
% P12 = P(1:8,9:end);
% P21 = P(9:end,1:8);
% P22 = P(9:end,9:end);
% 
% N = P11 + P12*Kn*inv(eye(2) - P22*Kn)*P21

M = N(1:4,1:4);
p = pole(M)

% Blockify incertainty
H = ultidyn('H',[1,1], "Bound", 1 );
H = [H, 0; 0, H];

delta = [Wi*H, zeros(2); zeros(2), Wo*H];

% generalized Nyquist
n = 100;
sample = 100;

f = logspace(-4,4,n); % frequency space
w = 2*pi*1i*f; % frequency omega j (rad/s)
S = zeros(4, n*sample); % L is 2x2 MIMO
Md = usample(M*delta, sample);

ec = 0;

for i=1:sample
    for k=1:n
        Lk = det(eye(4) - evalfr(Md(:,:,i,1),w(k)));
        if Lk == 0
            ec = ec + 1
        end

        S(:,k) = Lk;
    end
    % plot(real(S(2,:)), imag(S(2,:)), 'b')
    % plot(real(S(2,:)), -imag(S(2,:)), 'g')
end

figure(1)
plot(0,0,'x')
hold on
plot(real(S(2,:)), imag(S(2,:)), 'b')
plot(real(S(2,:)), -imag(S(2,:)), 'g')
hold off
axis equal
grid on
title('Generalized Nyquist plot')
xlabel('Re')
ylabel('Im')

% NS, N is internally stable
% NP, 

%% next 5)
N = lft(P,Kn);
N22 = N(4:8,4:8);
N22inf = norm(N22, Inf);

% not NP

n = 100;
omega = logspace(-4,4,n); % frequency space
Nf = frd(N, omega);

% for NP:
Nnp = Nf(4:8,4:8);
blk = [1 1; 1 1; 1 1; 1 1];
[mubnds,muinfo] = mussv(Nnp, blk, 'c');
[VDelta, VSigma, VLMI] = mussvextract(muinfo)

%% for NS
S1 = -inv(eye(size(G*Kn)) + G*Kn);
nRHP_S1 = numberRHP(S1)
S2 = Kn*inv(eye(size(G*Kn)) + G*Kn);
nRHP_S2 = numberRHP(S2)
S3 = -Kn*inv(eye(size(G*Kn)) + G*Kn);
nRHP_S3 = numberRHP(S3)
S4 = G*inv(eye(size(G*Kn)) + G*Kn);
nRHP_S4 = numberRHP(S4)
S5 = -Kn*inv(eye(size(G*Kn)) + G*Kn)*G;
nRHP_S5 = numberRHP(S5)
S6 = -G*Kn*inv(eye(size(G*Kn)) + G*Kn);
nRHP_S6 = numberRHP(S6)
S7 = -Kn*inv(eye(size(G*Kn)) + G*Kn)*Gd;
nRHP_S7 = numberRHP(S7)
S8 = -G*Kn*inv(eye(size(G*Kn)) + G*Kn)*Gd;
nRHP_S8 = numberRHP(S8)
S9 = inv(eye(size(G*Kn)) + G*Kn)*Gd;
nRHP_S9 = numberRHP(S9)
S10 = G*Kn*inv(eye(size(G*Kn)) + G*Kn);
nRHP_S10 = numberRHP(S10)

genNyquist(S10, 2)
%% for NP
N = lft(P,Kn);
n = 100;
omega = logspace(-4,4,n); % frequency space
Nf = frd(N, omega);

Nnp = Nf(5:8,5:8);

blk = [1 1; 1 1; 1 1; 1 1];

[mubnds,muinfo] = mussv(Nnp, blk, 'c');
muNP = mubnds(:,1);
[muNPinf,muNPw] = norm(muNP,inf)

%% for RS
N = lft(P,Kn);
n = 100;
omega = logspace(-4,4,n); % frequency space
Nf = frd(N, omega);

Nrs = Nf(1:4,1:4);

blk = [1 1; 1 1; 1 1; 1 1];

[mubnds,muinfo] = mussv(Nrs, blk, 'c');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf)

%% for Rp
N = lft(P,Kn);
n = 100;
omega = logspace(-4,4,n); % frequency space
Nf = frd(N, omega);

blk = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1];
% Nnp = Nf(4:8,4:8);

[mubnds,muinfo] = mussv(Nf, blk, 'c')
muRP = mubnds(:,1);
[muRPinf,muRPw] = norm(muRP,inf)

%% worst case sensitivity
delta = [ultidyn('del1',[1,1]), 0; 0, ultidyn('del1',[1,1])];
Np = lft(delta,N);
opt = wcgopt('Abadthreshold', 100);
Npw = wcgain(Np,opt)

%% functions
function nRHP = numberRHP(sys)
    nRHP = 0;
    p = pole(sys);
    n = size(p);

    for i = 1:n
        if p(i) > 0
            nRHP = nRHP + 1;
        end
    end
end

function genNyquist(sys, fignumber)
    s = size(sys, 1)
    L = tf(sys);

    n = 10000;
    f = logspace(-4,4,n); % frequency space
    w = 2*pi*1i*f; % frequency omega j (rad/s)
    M = zeros(s,n); % L is 2x2 MIMO
    
    for k=1:n
        Lk = det(eye(2)+evalfr(L,w(k)));
        M(:,k) = Lk;
    end
    
    % plot both parts of the semi circle
    
    figure(fignumber)
    plot(real(M(2,:)), imag(M(2,:)))
    hold on
    plot(real(M(2,:)), -imag(M(2,:)), 'g')
    plot(0,0,'x')
    hold off
    axis equal
    grid on
    title('Generalized Nyquist plot')
    xlabel('Re')
    ylabel('Im')
end





