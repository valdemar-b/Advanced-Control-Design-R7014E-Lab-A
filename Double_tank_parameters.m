clc
close all

fSamplingFrequency = 80;
fSamplingTime = 1/fSamplingFrequency;
A = 33; % 33 square centimetres
a = 0.16; % 0.16 square centimetres
T = 35; % 35 seconds time constant
K = 5; % metres per volt, motor gain
h = 15;
h0 = 3.2;
g = 982;

alpha = g*a/(A*sqrt(2*g*(h+h0)));
% given pump values
% fu0 = [16.47, 20.75, 24.52, 27.63, 30.06, 32.49, 34.80, 40.28, 48.24, 52.43];
% u0 = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 5, 8, 10];

% measured pump values:
fu0 = [15.7913, 21.36471, 25.94286, 30.26667, 33.01818, 36.32, 38.6383];
u0 = [0.5, 1, 1.5, 2, 2.5, 3, 3.5];

figure;
plot(u0, fu0, 'k-', 'LineWidth', 2)
p1 = polyfit(u0, fu0, 3);
px1 = min(u0) :.1: max(u0);
py1 = polyval(p1, px1);
hold on
grid on
plot(px1, py1, 'r-', 'LineWidth', 2)
legend('Flow', 'Polyfitted flow')
xlabel('Voltage')
ylabel('Flow (cm^3/s')

u_h = 0;
p1_alt = p1;
p1_alt(4) = p1(4)-a*sqrt(2*g*(h+h0));
r = roots(p1_alt);
%H = fu/A * u_h
for n = 1:length(r)
    if abs(r(n)) == r(n)
        u_h = r(n);
        disp(u_h)
    end
end

fu = p1(1)*u_h^3 + p1(2)*u_h^2 + p1(3)*u_h + p1(4);

fuprime = p1(1)*3*u_h^2 + p1(2)*2*u_h + p1(3);
beta = fuprime/A;
zeta = sqrt(2*g*(h+h0))/A;

Alin = [-alpha, 0; alpha, -alpha];
Blin = [beta; 0];
Clin = [0, 1];
Dlin = 0;

Aa = [-alpha, 0, zeta; alpha, -alpha, 0; 0, 0, 0];
Bb = [beta; 0; 0];
Nn = [0; 0; 1];
Cc = [0, 1, 0];

gH = [0.56, 2.18, 2.9, 3.7, 4.5, 6.1, 6.8, 7.6, 9.2, 9.9];
H = [3, 5, 6, 7, 8, 10, 11, 12, 14, 15];

figure;
plot(gH, H, 'k-', 'LineWidth', 2)
p2 = polyfit(gH, H, 1);
px2 = min(gH): .1: max(gH);
py2 = polyval(p2, px2);
hold on
grid on
plot(px2, py2, 'r-', 'LineWidth', 2)
xlabel('Voltage')
ylabel('Height')



SS = ss(Alin, Blin, Clin, Dlin);
gs = tf(SS);
[num, den] = tfdata(gs, 'v');
SSd = c2d(SS, fSamplingTime);
[Ad, Bd, Cd, Dd] = ssdata(SSd);
R1 = eye(2)*0.01;
R2 = 0.01;
SSs = ss(Aa, Bb, Cc, Dd);
SSsd = c2d(SSs, fSamplingTime);
[Aad, Bbd, Ccd, Ddd] = ssdata(SSsd);


Qd = R1; % det var typ R1-matrisen i think
Rd = R2; % varians för noise
Sd = 0; % spelar inte så mycket roll
Ed = eye(2,2); % ej heller så viktig
N = Ed;



[P, Kgain, L] = idare(Ad', Cd', Qd, Rd, [], []);
[PC, KgainC, LC] = icare(Alin', Clin', Qd, Rd, [], []);
%Kgain = Kgain'
Kgain = Cd'*(Cd*P* Cd' + R2)

Qqd = eye(3,3)*0.01;
Rrd = R2;
NN = eye(3,3);
[Pp, Kkgain, Ll] = idare(Aad', Ccd', Qqd, Rrd, [], []);
Kkgain = Ccd'*(Ccd*Pp* Ccd' + R2)


%% Estimation error for non-linear model
% skapa en loop och optimera denna del pls så den bara kan iterera över all
% data automatiskt
samples = 8000; %to make vectors of equal length, to remove 1-5 elements in the end
time = out.h1hat(1:samples, 1);
h1 = out.h1(1:samples, 2);
h2 = out.h2(1:samples, 2);

% Stationary Kalman filter
h1hat = out.h1hat(1:samples, 2);
h2hat = out.h2hat(1:samples, 2);

% Non-stationary Kalman filter
h1hatn = out.h1hatn(1:samples, 2);
h2hatn = out.h2hatn(1:samples, 2);

% RMS error for stationary
RMSerror_h1 = rmse(h1, h1hat);
RMSerror_h2 = rmse(h2, h2hat);

% RMS error for non-stationary
RMSerror_h1n = rmse(h1, h1hatn);
RMSerror_h2n = rmse(h2, h2hatn);

RMS_h1 = [RMSerror_h1, RMSerror_h1n]
RMS_h2 = [RMSerror_h2, RMSerror_h2n]

% determining settling time for 2% from data

% settling time for simulated data
for i = length(h1):-1:1
    if h1(i) <= h*0.98 || h1(i) >= h*1.02
        Ts1 = time(i);
        h1_s = h1(i);
        break
    end
end
for i = length(h2):-1:1
    if h2(i) <= h*0.98 || h2(i) >= h*1.02
        Ts2 = time(i);
        h2_s = h2(i);
        break
    end
end

% settling time for estimation, stationary KF
for i = length(h1hat):-1:1
    if h1hat(i) <= h*0.98 || h1hat(i) >= h*1.02
        Ts1hat = time(i);
        h1hat_s = h1hat(i);
        break
    end
end
for i = length(h2hat):-1:1
    if h2hat(i) <= h*0.98 || h2hat(i) >= h*1.02
        Ts2hat = time(i);
        h2hat_s = h2hat(i);
        break
    end
end

% settling time for estimation, non-stationary KF

for i = length(h1hatn):-1:1
    if h1hatn(i) <= h*0.98 || h1hatn(i) >= h*1.02
        Ts1hatn = time(i);
        h1hatn_s = h1hatn(i);
        break
    end
end
for i = length(h2hatn):-1:1
    if h2hatn(i) <= h*0.98 || h2hatn(i) >= h*1.02
        Ts2hatn = time(i);
        h2hatn_s = h2hatn(i);
        break
    end
end



settling_time_h1 = [Ts1, Ts1hat, Ts1hatn]
settling_time_h2 = [Ts2, Ts2hat, Ts2hatn]


bias1 = zeros(1,length(h1));
bias2 = zeros(1,length(h2));
bias1n = zeros(1,length(h1hatn));
bias2n = zeros(1,length(h2hatn));
for i = length(h1)
    bias1(i) = h1(i)-h1hat(i);
    bias2(i) = h2(i)-h2hat(i);
    bias1n(i) = h1(i)-h1hatn(i);
    bias2n(i) =  h2(i)-h2hatn(i);
end

bias1avg = sum(bias1)/length(bias1)
bias2avg = sum(bias2)/length(bias2)
bias1navg = sum(bias1n)/length(bias1n)
bias2navg = sum(bias2n)/length(bias2n)
bias1std = std(bias1)
bias2std = std(bias2)
bias1nstd = std(bias1n)
bias2nstd = std(bias2n)
