%% EX Basics
pm = 50;
wc = 0.4;
r = 0.00001;
ti = 10/wc; %%Tau_i
% EX1
s = tf('s');
G = 3*(-s+1)/((5*s+1)*(10*s+1));
F_lag = (ti*s+1)/(ti*s+r);
[mag, p_G] = bode(G, wc);
[mag, p_lag] = bode(F_lag, wc);
p_lead = pm -(p_G+p_lag-360+180);
p_lead = deg2rad(p_lead);
b = (1-sin(p_lead))/(1+sin(p_lead));
td = 1/(wc*sqrt(b));
K = 1/abs(evalfr(G*F_lag, j*wc))/abs(evalfr((td*s+1)/(b*td*s+1), j*wc))
F_lead = K*(td*s+1)/(b*td*s+1);
F = F_lead*F_lag;

figure(1)
margin(G)
hold on
margin(G*F)
legend('original','lag-lead')

% EX2
figure(2)
CL = G*F/(1+G*F);
step(CL)
bd = bandwidth(CL);
gpeak = getPeakGain(CL);
S = stepinfo(CL);
rt = S.RiseTime;
overshoot = S.Overshoot;
disp('=== system properties ===');
disp(['Bandwidth: ' num2str(bd)]);
disp(['Resonance Peak: ' num2str(gpeak)]);
disp(['Rise Time: ' num2str(rt)]);
disp(['Overshoot: ' num2str(overshoot) '%']);
%% EX Disturbance
%% EX1
s = tf('s');
G = 20/((s + 1)*((s/20)^2) + s/20 + 1);
Gd = 10/(s + 1);

[~, ~, ~, wc] = margin(Gd);
L = wc/s;
Fy = L/G;

% approximation to make controller proper
p1 = 100;
p2 = 100;
Fy_prop = p1*p2*Fy/((s + p1)*(s + p2));

% plot the result
figure(1)
margin(minreal(Fy*G));
hold on;
margin(minreal(Fy_prop*G));
legend('origin','adding poles')

figure(2)
bode(minreal(Gd/(1 + Fy*G)));
hold on;
bode(minreal(Gd/(1 + Fy_prop*G)));
legend('origin','adding poles')

figure(3)
step(minreal(Gd/(1 + Fy*G)));
hold on;
step(minreal(Gd/(1 + Fy_prop*G)));
legend('origin','adding poles')
%% EX2
c = 8;
s = tf('s');
G = 20/((s + 1)*((s/20)^2) + s/20 + 1);
Gd = 10/(s + 1);
Fy = (s+c)/s*Gd/G;
p1 = 100;
p2 = 100;
Fy_prop = p1*p2*Fy/((s + p1)*(s + p2));
figure(1)
bode(minreal(Gd/(1 + Fy*G)));
hold on;
bode(minreal(Gd/(1 + Fy_prop*G)));
legend('origin','adding poles')

figure(2)
step(minreal(Gd/(1 + Fy*G)));
hold on;
step(minreal(Gd/(1 + Fy_prop*G)));
legend('origin','adding poles')

% ex3
tau = 0.1;
p_lead = 15;
wc = 15;
p_lead = deg2rad(p_lead);
b = (1-sin(p_lead))/(1+sin(p_lead));
td = 1/(wc*sqrt(b));
K = 1/abs(evalfr(G*Fy_prop, 1i*wc))/abs(evalfr((td*s+1)/(b*td*s+1), j*wc));
Fy_lead = Fy_prop*K*(td*s+1)/(b*td*s+1);
Fr = 1/(1+tau*s);
S = 1/(1+G*Fy_lead);

figure(3)
margin(G*Fy_lead)
figure(4)
step(Fy_lead*Fr*G/(1+Fy_lead*G));
sys = Fy_lead*Fr*G/(1+Fy_lead*G);
Sys = stepinfo(sys);
rt = Sys.RiseTime;
overshoot = Sys.Overshoot;
disp(['Rise Time: ' num2str(rt)]);
disp(['Overshoot: ' num2str(overshoot) '%']);

figure(5)
step(Gd*S);

figure(6)
bode(S);
hold on
bode(G*Fy_lead*S);
legend('Sensitivity','Complementary Sensitivity');

figure(7)
subplot(2,1,1)
step(minreal(Fy_lead*Fr/(1 + Fy_lead*G)));
legend('r to y')
subplot(2,1,2)
step(minreal(Fy_lead*Gd/(1 + Fy_lead*G)));
legend('d to y');