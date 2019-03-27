close all; clear; clc;

s = tf('s');
G = 3*(-s+1)/((5*s+1)*(10*s+1));
F = 1;

% S = feedback(eye(2), G*F);
% T = feedback(G*F, eye(2));

S = 1/(1+G*F);
T = G*F/(1+G/F);

figure
bode(G)

[m,p] = bode(G,3.1415/2);

[Gm,Pm,wp,wc] = margin(G*F);

figure
step(G)

figure
step(T)
