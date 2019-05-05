function [K, T] = PI(PM, Wc, G)
s = tf('s');
[~, phase_degree] = bode(G, Wc);
phase = phase_degree / 180 * pi;
T = tan(-pi/2 + PM - phase) / Wc;
F = 1 + 1/(T*s);
L = G * F;
[mag, ~] = bode(L, Wc);
K = 1 / mag;
end