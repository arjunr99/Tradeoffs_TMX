function H = rot_im(theta)

t = theta(1);
phi = theta(2);
phi1 = theta(3);
phi2 = theta(4);
rot = [cos(t) -sin(t)*exp(1j*phi); sin(t)*exp(-1j*phi) cos(t)];
rot2 = diag([exp(1j*phi1) exp(1j*phi2)]);

H = rot2 * rot;
end