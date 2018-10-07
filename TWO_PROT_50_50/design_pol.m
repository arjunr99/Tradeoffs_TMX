clear all;
close all;
clc;

N = 9;                      % Number of carriers
K = 8;                      % Upsampling

d = lt_poly({1},-72);
d1 = lt_poly({1;[0 1];[0 0 1]; [zeros(1,3) 1];[zeros(1,4) 1]; [zeros(1,5) 1];[zeros(1,6) 1]; [zeros(1,7) 1]; [zeros(1,8) 1]},0);
D = lt_poly({ 1; 1; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; [zeros(1,72) 1]},0);
D = diag(D);

d2 = lt_poly({1; [zeros(1,9) 1]; [zeros(1,18) 1]; [zeros(1,27) 1]; [zeros(1,36) 1] ;[zeros(1,45) 1] ; [zeros(1,54) 1] ; [zeros(1,63) 1]},0);
%d2 = diag(d2);
d3 = lt_poly({1 [zeros(1,63) 1] [zeros(1,54) 1] [zeros(1,45) 1] [zeros(1,36) 1] [zeros(1,27) 1] [zeros(1,18) 1] [zeros(1,9) 1] 1},0);
d3 = diag(d3);

W = dftmtx(8);
p = [1 zeros(1,7)]';
v0 = rand(1,8);
v0 = v0/sqrt((v0*v0'));

fun = @solv_rot_in;
options = optimoptions('fsolve','Display','iter');
options.MaxFunEvals = 57600;
options.MaxIter = 4800;
options.TolFun = 1e-10;
%options.TolX = 1e-10;
v = fsolve(fun,v0,options);

h1 = eye(8) - 2*(v')*v/(v*v');
h(1:K,1:K) = h1;

fun1 = @solv_rot;
options = optimoptions('fsolve','Display','iter');
options.MaxFunEvals = 57600;
options.MaxIter = 4800;
options.TolFun = 1e-10;
%options.TolX = 1e-10;
v0 = rand(1,9);
v0 = v0/sqrt((v0*v0'));
v = fsolve(fun1,v0,options);
H2 = eye(9) - 2*(v')*v/(v*v');

h2 = H2 * D * h ;

Y = d1' * d3 * h2 * d2;
coef = GetCoefs(Y);
a = coef{1};
%len = length(a)
x = fftshift(fft(a,1024));
plot(x);