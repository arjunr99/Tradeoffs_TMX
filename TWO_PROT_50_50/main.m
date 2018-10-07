%%% ------ 4/4 Transmultiplexer ----- %%%

clear;
close all;
clc;
N = 4;                      % Number of carriers
K = 4;                      % Upsampling

global h;
d1 = lt_poly({1;[0 1];[0 0 1]; [0 0 0 1]},0);

%d = lt_poly({1},-4);

%%%%%%%-----test-test -------%%%%%%%
v0 = rand(1,4) + 1j*rand(1,4);

fun = @solv_rot_in;
options = optimoptions('fminunc','Display','iter');
options.MaxFunEvals = 10000;
options.MaxIter = 4800;
options.TolFun = 1e-10;
options.TolX = 1e-10;
% v = fminunc(fun,v0,options);
% 
% h1 = unitary_mat(v(1:2),2);
% h2 = unitary_mat(v(3:4),2);

h1 = 1/sqrt(2)*[1 1; 1 -1];
h2 = 1/sqrt(2)*[1 1j; 1 -1j];

h = zeros(4,4);
h(1:2,1:2) = h1;
h(3:4,3:4) = h2;

w = dftmtx(2)/sqrt(2);
df = zeros(4,4);
df(1:2,1:2) = w;
df(3:4,3:4) = w;

v10 = rand(1,48); %+ 1j*rand(1,20);

fun = @solv_rot;
options = optimoptions('fminunc','Display','iter');
options.MaxFunEvals = 4000;
options.MaxIter = 4800;
options.TolFun = 1e-10;
options.TolX = 1e-10;
v1 = fminunc(fun,v10,options);

%v1 = v10;
lenv = length(v1);
U1 = cascade_rc(v1(1:lenv/2),4);
U2 = cascade_rc(v1(lenv/2+1:end),4);

U = diag(lt_poly({1,1,1,1},0));

U(1:2,1:2) = U1;
U(3:4,3:4) = U2;

P = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];

Y = d1'*P*U*h*P*df;

coef = GetCoefs(Y);

a = coef{1};
b = coef{2};
c = coef{3};
d = coef{4};
