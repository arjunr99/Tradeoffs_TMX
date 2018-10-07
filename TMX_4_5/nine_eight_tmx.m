clear;
close all;
clc;
N = 4;                      % Number of carriers
K = 5;                      % Upsampling

d1 = lt_poly({1;[0 1];[0 0 1]; [zeros(1,3) 1];[zeros(1,4) 1]},0);
%D = lt_poly({ 1; 1; 1 ; 1 ; [zeros(1,20) 1]},0);
%D = diag(D);

d2 = lt_poly({1 , [zeros(1,5) 1] ,  [zeros(1,10) 1] ,  [zeros(1,15) 1]},0);
d21 = diag(d2);
d3 = lt_poly({1 , [zeros(1,15) 1] , [zeros(1,10) 1] , [zeros(1,5) 1] , 1},0);
d31 = diag(d3);

v = rand(1,14);

fun = @solv_rot;
options = optimoptions('fminunc','Display','iter');
options.MaxFunEvals = 4000;
options.MaxIter = 4800;
options.TolFun = 1e-10;
options.TolX = 1e-10;
v1 = fminunc(fun,v,options);

W = dftmtx(4);

H1 = unitary_mat(v1);

Y = d1'*d31*H1*d21*W;
coef = GetCoefs(Y);
a = coef{1};
b = coef{2};
c = coef{3};
d = coef{4};

plot(20*log10(abs(fft([a;b;c;d].',1024)))-max(20*log10(abs(fft([a].',1024)))));
