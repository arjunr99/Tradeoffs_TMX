clear;
close all;
clc;
N = 4;                      % Number of carriers
K = 5;                      % Upsampling

d1 = lt_poly({1;[0 1];[0 0 1]; [zeros(1,3) 1];[zeros(1,4) 1]; [zeros(1,5) 1]; [zeros(1,6) 1]; [zeros(1,7) 1]},0);

v = rand(1,16);

fun = @solv_rot;
options = optimoptions('fminunc','Display','iter');
options.MaxFunEvals = 4000;
options.MaxIter = 4800;
options.TolFun = 1e-12;
options.TolX = 1e-12;
v1 = fminunc(fun,v,options);

%v1 = rand(1,16);

W = dftmtx(4);
len = length(v1);

h1 = cascade_r(v1(1:len/4),8)*[1/sqrt(2);1/sqrt(2)];
h2 = cascade_r(v1(len/4+1:len/2),8)*[1/sqrt(2);1/sqrt(2)];
h3 = cascade_r(v1(len/2 + 1: 3*len/4),8)*[1/sqrt(2);1/sqrt(2)];
h4 = cascade_r(v1(3*len/4 +1 :len),8)*[1/sqrt(2);1/sqrt(2)];

z = lt_poly({0},0);

H1 = [h1{1} z z z; z  h2{1} z z ; z z h3{1} z ; z z z h4{1} ; h1{2} z z z; z h2{2} z z ; z z h3{2} z; z z z h4{2}];


Y = d1' * H1 * W;
coef = GetCoefs(Y);
a = coef{1};
b = coef{2};
c = coef{3};
d = coef{4};

plot(20*log10(abs(fft([a;b;c;d].',1024))));
