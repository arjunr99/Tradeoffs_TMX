function F = solv_rot_in(v)

% N = 4;                      % Number of carriers
% K = 4;                      % Upsampling

d1 = lt_poly({1;[0 1];[0 0 1]; [0 0 0 1]},0);

h1 = unitary_mat(v(1:2),2);
h2 = unitary_mat(v(3:4),2);
h = zeros(4,4);
h(1:2,1:2) = h1;
h(3:4,3:4) = h2;

w = dftmtx(2)/sqrt(2);
df = zeros(4,4);
df(1:2,1:2) = w;
df(3:4,3:4) = w;

P = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];

Y = d1'*P*h*P*df;

coef = GetCoefs(Y);
a = coef{1};
c = coef{3};

af = fft(a,128);
cf = fft(c,128);

%x1 = sum(abs(af([32,64,96])).^2);
x2 = sum(abs(cf([1,64,96])).^2);

F = x2;
%F = [af([32,64,96]) cf([1,64,96])];

end