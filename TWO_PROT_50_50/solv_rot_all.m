function F = solv_rot_all(v)

epsilon = 32;

% N = 4;                      % Number of carriers
% K = 4;                      % Upsampling

lenv = length(v);
d1 = lt_poly({1;[0 1];[0 0 1]; [0 0 0 1]},0);

h1 = unitary_mat(v(1:2),2);
h2 = unitary_mat(v(3:4),2);
h = zeros(4,4);
h(1:2,1:2) = h1;
h(3:4,3:4) = h2;

U1 = cascade(v(5:4+lenv/2),2,4);
U2 = cascade(v(5+lenv/2:end),2,4);

U = diag(lt_poly({1,1,1,1},0));

U(1:2,1:2) = U1;
U(3:4,3:4) = U2;
w = dftmtx(2)/sqrt(2);

df = zeros(4,4);
df(1:2,1:2) = w;
df(3:4,3:4) = w;

P = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];

Y = d1'*P*U*h*P*df;

coef = GetCoefs(Y);
a = coef{1};
c = coef{3};

af = fft(a,1024);
cf = fft(c,1024);

 x1 = sum(abs(af(128+epsilon: 896-epsilon)).^2)/1024;
 x2 = sum(abs(cf([1:128-epsilon,384+epsilon:end])).^2)/1024;

%F = [af(128+epsilon: 896-epsilon) cf([1:128-epsilon,384+epsilon:end])];

F = x1+x2;
end