function F = solv_rot(v)

epsilon = 32;

d1 = lt_poly({1;[0 1];[0 0 1]; [zeros(1,3) 1];[zeros(1,4) 1]},0);
%D = lt_poly({ 1; 1; 1 ; 1 ; [zeros(1,20) 1]},0);
%D = diag(D);

d2 = lt_poly({1 , [zeros(1,5) 1] ,  [zeros(1,10) 1] ,  [zeros(1,15) 1]},0);
d21 = diag(d2);
d3 = lt_poly({1 , [zeros(1,15) 1] , [zeros(1,10) 1] , [zeros(1,5) 1] , 1},0);
d31 = diag(d3);

W = dftmtx(4);

H1 = unitary_mat(v);

Y = d1'*d31*H1*d21*W;
coef = GetCoefs(Y);
a = coef{1};

fa = fft(a,1024);
x1 = sum(abs(fa(128+epsilon: 896-epsilon)).^2)/1024;

F = x1;

end