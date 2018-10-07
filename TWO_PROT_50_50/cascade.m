function H = cascade(v,M)
N=2;
st_num = length(v)/N;

H = eye(N);

for l=1:st_num
    
    v1 = v(((l-1)*N + 1) : l*N).';
    H = H * (eye(N) - ((v1*v1')/(v1'*v1))+ lt_poly({1},-M) * ((v1*v1')/(v1'*v1)));
end

end