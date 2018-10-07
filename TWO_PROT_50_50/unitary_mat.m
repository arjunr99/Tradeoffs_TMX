function [mat] = unitary_mat(v,N)

H = eye(N);
temp = eye(N);

mat = H;

for i = 1:N-1
        v1 = v(i*(i+1)/2 : (i+1)*(i+2)/2 - 1);
        temp(N-i:end,N-i:end) = eye(i+1) - 2*(v1')*v1/(v1*v1');
        mat = temp*mat;
end
        
end