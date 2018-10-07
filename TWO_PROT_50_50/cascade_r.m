function H = cascade_r(theta,M)
st_num = length(theta);

H = rot(theta(1));

%DIST = lt_poly({[1 zeros(1,M-1) theta(st_num)] , 0 ; 0 [theta(st_num) zeros(1,M-1) 1] },0);

D = diag([lt_poly({1},0) lt_poly({1},-M)]);


if(length(theta)>1)

for l = 2:st_num
    
    H = rot(theta(l)) * D * H ;
    
end
end

end