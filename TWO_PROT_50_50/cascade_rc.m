function H = cascade_rc(theta,M)
st_num = length(theta);

H = rot_im(theta(1:4));

%DIST = lt_poly({[1 zeros(1,M-1) theta(st_num)] , 0 ; 0 [theta(st_num) zeros(1,M-1) 1] },0);

D = diag([lt_poly({1},0) lt_poly({1},-M)]);


if(length(theta)>4)

for l = 5:4:st_num
    
    H = rot_im(theta(l:l+3)) * D * H ;
    
end
end

end