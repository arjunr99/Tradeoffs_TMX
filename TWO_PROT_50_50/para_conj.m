function H = para_conj(M)

b = GetCoefs(M);
len = length(b{1,1});
[s0,s1] = size(b);

for l = 1:s0
    for m = 1:s1
       c{l,m} = fliplr((b{l,m}));
    end
end

H = lt_poly(c,len-1);
H = H';

end