function SKm = SKm_Gen(index_G, md, len_skm, S, T, L, key_skm)

rng(key_skm);
SKm = zeros(md*S,1);
section = floor(md/len_skm);
flag = 1;
h = 1;
for i = 1 : md
    if i <= section * len_skm
        source = find(index_G == h);
        ind = randperm(length(source));
        for j = 1:S
            SKm((i-1)*S+j) = source(ind(j));
        end
        if h <= L && mod(i,section) == 0 && flag == 1
            h = h + 1;
        end
        if h > L || flag == 0
            h = ceil(rand()*L);
            flag = 0;
        end
    else
        SKm((i-1)*S+1:i*S) = randperm(T, S);
    end
end

end