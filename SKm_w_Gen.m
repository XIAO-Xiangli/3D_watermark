function SKm_w = SKm_w_Gen(index_G, md, len_skm, S, index_sd, L, key_skm)

rng(key_skm);
SKm_w = zeros(len_skm*S,1);
section = floor(md/len_skm);
flag = 1;
h = 1;
k = 1;
for i = 1 : md
    if i <= section * len_skm
        source = find(index_G == h);
        ind = randperm(length(source));
        if index_sd(k) == i
            for j = 1:S
                % SKm((i-1)*S+j) = source(ind(j));
                SKm_w((k-1)*S+j) = source(ind(j));
            end
            if k < len_skm
                k = k+1;
            end
        end
        if h <= L && mod(i,section) == 0 && flag == 1
            h = h + 1;
        end
        if h > L || flag == 0
            h = ceil(rand()*L);
            flag = 0;
        end
%     else
%         SKm((i-1)*S+1:i*S) = randperm(T, S);
    end
end

end