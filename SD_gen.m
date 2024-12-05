function [SD_w, A, index] = SD_gen(num_vertex_w, md, len_skm)

SD_w = rand(num_vertex_w, md+1);
SD_w = orth(SD_w);
SD_w(:,1) = [];

section = floor(md/len_skm);
A = zeros(num_vertex_w,len_skm);
index = zeros(1,len_skm);
j = 1;
for i = 1 : len_skm
    index_i = randperm(section,1);
    A(:,j) = SD_w(:,section*(i-1)+index_i);
    index(j) = index_i + section * (i-1);
    j = j + 1;
end
% index = sort(index);
% A = SD_w(:,index);

end

