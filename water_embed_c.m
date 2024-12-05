function r_vec = water_embed_c(r_i, A, Delta, delta)

pho = A' * r_i;
pho_q = Delta * round((pho + delta) / Delta) - delta;
% pho_img = A * (pho_q - pho);

% pho = zeros(L,0);
% pho_q = pho;
% for i = 1 : L
%     pho(i) = dot(orginal, A(:,i));
%     pho_q(i) = (Delta * round((pho(i) - delta(i)) / Delta) + delta(i));
% end
% 
% 
[~,len_skm] = size(A);
[M,N] = size(r_i);
r_vec = zeros(M,N);
for i = 1 : len_skm
    r_vec = r_vec - pho(i)*A(:,i) + pho_q(i)*A(:,i);
end

%water_img = mod(water_img, 256);



end