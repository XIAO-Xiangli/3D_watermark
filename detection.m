function [w_aa,w_a] = detection(Points_w,Delta_num,key_skm,A,index_G,G,md,len_skm,S,T,index_sd,L,ConnectivityList,num_vertex_w)

[num,~] = size(Points_w);
Points_w_cen = sum(Points_w)/num;
Points_w_amend = Points_w - Points_w_cen;
Points_w_x = Points_w_amend(:,1);
Points_w_y = Points_w_amend(:,2);
Points_w_z = Points_w_amend(:,3);
C = [sum(Points_w_x.*Points_w_x),sum(Points_w_x.*Points_w_y),sum(Points_w_x.*Points_w_z);...
    sum(Points_w_x.*Points_w_y),sum(Points_w_y.*Points_w_y),sum(Points_w_y.*Points_w_z);...
    sum(Points_w_x.*Points_w_z),sum(Points_w_y.*Points_w_z),sum(Points_w_z.*Points_w_z)];
[e_vector, ~] = eig(C);
Points_w_amend_new = [Points_w_amend*e_vector(:,1),...
    Points_w_amend*e_vector(:,2),Points_w_amend*e_vector(:,3)];
[Points_w_amend_new(:,1),Points_w_amend_new(:,2),Points_w_amend_new(:,3)] = ...
    Dire_corr(Points_w_amend_new(:,1),Points_w_amend_new(:,2),Points_w_amend_new(:,3));
xyz_min = min(Points_w_amend_new);
xyz_max = max(Points_w_amend_new);
Points_w_amend_new = Points_w_amend_new - xyz_min + xyz_max;

% Spherical coordinates
Points_w_x_new = Points_w_amend_new(:,1);
Points_w_y_new = Points_w_amend_new(:,2);
Points_w_z_new = Points_w_amend_new(:,3);
% patch('vertices',Points_w_amend_new,'faces',ConnectivityList,'edgecolor','none',...
%     'facecolor',[1 1 1],'facelighting','phong')
% light
% axis equal off

r_w = sqrt(Points_w_x_new.*Points_w_x_new+Points_w_y_new.*Points_w_y_new+Points_w_z_new.*Points_w_z_new);
theta = acos(Points_w_z_new./r_w);

% Standardized sizing
r_w = r_w/(mean(r_w)/1000);

% Extract
length_vertex = length(theta);
repeat_t = floor(length_vertex/num_vertex_w);
rng(key_skm(repeat_t+2));
index = randperm(length_vertex);
% num_vertex_w = floor(length_vertex/repeat_t);
j = 1;
w_aa = zeros(L,repeat_t);
Delta = 10/Delta_num;
rng(key_skm(repeat_t+3));
delta = Delta * rand(len_skm,repeat_t) - Delta/2;
for i = 0:num_vertex_w:(repeat_t-1)*num_vertex_w
    index_i = index(i+1:i+num_vertex_w);
    r_i = r_w(index_i);
    SKm_w = SKm_w_Gen(index_G, md, len_skm, S, index_sd, L, key_skm(j));
    Bm = Bm_Gen(len_skm, T, S, SKm_w);
    G = sparse(G);
    GG = full(Bm * G);
    rw = A' * r_i;
    rw_q = Delta * round((rw + delta(:,j)) / Delta) - delta(:,j);
    diff = rw - rw_q;
    ee = GG' * diff;
    arb = sign(ee);
    w_aa(:,j) = arb;
    j = j + 1;
end
w_aa = w_aa';
w_a = sum(w_aa);
w_a = w_a';
for i = 1: L
    if w_a(i) < 0
        w_a(i) = -1;
    else
        w_a(i) = 1;
    end
end



end