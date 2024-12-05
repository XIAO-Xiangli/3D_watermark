clc;clear;
L = 128;
Delta_num = 1000;
sigma_E = 0.5;
T = 3000; S = 4;
md = 256;
len_skm = 128;
num_vertex_w = 640;

% load the file
F = stlread('xiaomei.stl');
Points = F.Points;
ConnectivityList = F.ConnectivityList;
% patch('vertices',Points,'faces',ConnectivityList,'edgecolor','none',...
%     'facecolor',[0.99 0.99 0.99],'facelighting','phong')
% light
% axis equal off

% PCA
tic;
[num,coor] = size(Points);
Points_cen = sum(Points)/num;
Points_amend = Points - Points_cen;
Points_x = Points_amend(:,1);
Points_y = Points_amend(:,2);
Points_z = Points_amend(:,3);
C = [sum(Points_x.*Points_x),sum(Points_x.*Points_y),sum(Points_x.*Points_z);...
    sum(Points_x.*Points_y),sum(Points_y.*Points_y),sum(Points_y.*Points_z);...
    sum(Points_x.*Points_z),sum(Points_y.*Points_z),sum(Points_z.*Points_z)];
[e_vector, ~] = eig(C);
Points_amend_new = [Points_amend*e_vector(:,1),...
    Points_amend*e_vector(:,2),Points_amend*e_vector(:,3)];
[Points_amend_new(:,1),Points_amend_new(:,2),Points_amend_new(:,3)] = ...
    Dire_corr(Points_amend_new(:,1),Points_amend_new(:,2),Points_amend_new(:,3));
xyz_min = min(Points_amend_new);
xyz_max = max(Points_amend_new);
Points_amend_new = Points_amend_new - xyz_min + xyz_max;

% %display
% R_xyz = rotx(-60)*roty(100)*rotz(-18);
% Points_amend_new = Points_amend_new * R_xyz;
% patch('vertices',Points_amend_new,'faces',ConnectivityList,'edgecolor','none',...
%     'facecolor',[0.7 0.7 0.7],'facelighting','phong')
% light
% axis equal off
% F_r = triangulation(ConnectivityList,Points_amend_new);
% stlwrite(F_r,'shenli_o.stl');

% Spherical coordinates
Points_x_new = Points_amend_new(:,1);
Points_y_new = Points_amend_new(:,2);
Points_z_new = Points_amend_new(:,3);
r = sqrt(Points_x_new.*Points_x_new+Points_y_new.*Points_y_new+Points_z_new.*Points_z_new);
theta = acos(Points_z_new./r);
pha = atan(Points_y_new./Points_x_new);

% Standardized sizing
r = r/(mean(r)/1000);
toc;

tic;
%generate ELUT, SKm, G, SD, A, and delta
G_zero = 1;
index_G = zeros(T,1);
while isempty(G_zero) == 0
    G = zeros(T,L);
    for i=1:T
        j=randi([1,L]);
        G(i,j)=1;
        index_G(i) = j;
    end
    sum_G =sum(G);
    G_zero = find(sum_G<S);
end

Delta = 10/Delta_num;
length_vertex = length(theta);
%num_vertex_w = floor(length_vertex/repeat_t);
repeat_t = floor(length_vertex/num_vertex_w);
[SD_w, A, index_sd] = SD_gen(num_vertex_w, md, len_skm);
key_skm = randi([1,1000000],repeat_t+3,1);
rng(key_skm(repeat_t+3));
delta = Delta * rand(len_skm,repeat_t) - Delta/2;
toc;

tic;
rng( "shuffle");
w = sign(randi([0,1],L,1)-1/2);
WLUT =  Delta/(4*S) * G * w;


% Encryption
r_e = r;
rng(key_skm(repeat_t+2));
index_th = randperm(length_vertex);
for k = 0:num_vertex_w:(repeat_t-1)*num_vertex_w
    index_i = index_th(k+1:k+num_vertex_w);
    r_i = r(index_i);
    r_vec = water_embed_c(r_i, A, Delta, delta(:,k/num_vertex_w+1));
    SKm = SKm_Gen(index_G, md, len_skm, S, T, L, key_skm(k/num_vertex_w+1));
    pad_e = zeros(md,1);
    for i = 1:md
        for j = 1:S
            pad_e(i) = pad_e(i) + WLUT(SKm((i-1)*S+j));
        end
        r_i = r_i + pad_e(i) * SD_w(:,i);
    end
    r_e(index_i) = r_i + r_vec;
end
index_i = index_th(k+num_vertex_w+1:length_vertex);
r_i = r(index_i);
SKm = SKm_Gen(index_G, md, len_skm, S, T, L, key_skm(k/num_vertex_w+2));
pad_e = zeros(md,1);
for i = 1:md
    for j = 1:S
        pad_e(i) = pad_e(i) + WLUT(SKm((i-1)*S+j));
    end
    r_i = r_i + pad_e(i) * SD_w(1:length(r_i),i);
end
r_e(index_i) = r_i;
toc;

% Display
x_e = r_e .* sin(theta) .* cos(pha);
y_e = r_e .* sin(theta) .* sin(pha);
z_e = r_e .* cos(theta);
Points_e = [x_e,y_e,z_e];
% %display
% R_xyz = rotx(-60)*roty(100)*rotz(-18);
% Points_e = Points_e * R_xyz;
% patch('vertices',Points_e,'faces',ConnectivityList,'edgecolor','none',...
%     'facecolor',[0.7 0.7 0.7],'facelighting','phong')
% light
% axis equal off

% Attack
R_xyz = rotx(60)*roty(20)*roty(30);
Points_d = Points_e * R_xyz;
Points_d = 1.256 * Points_d + [100,-100,0];
% patch('vertices',Points_dd,'faces',ConnectivityList,'edgecolor','none',...
%     'facecolor',[1 1 0],'facelighting','phong')
% light
% axis equal off

% Extract the watermark
tic;
[w_aa,w_a] = detection(Points_d,Delta_num,key_skm,A,index_G,...
    G,md,len_skm,S,T,index_sd,L,ConnectivityList,num_vertex_w);
length(find(w~=w_a));
toc;











