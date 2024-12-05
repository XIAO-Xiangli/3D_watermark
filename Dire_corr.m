function [x_new,y_new,z_new] = Dire_corr(x,y,z)

x_mean = mean(x);
y_mean = mean(y);
z_mean = mean(z);
a_x = length(find(x>x_mean));
b_x = length(find(x<=x_mean));
if a_x > b_x
    x_new = -x;
else
    x_new = x;
end
a_y = length(find(y>y_mean));
b_y = length(find(y<=y_mean));
if a_y > b_y
    y_new = -y;
else
    y_new = y;
end
a_z = length(find(z>z_mean));
b_z = length(find(z<=z_mean));
if a_z > b_z
    z_new = -z;
else
    z_new = z;
end

end


