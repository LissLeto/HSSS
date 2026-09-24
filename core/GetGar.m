function [Gar,Gx,Gy]=GetGar(fig_ori)
fig_gray=rgb2gray(fig_ori);
Gx = imfilter(double(fig_gray), fspecial('sobel')',"replicate"); % 水平方向梯度
Gy = imfilter(double(fig_gray), fspecial('sobel'),"replicate");  % 垂直方向梯度
Gar = sqrt(Gx.^2 + Gy.^2); % 梯度幅值
end