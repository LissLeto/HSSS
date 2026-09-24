function [Gar,Gx,Gy]=GetGar(fig_ori)
fig_gray=rgb2gray(fig_ori);
Gx = imfilter(double(fig_gray), fspecial('sobel')',"replicate");
Gy = imfilter(double(fig_gray), fspecial('sobel'),"replicate");
Gar = sqrt(Gx.^2 + Gy.^2);
end
