function [labels,psi,seeds,Gar]=fundationSP_accelerated( ...
    img,~,Gar,Gar_median,max_l,spatial_weight,search_radius_multiplier)
%FUNDATIONSP Generate the HSSS foundation segmentation.
%
% The last two arguments are optional. Their defaults reproduce the
% original implementation:
%   spatial_weight           = 0.5
%   search_radius_multiplier = 5

if nargin < 6 || isempty(spatial_weight)
    spatial_weight = 0.5;
end
if nargin < 7 || isempty(search_radius_multiplier)
    search_radius_multiplier = 5;
end
validateattributes(spatial_weight, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, ...
    mfilename, 'spatial_weight', 6);
validateattributes(search_radius_multiplier, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, ...
    mfilename, 'search_radius_multiplier', 7);

[L,A,B]=imsplit(img);
para = Getpara2(img,max_l);
%[Gar,~,~]=GetGar(fig_ori);

psi = psi_pixel_fast(img,para);

[M, N] = size(img,1:2);
labels = zeros(M,N);%labelCount = 0;
[seeds,range]=seed_elect_fast(psi,Gar,Gar_median);
seeds(1,:)=false;seeds(end,:)=false;seeds(:,1)=false;seeds(:,end)=false;

seedIdx = find(seeds);
labels(seedIdx) = 1:numel(seedIdx);

seeds_L=L(seeds);seeds_A=A(seeds);seeds_B=B(seeds);
[seed_rows,seed_cols]=find(seeds);
Seed=[seed_rows/M,seed_cols/N,seeds_L,seeds_A,seeds_B];
w_S=spatial_weight;%距离权重
% search_R=range*5/max(max([seed_rows,seed_cols]));
search_R = range * search_radius_multiplier / ...
    max(max([seed_rows; seed_cols]), 1);
weight=[w_S.^2,w_S.^2,para]';
nonSeedIdx = find(~seeds);
newLabels = assign_pixels_grid_mex( ...
    double(L), double(A), double(B), ...
    double(seed_rows), double(seed_cols), ...
    double(seeds_L), double(seeds_A), double(seeds_B), ...
    double(nonSeedIdx), double(search_R), double(weight), ...
    double(M), double(N));
labels(nonSeedIdx) = newLabels;
end
%%%%%%%%%%
function psi = psi_pixel_fast(LAB, para)
[m,n,~] = size(LAB);

L = LAB(:,:,1);
A = LAB(:,:,2);
B = LAB(:,:,3);

alpha = para(1);
beta  = para(2);
gamma = para(3);

Lacc = zeros(m,n);
Aacc = zeros(m,n);
Bacc = zeros(m,n);

for dr = -1:1
    for dc = -1:1
        if dr == 0 && dc == 0
            continue;
        end

        r0 = max(1,1-dr):min(m,m-dr);   % 当前中心像素位置
        c0 = max(1,1-dc):min(n,n-dc);
        r1 = r0 + dr;                   % 邻居位置
        c1 = c0 + dc;

        dL = L(r1,c1) - L(r0,c0);
        dA = A(r1,c1) - A(r0,c0);
        dB = B(r1,c1) - B(r0,c0);

        Lacc(r0,c0) = Lacc(r0,c0) + dL.^2;
        Aacc(r0,c0) = Aacc(r0,c0) + dA.^2;
        Bacc(r0,c0) = Bacc(r0,c0) + dB.^2;
    end
end

neighborCount = 8 * ones(m,n);
neighborCount(1,:)   = 5;
neighborCount(end,:) = 5;
neighborCount(:,1)   = 5;
neighborCount(:,end) = 5;
neighborCount(1,1)     = 3;
neighborCount(end,1)   = 3;
neighborCount(1,end)   = 3;
neighborCount(end,end) = 3;

Difference = sqrt(alpha*Lacc + beta*Aacc + gamma*Bacc);
psi = neighborCount ./ (Difference + 0.1);

max_psi_achieve = 10 * neighborCount;
XXX = find(psi == max_psi_achieve);
if ~isempty(XXX)
    XXX = XXX(randperm(length(XXX)));
    for kk = 1:numel(XXX)
        psi(XXX(kk)) = 80 + 0.0001 * kk;
    end
end
end
%%
function [seeds, range_big] = seed_elect_fast(psi, Gar,Gar_median)
[M, N] = size(psi);

%Gar_median = prctile(Gar(:), 95);
%Gar_median=mean(Gar(:));
range_big    = floor(max(M,N)/100)+1;
range_little = floor(min(M,N)/200);
range_little=max(1,range_little);

seeds_big = local_unique_max(psi, range_big);

if range_little == range_big
    seeds = seeds_big;
else
    seeds_little = local_unique_max(psi, range_little);
    seeds = seeds_big;
    mask_small = (Gar >= Gar_median);
    seeds(mask_small) = seeds_little(mask_small);
end
end

function s = local_unique_max(psi, r)
if r <= 0
    s = true(size(psi));
    return;
end

win = true(2*r+1);

% 局部最大值
localMax = ordfilt2(psi, numel(win), win);

% 统计窗口中有多少个位置等于局部最大值
countMax = conv2(double(psi == localMax), double(win), 'same');

% 唯一最大值才保留
s = (psi == localMax) & (countMax == 1);
end
