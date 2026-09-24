function [labels,labelsj,wantedNum,para] = series_SP_fast(labels,img,wantedNum,para,level,labelsj,D_gap)

spNum = max(labels(:));
[M,N] = size(labels);

adj = double(GetAdj_mex(uint32(labels), spNum));
% adj = sparse(double(GetAdj_mex(uint32(labels), spNum)));

[L,A,B] = imsplit(img);
labels_vec = labels(:);

count_vec = accumarray(labels_vec, 1, [spNum,1], @sum, 0);

L_sum = accumarray(labels_vec, L(:), [spNum,1], @sum, 0);
A_sum = accumarray(labels_vec, A(:), [spNum,1], @sum, 0);
B_sum = accumarray(labels_vec, B(:), [spNum,1], @sum, 0);

L_mean_exact = L_sum ./ count_vec;
A_mean_exact = A_sum ./ count_vec;
B_mean_exact = B_sum ./ count_vec;

LAB = zeros(size(img), 'like', img);
LAB(:,:,1) = reshape(L_mean_exact(labels_vec), size(L));
LAB(:,:,2) = reshape(A_mean_exact(labels_vec), size(A));
LAB(:,:,3) = reshape(B_mean_exact(labels_vec), size(B));

%para = Getpara2(LAB);

L_mean_sp = L_sum ./ (count_vec + 0.0001);
A_mean_sp = A_sum ./ (count_vec + 0.0001);
B_mean_sp = B_sum ./ (count_vec + 0.0001);

sp_feature = [L_mean_sp, A_mean_sp, B_mean_sp, log(double(count_vec))];

psi_sp = Get_psi_sp_fast(adj, sp_feature, para);
[psi_sorted, order_idx] = sort(psi_sp, 'ascend');
psi_sp_order = [psi_sorted, order_idx];

sp_seed = find_sp_seed_fast(psi_sp, adj);

shortest_Distance = saveDshort_fast(adj, sp_feature, para);
median_D = median(shortest_Distance, 'omitnan');
if median_D == 0 || isnan(median_D)
    median_D = mean(shortest_Distance, 'omitnan');
end

%D_gap = 1;

seed_Cluster = cell(0,1);
ori_Cluster = num2cell((1:spNum).');
used = false(spNum,1);

for i = 1:spNum
    self = psi_sp_order(i,:);   % self(1)=psi, self(2)=label
    search = [];
    clusterII = [];
    shortestD = [];

    if used(self(2))
        continue
    end

    link_Vector = adj(self(2),:);

    while true
        search(end+1) = self(2);
        used(self(2)) = true;

        if size(seed_Cluster,1) + sum(~cellfun('isempty', ori_Cluster)) - length(search) + 1 == wantedNum(1)
            labels_turn = savelabels_fast(search,clusterII,seed_Cluster,ori_Cluster,M,N,labels,1);
            labelsj(:,:,end+1) = labels_turn;
            wantedNum(1) = [];
            if isempty(wantedNum)
                break
            end
        end

        if sp_seed(self(2)) == 1
            ori_Cluster(search) = {[]};
            seed_Cluster(end+1,1) = {search};
            break

        else
            linked_SP = find(link_Vector ~= 0);

            if numel(clusterII) <= 1
                if level < 5
                    Max_innerD = median_D;
                else
                    Max_innerD = (log(level)+1) * median_D;
                end
            else
                Max_innerD = max(shortestD);
            end

            clusterII = unique([clusterII, self(2)]);

            linked_SP_psi = psi_sp(linked_SP);
            higherN = linked_SP(linked_SP_psi > self(1));

            self_LAB = sp_feature(self(2),1:3);
            higherN_LAB = sp_feature(higherN,1:3);

            D2_sp_LAB = (higherN_LAB - self_LAB).^2 * para(:);
            [minD2, father_pos] = min(D2_sp_LAB);
            minD = sqrt(minD2);
            father = higherN(father_pos);

            if minD > D_gap * Max_innerD
                seed_Cluster(end+1,1) = {search};
                ori_Cluster(search) = {[]};
                break

            elseif used(father)
                for ii = 1:length(seed_Cluster)
                    seed_Cluster_ii = seed_Cluster{ii};
                    if any(seed_Cluster_ii == father)
                        seed_Cluster{ii} = [seed_Cluster{ii}, clusterII];
                        ori_Cluster(clusterII) = {[]};
                        break
                    end
                end
                break

            else
                clusterII = [clusterII, father];
                shortestD = [shortestD, minD];

                current_label = self(2);
                link_Vector(father) = 0;
                highlink_V = adj(father,:);
                highlink_V(current_label) = 0;
                link_Vector = link_Vector + highlink_V;

                self(2) = father;
                self(1) = psi_sp(father);
            end
        end
    end

    if isempty(wantedNum)
        break;
    end

    if size(seed_Cluster,1) + sum(~cellfun('isempty', ori_Cluster)) <= wantedNum(1)
        labels_turn = savelabels_fast(search,clusterII,seed_Cluster,ori_Cluster,M,N,labels,0);
        labelsj(:,:,end+1) = labels_turn;
        wantedNum(1) = [];
        if isempty(wantedNum)
            break;
        end
    end
end

labels = savelabels_fast(search,clusterII,seed_Cluster,ori_Cluster,M,N,labels,0);
end
function labels_turn = savelabels_fast(search,clusterII,seed_Cluster,ori_Cluster,M,N,labels,kind)

if kind == 1
    seed_Cluster(end+1,1) = {search};
    ori_Cluster(clusterII) = {[]};
end

new_Cluster = ori_Cluster(~cellfun('isempty', ori_Cluster));
seed_Cluster = vertcat(seed_Cluster, new_Cluster);

for ii = 1:numel(seed_Cluster)
    seed_Cluster{ii} = unique(seed_Cluster{ii});
end

maxLabel = max(labels(:));
map = zeros(maxLabel,1,'uint32');

for k = 1:numel(seed_Cluster)
    map(seed_Cluster{k}) = k;
end

labels_turn = reshape(map(labels(:)), M, N);
labels_turn = double(labels_turn);
end
%%
function psi_sp = Get_psi_sp_fast(adj, sp_feature, para)

spNum = size(sp_feature,1);
[row_idx, col_idx, w] = find(adj);

if isempty(row_idx)
    psi_sp = zeros(spNum,1);
    return;
end

diffLAB = sp_feature(row_idx,1:3) - sp_feature(col_idx,1:3);

D = sqrt(sum((diffLAB.^2) .* ((para(:)').^2), 2));

boundary_sum = full(sum(adj,2));
neighbor_num = full(sum(adj~=0,2));

weighted_D = D .* (w ./ boundary_sum(row_idx));
D_final = accumarray(row_idx, weighted_D, [spNum,1], @sum, 0);

psi_sp = zeros(spNum,1);
valid = (neighbor_num > 0) & (D_final > 0);
psi_sp(valid) = neighbor_num(valid) ./ D_final(valid);
end
%%
function sp_seed = find_sp_seed_fast(psi_sp, adj)

spNum = numel(psi_sp);
[row_idx, col_idx] = find(adj);

sp_seed = ones(spNum,1);

if isempty(row_idx)
    return;
end

has_higher_neighbor = psi_sp(col_idx) > psi_sp(row_idx);
sp_seed(unique(row_idx(has_higher_neighbor))) = 0;
end
%%
function shortestDistance = saveDshort_fast(linkMatrix, SP_featureMatrix, para)

SPNumber = size(SP_featureMatrix,1);
[row_idx, col_idx] = find(linkMatrix);

if isempty(row_idx)
    shortestDistance = nan(SPNumber,1);
    return;
end

diffLAB = SP_featureMatrix(row_idx,1:3) - SP_featureMatrix(col_idx,1:3);

alpha = para(1);
beta  = para(2);
gamma = para(3);

D = sqrt(alpha*(diffLAB(:,1).^2) + ...
         beta *(diffLAB(:,2).^2) + ...
         gamma*(diffLAB(:,3).^2));

shortestDistance = accumarray(row_idx, D, [SPNumber,1], @min, NaN);
end
%%
function sp_feature = Get_sp_feature_fast(LAB, spNum, labels)
L = LAB(:,:,1); A = LAB(:,:,2); B = LAB(:,:,3);
labels_vec = labels(:);

count_vec = accumarray(labels_vec, 1, [spNum,1], @sum, 0);
L_sum = accumarray(labels_vec, L(:), [spNum,1], @sum, 0);
A_sum = accumarray(labels_vec, A(:), [spNum,1], @sum, 0);
B_sum = accumarray(labels_vec, B(:), [spNum,1], @sum, 0);

L_mean = L_sum ./ (count_vec + 0.0001);
A_mean = A_sum ./ (count_vec + 0.0001);
B_mean = B_sum ./ (count_vec + 0.0001);

sp_feature = [L_mean, A_mean, B_mean, log(double(count_vec))];
end
