function [labels, labelsj, wantedNum] = series_SP_fast_size( ...
    labels, img, wantedNum, para, level, labelsj, D_gap, ...
    lambda_size, size_tau)
%SERIES_SP_FAST_SIZE Experimental HSSS merge with a soft size constraint.
%
% This file is intentionally isolated from the original series_SP_fast.m.
% The original Lab-based hierarchy is retained. Candidate merges whose
% estimated area exceeds size_tau times the target average area receive
% an additive penalty scaled by the median neighboring Lab distance.

% lambda_size = 0 reproduces the original candidate score.

% The size penalty is
%   lambda_size * median_D * max(0, merged_area/target_area-size_tau)^2.

% It affects candidate selection and the merge stopping test, but it does
% not alter the Lab statistics or the representative-density definition.


spNum = max(labels(:));
[M, N] = size(labels);
adj = double(GetAdj_mex(uint32(labels), spNum));

[L, A, B] = imsplit(img);
labels_vec = labels(:);
count_vec = accumarray(labels_vec, 1, [spNum, 1], @sum, 0);

L_sum = accumarray(labels_vec, L(:), [spNum, 1], @sum, 0);
A_sum = accumarray(labels_vec, A(:), [spNum, 1], @sum, 0);
B_sum = accumarray(labels_vec, B(:), [spNum, 1], @sum, 0);

L_mean = L_sum ./ (count_vec + 0.0001);
A_mean = A_sum ./ (count_vec + 0.0001);
B_mean = B_sum ./ (count_vec + 0.0001);
sp_feature = [L_mean, A_mean, B_mean, log(double(count_vec))];

psi_sp = get_psi_sp(adj, sp_feature, para);
[psi_sorted, order_idx] = sort(psi_sp, 'ascend');
psi_sp_order = [psi_sorted, order_idx];
sp_seed = find_sp_seed(psi_sp, adj);

shortest_distance = save_shortest_distance(adj, sp_feature, para);
median_D = median(shortest_distance, 'omitnan');
if median_D == 0 || isnan(median_D)
    median_D = mean(shortest_distance, 'omitnan');
end
if median_D == 0 || isnan(median_D)
    median_D = eps;
end

target_area = (M * N) / wantedNum(1);

seed_Cluster = cell(0, 1);
ori_Cluster = num2cell((1:spNum).');
used = false(spNum, 1);

% For unused nodes this is their own area. Once a cluster is finalized,
% every member stores that cluster's total area so that merging into an
% already-used node can still be penalized correctly.
region_area = double(count_vec);

for i = 1:spNum
    self = psi_sp_order(i, :);
    search = [];
    clusterII = [];
    shortestD = [];

    if used(self(2))
        continue
    end

    link_Vector = adj(self(2), :);

    while true
        search(end+1) = self(2); %#ok<AGROW>
        used(self(2)) = true;

        if size(seed_Cluster, 1) + ...
                sum(~cellfun('isempty', ori_Cluster)) - length(search) + 1 ...
                == wantedNum(1)
            labels_turn = save_labels(search, clusterII, seed_Cluster, ...
                ori_Cluster, M, N, labels, 1);
            labelsj(:, :, end+1) = labels_turn;
            wantedNum(1) = [];
            if isempty(wantedNum)
                break
            end
            target_area = (M * N) / wantedNum(1);
        end

        if sp_seed(self(2)) == 1
            members = unique(search);
            ori_Cluster(members) = {[]};
            seed_Cluster(end+1, 1) = {members}; %#ok<AGROW>
            total_area = sum(count_vec(members));
            region_area(members) = total_area;
            break
        end

        linked_SP = find(link_Vector ~= 0);

        if numel(clusterII) <= 1
            if level < 5
                Max_innerD = median_D;
            else
                Max_innerD = (log(level) + 1) * median_D;
            end
        else
            Max_innerD = max(shortestD);
        end

        clusterII = unique([clusterII, self(2)]);
        linked_SP_psi = psi_sp(linked_SP);
        higherN = linked_SP(linked_SP_psi > self(1));

        if isempty(higherN)
            members = unique(search);
            ori_Cluster(members) = {[]};
            seed_Cluster(end+1, 1) = {members}; %#ok<AGROW>
            total_area = sum(count_vec(members));
            region_area(members) = total_area;
            break
        end

        self_LAB = sp_feature(self(2), 1:3);
        higherN_LAB = sp_feature(higherN, 1:3);
        color_D = sqrt((higherN_LAB - self_LAB).^2 * para(:));

        growing_area = sum(count_vec(unique(search)));
        merged_ratio = (growing_area + region_area(higherN)) / target_area;
        size_penalty = max(0, merged_ratio - size_tau).^2;
        effective_D = color_D + lambda_size * median_D * size_penalty;

        [min_effective_D, father_pos] = min(effective_D);
        father = higherN(father_pos);
        min_color_D = color_D(father_pos);

        if min_effective_D > D_gap * Max_innerD
            members = unique(search);
            seed_Cluster(end+1, 1) = {members}; %#ok<AGROW>
            ori_Cluster(members) = {[]};
            total_area = sum(count_vec(members));
            region_area(members) = total_area;
            break
        elseif used(father)
            for ii = 1:length(seed_Cluster)
                if any(seed_Cluster{ii} == father)
                    members = unique([seed_Cluster{ii}, clusterII]);
                    seed_Cluster{ii} = members;
                    ori_Cluster(clusterII) = {[]};
                    total_area = sum(count_vec(members));
                    region_area(members) = total_area;
                    break
                end
            end
            break
        else
            clusterII = [clusterII, father]; %#ok<AGROW>
            shortestD = [shortestD, min_color_D]; %#ok<AGROW>

            current_label = self(2);
            link_Vector(father) = 0;
            highlink_V = adj(father, :);
            highlink_V(current_label) = 0;
            link_Vector = link_Vector + highlink_V;

            self(2) = father;
            self(1) = psi_sp(father);
        end
    end

    if isempty(wantedNum)
        break
    end

    if size(seed_Cluster, 1) + sum(~cellfun('isempty', ori_Cluster)) ...
            <= wantedNum(1)
        labels_turn = save_labels(search, clusterII, seed_Cluster, ...
            ori_Cluster, M, N, labels, 0);
        labelsj(:, :, end+1) = labels_turn;
        wantedNum(1) = [];
        if isempty(wantedNum)
            break
        end
        target_area = (M * N) / wantedNum(1);
    end
end

labels = save_labels(search, clusterII, seed_Cluster, ...
    ori_Cluster, M, N, labels, 0);
end

function labels_turn = save_labels(search, clusterII, seed_Cluster, ...
    ori_Cluster, M, N, labels, kind)
if kind == 1
    seed_Cluster(end+1, 1) = {search};
    ori_Cluster(clusterII) = {[]};
end
new_Cluster = ori_Cluster(~cellfun('isempty', ori_Cluster));
seed_Cluster = vertcat(seed_Cluster, new_Cluster);
for ii = 1:numel(seed_Cluster)
    seed_Cluster{ii} = unique(seed_Cluster{ii});
end
map = zeros(max(labels(:)), 1, 'uint32');
for k = 1:numel(seed_Cluster)
    map(seed_Cluster{k}) = k;
end
labels_turn = reshape(map(labels(:)), M, N);
labels_turn = double(labels_turn);
end

function psi_sp = get_psi_sp(adj, sp_feature, para)
spNum = size(sp_feature, 1);
[row_idx, col_idx, w] = find(adj);
if isempty(row_idx)
    psi_sp = zeros(spNum, 1);
    return
end
diffLAB = sp_feature(row_idx, 1:3) - sp_feature(col_idx, 1:3);
D = sqrt(sum((diffLAB.^2) .* ((para(:)').^2), 2));
boundary_sum = full(sum(adj, 2));
neighbor_num = full(sum(adj ~= 0, 2));
weighted_D = D .* (w ./ boundary_sum(row_idx));
D_final = accumarray(row_idx, weighted_D, [spNum, 1], @sum, 0);
psi_sp = zeros(spNum, 1);
valid = (neighbor_num > 0) & (D_final > 0);
psi_sp(valid) = neighbor_num(valid) ./ D_final(valid);
end

function sp_seed = find_sp_seed(psi_sp, adj)
spNum = numel(psi_sp);
[row_idx, col_idx] = find(adj);
sp_seed = ones(spNum, 1);
if isempty(row_idx)
    return
end
has_higher_neighbor = psi_sp(col_idx) > psi_sp(row_idx);
sp_seed(unique(row_idx(has_higher_neighbor))) = 0;
end

function shortestDistance = save_shortest_distance( ...
    linkMatrix, SP_featureMatrix, para)
SPNumber = size(SP_featureMatrix, 1);
[row_idx, col_idx] = find(linkMatrix);
if isempty(row_idx)
    shortestDistance = nan(SPNumber, 1);
    return
end
diffLAB = SP_featureMatrix(row_idx, 1:3) - ...
    SP_featureMatrix(col_idx, 1:3);
D = sqrt(para(1) * (diffLAB(:, 1).^2) + ...
    para(2) * (diffLAB(:, 2).^2) + ...
    para(3) * (diffLAB(:, 3).^2));
shortestDistance = accumarray(row_idx, D, [SPNumber, 1], @min, NaN);
end
