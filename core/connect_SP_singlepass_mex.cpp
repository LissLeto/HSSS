#include "mex.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace {

struct Component {
    bool active = true;
    int label = 0;
    mwIndex min_index = 0;
    std::vector<mwIndex> pixels;
};

class DisjointSet {
public:
    explicit DisjointSet(mwSize count) : parent_(count), size_(count, 1) {
        for (mwIndex i = 0; i < count; ++i) {
            parent_[i] = i;
        }
    }

    mwIndex find(mwIndex value) {
        mwIndex root = value;
        while (parent_[root] != root) {
            root = parent_[root];
        }
        while (parent_[value] != value) {
            const mwIndex next = parent_[value];
            parent_[value] = root;
            value = next;
        }
        return root;
    }

    void unite(mwIndex left, mwIndex right) {
        left = find(left);
        right = find(right);
        if (left == right) {
            return;
        }
        if (size_[left] < size_[right]) {
            std::swap(left, right);
        }
        parent_[right] = left;
        size_[left] += size_[right];
    }

private:
    std::vector<mwIndex> parent_;
    std::vector<mwSize> size_;
};

void validate_inputs(int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("HSSS:connectSinglepass:nrhs",
            "Usage: labels_out = connect_SP_singlepass_mex(labels, LAB, para)");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("HSSS:connectSinglepass:labels",
            "labels must be a real double 2-D matrix.");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
        mxGetNumberOfDimensions(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("HSSS:connectSinglepass:LAB",
            "LAB must be a real double M-by-N-by-3 array.");
    }
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        mxGetNumberOfElements(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("HSSS:connectSinglepass:para",
            "para must contain three real double values.");
    }
}

}  // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    validate_inputs(nrhs, prhs);
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("HSSS:connectSinglepass:nlhs",
            "One output is supported.");
    }

    const mwSize* label_dims = mxGetDimensions(prhs[0]);
    const mwSize* lab_dims = mxGetDimensions(prhs[1]);
    const mwSize rows = label_dims[0];
    const mwSize cols = label_dims[1];
    const mwSize pixel_count = rows * cols;
    if (lab_dims[0] != rows || lab_dims[1] != cols || lab_dims[2] != 3) {
        mexErrMsgIdAndTxt("HSSS:connectSinglepass:size",
            "LAB must have the same height and width as labels.");
    }

    const double* labels_in = mxGetPr(prhs[0]);
    const double* lab = mxGetPr(prhs[1]);
    const double* para = mxGetPr(prhs[2]);
    const double* L = lab;
    const double* A = lab + pixel_count;
    const double* B = lab + 2 * pixel_count;

    double min_label = std::numeric_limits<double>::infinity();
    double max_label = -std::numeric_limits<double>::infinity();
    for (mwIndex p = 0; p < pixel_count; ++p) {
        if (!std::isfinite(labels_in[p]) || labels_in[p] < 0.0 ||
            std::floor(labels_in[p]) != labels_in[p]) {
            mexErrMsgIdAndTxt("HSSS:connectSinglepass:labelValues",
                "labels must contain finite nonnegative integers.");
        }
        min_label = std::min(min_label, labels_in[p]);
        max_label = std::max(max_label, labels_in[p]);
    }

    const int shift = (min_label == 0.0) ? 1 : 0;
    const int start_label = shift ? 2 : 1;
    const int label_count = static_cast<int>(max_label) + shift;
    std::vector<int> current_label(pixel_count);
    for (mwIndex p = 0; p < pixel_count; ++p) {
        current_label[p] = static_cast<int>(labels_in[p]) + shift;
    }

    std::vector<double> sum_L(label_count + 1, 0.0);
    std::vector<double> sum_A(label_count + 1, 0.0);
    std::vector<double> sum_B(label_count + 1, 0.0);
    std::vector<mwSize> label_sizes(label_count + 1, 0);
    for (mwIndex p = 0; p < pixel_count; ++p) {
        const int label = current_label[p];
        sum_L[label] += L[p];
        sum_A[label] += A[p];
        sum_B[label] += B[p];
        ++label_sizes[label];
    }
    std::vector<double> mean_L(label_count + 1,
        std::numeric_limits<double>::quiet_NaN());
    std::vector<double> mean_A = mean_L;
    std::vector<double> mean_B = mean_L;
    for (int label = 1; label <= label_count; ++label) {
        if (label_sizes[label] > 0) {
            const double denom = static_cast<double>(label_sizes[label]);
            mean_L[label] = sum_L[label] / denom;
            mean_A[label] = sum_A[label] / denom;
            mean_B[label] = sum_B[label] / denom;
        }
    }

    // Build every initial 8-connected component in one image scan.
    DisjointSet dsu(pixel_count);
    for (mwIndex col = 0; col < cols; ++col) {
        for (mwIndex row = 0; row < rows; ++row) {
            const mwIndex p = row + col * rows;
            const int label = current_label[p];
            if (row > 0 && current_label[p - 1] == label) {
                dsu.unite(p, p - 1);
            }
            if (col > 0) {
                const mwIndex left = p - rows;
                if (current_label[left] == label) {
                    dsu.unite(p, left);
                }
                if (row > 0 && current_label[left - 1] == label) {
                    dsu.unite(p, left - 1);
                }
                if (row + 1 < rows && current_label[left + 1] == label) {
                    dsu.unite(p, left + 1);
                }
            }
        }
    }

    std::vector<int> root_to_component(pixel_count, -1);
    std::vector<int> pixel_component(pixel_count, -1);
    std::vector<Component> components;
    components.reserve(pixel_count / 4 + 1);
    for (mwIndex p = 0; p < pixel_count; ++p) {
        const mwIndex root = dsu.find(p);
        int component_id = root_to_component[root];
        if (component_id < 0) {
            component_id = static_cast<int>(components.size());
            root_to_component[root] = component_id;
            Component component;
            component.label = current_label[p];
            component.min_index = p;
            components.push_back(std::move(component));
        }
        components[component_id].pixels.push_back(p);
        pixel_component[p] = component_id;
    }

    using OrderedComponent = std::pair<mwIndex, int>;
    std::vector<std::set<OrderedComponent>> components_by_label(label_count + 1);
    for (int component_id = 0;
         component_id < static_cast<int>(components.size()); ++component_id) {
        const Component& component = components[component_id];
        components_by_label[component.label].insert(
            {component.min_index, component_id});
    }

    auto visit_neighbors = [rows, cols](mwIndex p, const auto& visitor) {
        const mwIndex row = p % rows;
        const mwIndex col = p / rows;
        const mwIndex row_begin = (row == 0) ? 0 : row - 1;
        const mwIndex row_end = std::min(rows - 1, row + 1);
        const mwIndex col_begin = (col == 0) ? 0 : col - 1;
        const mwIndex col_end = std::min(cols - 1, col + 1);
        for (mwIndex nc = col_begin; nc <= col_end; ++nc) {
            for (mwIndex nr = row_begin; nr <= row_end; ++nr) {
                const mwIndex q = nr + nc * rows;
                if (q != p) {
                    visitor(q);
                }
            }
        }
    };

    auto relabel_component = [&](int component_id, int new_label) {
        Component& component = components[component_id];
        const int old_label = component.label;
        components_by_label[old_label].erase(
            {component.min_index, component_id});

        std::set<int> adjacent_components;
        for (mwIndex p : component.pixels) {
            visit_neighbors(p, [&](mwIndex q) {
                if (current_label[q] == new_label) {
                    const int adjacent_id = pixel_component[q];
                    if (adjacent_id != component_id &&
                        components[adjacent_id].active) {
                        adjacent_components.insert(adjacent_id);
                    }
                }
            });
        }

        for (mwIndex p : component.pixels) {
            current_label[p] = new_label;
        }
        component.label = new_label;

        std::vector<int> merge_ids;
        merge_ids.reserve(adjacent_components.size() + 1);
        merge_ids.push_back(component_id);
        for (int adjacent_id : adjacent_components) {
            components_by_label[new_label].erase(
                {components[adjacent_id].min_index, adjacent_id});
            merge_ids.push_back(adjacent_id);
        }

        int base_id = merge_ids.front();
        for (int id : merge_ids) {
            if (components[id].pixels.size() > components[base_id].pixels.size()) {
                base_id = id;
            }
        }
        Component& base = components[base_id];
        base.label = new_label;
        for (int id : merge_ids) {
            if (id == base_id) {
                continue;
            }
            Component& other = components[id];
            base.min_index = std::min(base.min_index, other.min_index);
            for (mwIndex p : other.pixels) {
                base.pixels.push_back(p);
                pixel_component[p] = base_id;
            }
            other.pixels.clear();
            other.active = false;
        }
        // Keep a deterministic linear-index order after dynamic unions.
        std::sort(base.pixels.begin(), base.pixels.end());
        for (mwIndex p : base.pixels) {
            pixel_component[p] = base_id;
        }
        components_by_label[new_label].insert({base.min_index, base_id});
    };

    for (int label = start_label; label <= label_count; ++label) {
        std::vector<int> label_components;
        label_components.reserve(components_by_label[label].size());
        for (const OrderedComponent& item : components_by_label[label]) {
            label_components.push_back(item.second);
        }
        if (label_components.empty()) {
            continue;
        }

        if (label_components.size() == 1) {
            const int component_id = label_components.front();
            const Component& component = components[component_id];
            if (component.pixels.size() == 1) {
                const mwIndex p = component.pixels.front();
                const mwIndex row = p % rows;
                const mwIndex col = p / rows;
                int new_label = label;
                if (row + 1 < rows) {
                    new_label = current_label[p + 1];
                } else if (row > 0) {
                    new_label = current_label[p - 1];
                } else if (col + 1 < cols) {
                    new_label = current_label[p + rows];
                } else if (col > 0) {
                    new_label = current_label[p - rows];
                }
                if (new_label != label) {
                    relabel_component(component_id, new_label);
                }
            }
            continue;
        }

        mwSize max_area = 0;
        mwIndex max_position = 0;
        for (mwIndex pos = 0; pos < label_components.size(); ++pos) {
            const mwSize area = components[label_components[pos]].pixels.size();
            if (area > max_area) {
                max_area = area;
                max_position = pos;
            }
        }

        for (mwIndex pos = 0; pos < label_components.size(); ++pos) {
            if (pos == max_position && max_area > 8) {
                continue;
            }
            const int component_id = label_components[pos];
            const Component& component = components[component_id];
            if (!component.active || component.label != label) {
                continue;
            }

            std::set<int> neighbor_labels;
            double component_L = 0.0;
            double component_A = 0.0;
            double component_B = 0.0;
            for (mwIndex p : component.pixels) {
                component_L += L[p];
                component_A += A[p];
                component_B += B[p];
                visit_neighbors(p, [&](mwIndex q) {
                    const int neighbor_label = current_label[q];
                    if (neighbor_label != 0 && neighbor_label != label) {
                        neighbor_labels.insert(neighbor_label);
                    }
                });
            }
            if (neighbor_labels.empty()) {
                continue;
            }

            const double denom = static_cast<double>(component.pixels.size());
            component_L /= denom;
            component_A /= denom;
            component_B /= denom;
            double best_distance = std::numeric_limits<double>::infinity();
            int best_label = -1;
            for (int neighbor_label : neighbor_labels) {
                const double dL = mean_L[neighbor_label] - component_L;
                const double dA = mean_A[neighbor_label] - component_A;
                const double dB = mean_B[neighbor_label] - component_B;
                const double distance = dL * dL * para[0] +
                    dA * dA * para[1] + dB * dB * para[2];
                if (distance < best_distance) {
                    best_distance = distance;
                    best_label = neighbor_label;
                }
            }
            if (best_label > 0) {
                relabel_component(component_id, best_label);
            }
        }
    }

    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double* labels_out = mxGetPr(plhs[0]);
    for (mwIndex p = 0; p < pixel_count; ++p) {
        labels_out[p] = static_cast<double>(current_label[p] - shift);
    }
}
