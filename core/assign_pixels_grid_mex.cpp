#include "mex.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

/*
 * newLabels = assign_pixels_grid_mex( ...
 *      L, A, B, ...
 *      seed_rows, seed_cols, seeds_L, seeds_A, seeds_B, ...
 *      nonSeedIdx, search_R, weight, M, N)
 *
 * Inputs:
 *   L, A, B       : MxN double
 *   seed_rows     : nSeed x 1 double, 1-based
 *   seed_cols     : nSeed x 1 double, 1-based
 *   seeds_L/A/B   : nSeed x 1 double
 *   nonSeedIdx    : nNonSeed x 1 double, 1-based linear indices
 *   search_R      : scalar double
 *   weight        : 5x1 double = [wS^2, wS^2, para1, para2, para3]'
 *   M, N          : scalar double
 *
 * Output:
 *   newLabels     : nNonSeed x 1 double
 *
 * Notes:
 *   Output label is the 1-based seed index.
 *   If no candidate seed is found, output 0.
 */

static void checkInputs(int nrhs, const mxArray *prhs[])
{
    if (nrhs != 13) {
        mexErrMsgIdAndTxt("assign_pixels_grid_mex:nrhs",
                          "Expected 13 input arguments.");
    }

    for (int i = 0; i < 13; ++i) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt("assign_pixels_grid_mex:type",
                              "All inputs must be real double.");
        }
    }
}
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    checkInputs(nrhs, prhs);

    const mxArray *L_in         = prhs[0];
    const mxArray *A_in         = prhs[1];
    const mxArray *B_in         = prhs[2];
    const mxArray *seed_rows_in = prhs[3];
    const mxArray *seed_cols_in = prhs[4];
    const mxArray *seeds_L_in   = prhs[5];
    const mxArray *seeds_A_in   = prhs[6];
    const mxArray *seeds_B_in   = prhs[7];
    const mxArray *nonSeed_in   = prhs[8];
    const mxArray *searchR_in   = prhs[9];
    const mxArray *weight_in    = prhs[10];
    const mxArray *M_in         = prhs[11];
    const mxArray *N_in         = prhs[12];

    const double *L = mxGetPr(L_in);
    const double *A = mxGetPr(A_in);
    const double *B = mxGetPr(B_in);

    const double *seed_rows = mxGetPr(seed_rows_in);
    const double *seed_cols = mxGetPr(seed_cols_in);
    const double *seeds_L   = mxGetPr(seeds_L_in);
    const double *seeds_A   = mxGetPr(seeds_A_in);
    const double *seeds_B   = mxGetPr(seeds_B_in);

    const double *nonSeedIdx = mxGetPr(nonSeed_in);
    const double search_R = mxGetScalar(searchR_in);
    const double *weight = mxGetPr(weight_in);

    const mwSize M = static_cast<mwSize>(mxGetScalar(M_in));
    const mwSize N = static_cast<mwSize>(mxGetScalar(N_in));

    const mwSize nSeed = mxGetNumberOfElements(seed_rows_in);
    const mwSize nNonSeed = mxGetNumberOfElements(nonSeed_in);

    if (mxGetNumberOfElements(seed_cols_in) != nSeed ||
        mxGetNumberOfElements(seeds_L_in)   != nSeed ||
        mxGetNumberOfElements(seeds_A_in)   != nSeed ||
        mxGetNumberOfElements(seeds_B_in)   != nSeed) {
        mexErrMsgIdAndTxt("assign_pixels_grid_mex:size",
                          "seed_rows, seed_cols, seeds_L, seeds_A, seeds_B must have the same length.");
    }

    if (mxGetNumberOfElements(weight_in) != 5) {
        mexErrMsgIdAndTxt("assign_pixels_grid_mex:weight",
                          "weight must have 5 elements.");
    }

    plhs[0] = mxCreateDoubleMatrix(nNonSeed, 1, mxREAL);
    double *newLabels = mxGetPr(plhs[0]);

    if (!(search_R > 0.0) || !std::isfinite(search_R)) {
        mexErrMsgIdAndTxt("assign_pixels_grid_mex:searchRadius",
                          "search_R must be finite and positive.");
    }

    // Uniform grid in the same normalized coordinate system used by the
    // MATLAB implementation. A candidate with |delta| < search_R can only
    // be in the pixel's own grid cell or one of its eight neighboring cells.
    const mwSize grid_count =
        static_cast<mwSize>(std::floor(1.0 / search_R)) + 2;
    std::vector<std::vector<mwIndex>> grid(grid_count * grid_count);
    auto grid_coord = [grid_count, search_R](double coordinate) {
        mwIndex cell = static_cast<mwIndex>(std::floor(coordinate / search_R));
        return std::min<mwIndex>(cell, grid_count - 1);
    };
    for (mwIndex s = 0; s < nSeed; ++s) {
        const double seed_row_norm = seed_rows[s] / static_cast<double>(M);
        const double seed_col_norm = seed_cols[s] / static_cast<double>(N);
        const mwIndex grid_row = grid_coord(seed_row_norm);
        const mwIndex grid_col = grid_coord(seed_col_norm);
        grid[grid_row + grid_col * grid_count].push_back(s);
    }

    for (mwSize t = 0; t < nNonSeed; ++t) {
        mwIndex idx0 = static_cast<mwIndex>(nonSeedIdx[t] - 1.0);

        const double i = static_cast<double>((idx0 % M) + 1);
        const double j = static_cast<double>((idx0 / M) + 1);
        const double i_norm = i / static_cast<double>(M);
        const double j_norm = j / static_cast<double>(N);

        const double selfL = L[idx0];
        const double selfA = A[idx0];
        const double selfB = B[idx0];

        double bestD = std::numeric_limits<double>::infinity();
        mwIndex bestSeed = static_cast<mwIndex>(-1);

        const mwIndex pixel_grid_row = grid_coord(i_norm);
        const mwIndex pixel_grid_col = grid_coord(j_norm);
        const mwIndex row_begin = pixel_grid_row == 0 ? 0 : pixel_grid_row - 1;
        const mwIndex row_end = std::min<mwIndex>(grid_count - 1, pixel_grid_row + 1);
        const mwIndex col_begin = pixel_grid_col == 0 ? 0 : pixel_grid_col - 1;
        const mwIndex col_end = std::min<mwIndex>(grid_count - 1, pixel_grid_col + 1);
        std::vector<mwIndex> candidates;
        for (mwIndex grid_col = col_begin; grid_col <= col_end; ++grid_col) {
            for (mwIndex grid_row = row_begin; grid_row <= row_end; ++grid_row) {
                const std::vector<mwIndex>& cell =
                    grid[grid_row + grid_col * grid_count];
                candidates.insert(candidates.end(), cell.begin(), cell.end());
            }
        }
        std::sort(candidates.begin(), candidates.end());

        for (mwIndex s : candidates) {
            const double seed_row_norm = seed_rows[s] / static_cast<double>(M);
            const double seed_col_norm = seed_cols[s] / static_cast<double>(N);
            const double dr = seed_row_norm - i_norm;
            const double dc = seed_col_norm - j_norm;
            if (std::fabs(dr) >= search_R || std::fabs(dc) >= search_R) {
                continue;
            }
            const double dL = seeds_L[s] - selfL;
            const double dA = seeds_A[s] - selfA;
            const double dB = seeds_B[s] - selfB;

            const double D =
                dr * dr * weight[0] +
                dc * dc * weight[1] +
                dL * dL * weight[2] +
                dA * dA * weight[3] +
                dB * dB * weight[4];

            if (D < bestD) {
                bestD = D;
                bestSeed = s;
            }
        }

        if (bestSeed == static_cast<mwIndex>(-1)) {
            newLabels[t] = 0.0;
        } else {
            newLabels[t] = static_cast<double>(bestSeed + 1);
        }
    }
}
