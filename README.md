# Functional-Edged Network Modeling

This repository contains the MATLAB implementation accompanying the paper
**"Functional-Edged Network Modeling"**. It provides an end-to-end simulation
workflow for generating an irregularly observed functional-edged network,
estimating the proposed Functional-Edged Network (FEN) model by Riemannian
conjugate gradient optimization, and evaluating the fitted tensor and network
statistics. Running `main.m` produces the seven goodness-of-fit and residual
diagnostic panels underlying Figure 5.

## Scope of this release

The included driver focuses on the complete simulation workflow used for
Figure 5. The same model-fitting functions under `src/` implement the functional
Tucker model and optimization algorithm described in Sections 3 and 4 of the
paper. No external dataset is required for this simulation: the latent tensor,
irregular observation mask, and noisy observations are generated within each
Monte Carlo replication.

The current values in `main.m` use a reduced network size so users can check
their MATLAB installation more quickly. To run the large-scale Figure 5 setting,
change only `m`, `s1`, and `K` as described below.

## Repository structure

```text
FEN-code/
|-- main.m                         User-facing end-to-end driver
|-- README.md                      Documentation and reproduction instructions
|-- src/                           Backend FEN implementation
|   |-- data_generation/           Synthetic tensor and mask construction
|   |-- optimization/              FEN objective and Riemannian optimization
|   `-- evaluation/                Factor alignment and fit diagnostics
|-- external/
|   `-- tensor_toolbox/            Tensor Toolbox for MATLAB 3.2.1
`-- results/                       Outputs created by main.m
    `-- README.md                  Description of generated result files
```

Only `main.m` is intended to be executed directly. The script automatically
adds `src/` and the bundled Tensor Toolbox to the MATLAB path. Users therefore
do not need to run the backend files individually or manually configure paths.

## Software requirements

- MATLAB R2022a or later is recommended.
- Statistics and Machine Learning Toolbox, used by:
  - `unidrnd` and `normrnd` for data generation;
  - `boxplot` for Figure 5 panels;
  - `lillietest` for the residual normality diagnostic.
- Econometrics Toolbox, used by `lbqtest` for the Ljung-Box residual diagnostic.
- Tensor Toolbox for MATLAB 3.2.1. A compatible copy is bundled in
  `external/tensor_toolbox/` and is added to the MATLAB path by `main.m`.

The bundled Tensor Toolbox is distributed under its own 2-clause BSD license;
see `external/tensor_toolbox/LICENSE.txt` and its accompanying README.

## Quick start

1. Download and extract the repository without changing its internal folder
   structure.
2. Start MATLAB and set the current folder to the extracted repository root,
   i.e., the folder containing `main.m`.
3. Run:

   ```matlab
   main
   ```

4. The script prints the Monte Carlo mean and standard deviation of MSE and
   R-squared in the Command Window, opens seven diagnostic figures, and writes
   all reproducibility outputs to `results/`.

The full setting performs repeated tensor decompositions and may require
substantial computation time. For a quick installation check, users can set
`rep = 1` near the top of `main.m`; this changes only the number of Monte Carlo
replications.

## Reproducing Figure 5

Figure 5 uses the large-scale setting with 40% missing observations and noise
variance 0.2. In the user-adjustable parameter block of `main.m`, set:

```matlab
rep = 20;
m = 50;
n = m;
s1 = 15;
s2 = s1;
K = 25;
sigma = sqrt(0.2);
w = 0.4;
```

The remaining defaults already define 100 functional grid points on
`[-0.98,1]`, a smoothness coefficient of 0.1 for each functional component,
and the optimization controls used by the supplied implementation. In
particular:

```matlab
t_start = -0.98;
t_end = 1;
seq = 0.02;
alpha = 0.1 * ones(1,K);
max_m = 5;
tol = 1e-6;
maxiter = 1000;
```


The seven generated panels compare estimated and true values for:

1. reciprocity, corresponding to theta_R;
2. cyclic triads, corresponding to theta_CT;
3. in-two-stars, corresponding to theta_ITS;
4. out-two-stars, corresponding to theta_OTS;
5. overall edge density/intensity, corresponding to theta_E;
6. the Lilliefors residual-normality p-value, p_normal;
7. the lag-20 Ljung-Box residual diagnostic p-value, p_ind.

The MATLAB figures contain the underlying simulation output. Publication
typesetting, font selection, and final multi-panel assembly may be adjusted
separately without rerunning the model.

## End-to-end workflow

Each Monte Carlo replication follows the steps below.

### 1. Generate the latent functional adjacency tensor

`generate_data_true.m` constructs the functional Tucker model

```text
X = B x_1 Phi x_2 Phi x_3 G,
```

which corresponds to Equation (3) of the paper. The two node modes share the
same factor `Phi`, representing the community-related node basis. The third
factor contains Fourier basis functions on the functional domain.

The supporting functions are:

- `generate_basis_discrete.m`: generates the orthonormal node basis `Phi`;
- `generate_basis_continue.m`: generates the continuous Fourier functions;
- `generate_core.m`: generates the full-rank slices of the Tucker core `B`;
- `generate_data_true_funhandle.m`: evaluates an individual functional edge;
- `mu_true.m`: specifies the true mean function, currently zero.

### 2. Generate noisy, irregular observations

`generate_data_observe.m` evaluates all functional edges on the common grid,
adds independent Gaussian noise, and independently samples observation points
for every ordered node pair. It returns:

- `X_observe`: the observed tensor with zeros at missing entries;
- `Omiga`: the binary observation mask (the spelling is retained for backward
  compatibility with the original implementation);
- `X_observe_full`: the complete tensor including simulated noise;
- `X_observe_full_nornd`: the complete noiseless latent tensor;
- `seqtime`: the functional observation grid.

This implements the observation and masking construction in Equations (7) and
(8). With the default `w = 0.4`, each edge is observed at 60% of the grid points.

### 3. Initialize the FEN decomposition

`Initialization_tucker.m` uses the symmetry-preserving higher-order SVD
initialization. One common node factor is used in modes 1 and 2, so
`U{1} = U{2} = Phi`. The functional factor is initialized using the dominant
mode-3 singular vectors.

### 4. Optimize the FEN objective

The fitted tensor minimizes the masked least-squares loss plus the functional
smoothness penalty in Equations (8) and (13). The optimizer follows Algorithm 1:

- `target_f.m`: evaluates the objective function;
- `grad_Euc.m`: computes the Euclidean gradient (Equation (14));
- `Euc_to_Rie.m`: projects onto the shared-Phi tangent space (Equation (15));
- `vector_transport.m`: transports the previous direction (Equation (16));
- `beta.m`: computes the nonnegative Polak-Ribiere coefficient;
- `find_alpha.m`: proposes the step size;
- `find_m.m`: performs Armijo backtracking;
- `Retraction.m`: applies the symmetry-HOSVD retraction in Algorithm 2;
- `conjugate_dradient.m`: coordinates the Riemannian conjugate gradient
  iterations. The historical filename is retained to avoid changing the
  numerical implementation and call interface.

`dif_mat.m` creates the first-difference penalty matrix, while `innerP.m`
implements the tensor inner product used by the optimization routines.

### 5. Align factors and evaluate the fit

`adjust.m` applies orthogonal Procrustes rotations only to make estimated and
true factors comparable for reporting. This rotation does not alter the fitted
tensor. The prediction used for evaluation is the direct optimizer output
`X_next`.

`evlua.m` evaluates the estimated tensor against the noiseless truth. It
computes:

- MSE on the entries held out by `Omiga`;
- R-squared over the complete noiseless tensor;
- estimated and true reciprocity;
- estimated and true cyclic-triad strength;
- estimated and true in-two-star strength;
- estimated and true out-two-star strength;
- estimated and true overall edge density/intensity.

The residual diagnostics in `main.m` use the complete fitted residual

```matlab
res = double(X_output - X_observe_full);
res = res(:);
```

and report the Lilliefors and lag-20 Ljung-Box p-values. The former examines
the marginal normality assumption; the latter is used as a serial-correlation
diagnostic. As stated in the paper, a large p-value means that the diagnostic
does not provide evidence against the corresponding working assumption; it is
not interpreted as proof that the assumption holds.

## Function reference

| File | Purpose | Called by |
|---|---|---|
| `main.m` | End-to-end simulation, fitting, evaluation, plotting, and output saving | User |
| `generate_data_true.m` | Builds latent functional edges from the functional Tucker factors | `main.m` |
| `generate_basis_discrete.m` | Generates the orthonormal node/community basis | `generate_data_true.m` |
| `generate_basis_continue.m` | Defines the continuous Fourier basis functions | `generate_data_true.m` |
| `generate_core.m` | Generates full-rank slices of the true Tucker core | `generate_data_true.m` |
| `generate_data_true_funhandle.m` | Evaluates a functional edge from its basis coefficients | `generate_data_true.m` |
| `mu_true.m` | Defines the simulation mean function, currently zero | `generate_data_true_funhandle.m` |
| `generate_data_observe.m` | Adds noise and irregular missingness | `main.m` |
| `mu_hat.m` | Optional local-linear mean estimator for nonzero unknown means; not used in the centered Figure 5 simulation | Optional preprocessing |
| `kappa_1D.m` | Epanechnikov kernel used by `mu_hat.m` | `mu_hat.m` |
| `Initialization_tucker.m` | Computes the shared-Phi initial point | `main.m` |
| `conjugate_dradient.m` | Runs Riemannian conjugate gradient iterations | `main.m` |
| `grad_Euc.m` | Euclidean gradient of data-fit and smoothness terms | Optimizer |
| `Euc_to_Rie.m` | Shared-Phi tangent projection | Optimizer |
| `vector_transport.m` | Projection-based vector transport | Optimizer |
| `Retraction.m` | Symmetry-HOSVD fixed-rank retraction | Optimizer |
| `target_f.m` | Objective value | Optimizer, `main.m` |
| `find_alpha.m` | Closed-form proposal for the line-search step size | Optimizer |
| `find_m.m` | Armijo backtracking search | Optimizer |
| `beta.m` | Conjugate-direction coefficient | Optimizer |
| `innerP.m` | Frobenius inner product for two tensors | Optimizer |
| `dif_mat.m` | First-difference penalty matrix for functional smoothness | `main.m` |
| `adjust.m` | Orthogonal factor alignment for reporting | `main.m` |
| `evlua.m` | Prediction and network goodness-of-fit metrics | `main.m` |

## Inputs and parameter meanings

The Figure 5 workflow does not read external files. All inputs are specified in
the parameter block at the top of `main.m`.

| Variable | Meaning |
|---|---|
| `rep` | Number of independent Monte Carlo replications |
| `m`, `n` | Numbers of nodes in the two adjacency modes; the FEN setting uses `n=m` |
| `s1`, `s2` | Node-mode Tucker ranks; the shared-Phi model uses `s1=s2` |
| `K` | Functional-mode Tucker rank |
| `sigma` | Standard deviation of Gaussian observation noise |
| `w` | Proportion of functional grid entries treated as missing |
| `N` | Scale bound used when generating the Tucker core; not the number of samples |
| `t_start`, `t_end`, `seq` | Functional grid endpoints and spacing |
| `alpha` | Smoothness penalty coefficients for the functional basis vectors |
| `max_m` | Maximum number of backtracking reductions |
| `tol` | Absolute objective-change convergence tolerance |
| `maxiter` | Maximum number of conjugate-gradient iterations |

## Output files and workspace variables

Every run creates or overwrites the following files under `results/`:

- `figure5_results.mat`: contains `myoutput`, `settings`, and `summary`;
- `run_summary.txt`: records the dimensions, rank, missing rate, noise
  variance, and aggregate MSE/R-squared;
- `reciprocity.fig/.png`;
- `cyclic_triads.fig/.png`;
- `in_two_stars.fig/.png`;
- `out_two_stars.fig/.png`;
- `edge_density.fig/.png`;
- `residual_normality.fig/.png`;
- `residual_independence.fig/.png`.

The `.fig` files retain editable MATLAB graphics, while the `.png` files are
saved at 300 dpi for convenient inspection.

`myoutput(tt)` stores the following record for replication `tt`:

| Field | Contents |
|---|---|
| `G_true`, `U_true` | True Tucker core and factors |
| `Omiga` | Binary irregular-observation mask |
| `X_observe_full` | Complete noisy tensor |
| `X_observe_full_nornd` | Complete noiseless tensor |
| `X_observe_central` | Masked tensor supplied to the optimizer |
| `G_output`, `U_output` | Estimated factors after reporting alignment |
| `X_output` | Final fitted tensor returned by the optimizer |
| `MSE`, `R2` | Tensor prediction metrics |
| `Rec_esti`, `Rec_true` | Estimated and true reciprocity |
| `Cyc_esti`, `Cyc_true` | Estimated and true cyclic-triad strength |
| `InTwo_esti`, `InTwo_true` | Estimated and true in-two-star strength |
| `OutTwo_esti`, `OutTwo_true` | Estimated and true out-two-star strength |
| `Densi_esti`, `Densi_true` | Estimated and true edge density/intensity |
| `normal_p` | Lilliefors residual-normality p-value |
| `LjungBox` | Lag-20 Ljung-Box residual diagnostic p-value |
| `score_list` | Objective values recorded during optimization |

## Using the optimizer with another tensor

The simulation driver represents missing entries by zero and provides a binary
mask of the same size. To fit another three-dimensional tensor on a common
functional grid:

1. construct a Tensor Toolbox `tensor` object `X_observe_central` with zeros at
   unobserved entries;
2. construct the matching 0/1 Tensor Toolbox object `Omiga`;
3. set `r = [s,s,K]`, the difference matrix `Sigma = dif_mat(L)`, and the
   smoothness coefficients `alpha`;
4. follow the initialization and optimization calls in lines corresponding to
   Steps 3 and 4 of `main.m`.

The input must have equal first and second dimensions, and the first two ranks
must be equal because the FEN model uses one shared node factor in both modes.

## Troubleshooting

- **`tensor` or `ttm` is undefined:** run `main.m` from the repository root.
  The script adds `external/tensor_toolbox` automatically; do not move that
  folder independently.
- **`lillietest`, `normrnd`, or `boxplot` is undefined:** install or enable the
  Statistics and Machine Learning Toolbox.
- **`lbqtest` is undefined:** install or enable the Econometrics Toolbox.
- **A full Figure 5 run is slow:** first verify the installation with `rep=1`
  and the reduced dimensions, then restore `rep=20`, `m=50`, `s1=15`, and
  `K=25` for the reported setting.
- **Results differ after changing parameters:** confirm the parameter values
  stored in `results/figure5_results.mat` under `settings`. 
## Code availability

The implementation is organized so that the root directory contains only the
user-facing driver and documentation. Backend functions are separated by role,
and third-party software is isolated under `external/`, making the execution
path and the connection between the manuscript and implementation explicit.


Finally, `fitting_results_of_RSSI.zip` contains the full set of fitted results for the industrial IoT RSSI case study, including model outputs and intermediate fitting quantities used in the corresponding analysis.

