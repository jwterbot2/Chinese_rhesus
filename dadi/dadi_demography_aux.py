#################################################################################
#####                                                                       #####
#####           demographic modeling in dadi  (1d)               #####
#####                                                                       #####
#################################################################################
#!/usr/bin/env python3
##### ---------- ##### ---------- ##### ---------- ##### ---------- #####
# Dadi 1d pipeline 
##### ---------- ##### ---------- ##### ---------- ##### ---------- #####
"""
Script to test dadi's basic workflow for single population demographic inference
(1d)
"""
import numpy as np
import dadi
import matplotlib.pyplot as plt
import os
output_dir = "1d_models_results"
os.makedirs(output_dir, exist_ok=True)
data_fs = dadi.Spectrum.from_file("acalahor_data_folded_full.fs")
if not data_fs.folded:
    data_fs = data_fs.fold()
pts_l = [40, 50, 60] #grid points for extrapolation
# optimization runs for each model
n_runs = 20
max_iter = 100

# define the models to test
models = ["snm", "two_epoch", "growth", "bottlegrowth", "three_epoch"]
if "exponential_growth" in dir(dadi.Demographics1D):
    models.append("exponential_growth")
results = {}

# Model 1: Standard neutral model (snm)
func = dadi.Demographics1D.snm
func_ex = dadi.Numerics.make_extrap_log_func(func)
model = func_ex([], data_fs.sample_sizes, pts_l)
ll_snm = dadi.Inference.ll_multinom(model, data_fs)
theta_snm = dadi.Inference.optimal_sfs_scaling(model, data_fs)

results["snm"] = {
    "ll": ll_snm,
    "aic": -2*ll_snm, 
    "theta": theta_snm,
    "params": {}
}
plt.figure()
dadi.Plotting.plot_1d_comp_multinom(model, data_fs)
plt.title(f"Standard neutral model\nLL={ll_snm:.2f}, AIC={-2*ll_snm:.2f}")
plt.savefig(f"{output_dir}/snm_fit.png")
plt.close()

# Model 2: Two epoch model
func = dadi.Demographics1D.two_epoch
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = [1.0, 0.1]  # nu, T
lower_bound = [0.01, 0.001]
upper_bound = [100.0, 10.0]
param_names = ["nu", "T"]
best_params = []
best_ll = -np.inf
for i in range(n_runs):
    p0_i = dadi.Misc.perturb_params(p0, fold=2, lower_bound=lower_bound, upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0_i, data_fs, func_ex, pts_l,
                                     lower_bound=lower_bound,
                                     upper_bound=upper_bound,
                                     maxiter=max_iter)
    model = func_ex(popt, data_fs.sample_sizes, pts_l)
    ll = dadi.Inference.ll_multinom(model, data_fs)
    if ll > best_ll:
        best_ll = ll
        best_params = popt
model = func_ex(best_params, data_fs.sample_sizes, pts_l)
theta = dadi.Inference.optimal_sfs_scaling(model, data_fs)
aic = 2*len(best_params) - 2*best_ll
results["two_epoch"] = {
    "ll": best_ll,
    "aic": aic,
    "theta": theta,
    "params": dict(zip(param_names, best_params))
}
# plot the model fit
plt.figure()
dadi.Plotting.plot_1d_comp_multinom(model, data_fs)
param_str = ", ".join([f"{name}={val:.4f}" for name, val in zip(param_names, best_params)])
plt.title(f"Two epoch model\nLL={best_ll:.2f}, AIC={aic:.2f}\n{param_str}")
plt.savefig(f"{output_dir}/two_epoch_fit.png")
plt.close()

# Model 3: Growth model
func = dadi.Demographics1D.growth
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = [1.0, 0.1]  # nu, T
lower_bound = [0.01, 0.001]
upper_bound = [100.0, 10.0]
param_names = ["nu", "T"]
best_params = []
best_ll = -np.inf
for i in range(n_runs):
    p0_i = dadi.Misc.perturb_params(p0, fold=2, lower_bound=lower_bound, upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0_i, data_fs, func_ex, pts_l,
                                     lower_bound=lower_bound,
                                     upper_bound=upper_bound,
                                     maxiter=max_iter)
    model = func_ex(popt, data_fs.sample_sizes, pts_l)
    ll = dadi.Inference.ll_multinom(model, data_fs)
    if ll > best_ll:
        best_ll = ll
        best_params = popt
model = func_ex(best_params, data_fs.sample_sizes, pts_l)
theta = dadi.Inference.optimal_sfs_scaling(model, data_fs)
aic = 2*len(best_params) - 2*best_ll
results["growth"] = {
    "ll": best_ll,
    "aic": aic,
    "theta": theta,
    "params": dict(zip(param_names, best_params))
}
# plot the model fit
plt.figure()
dadi.Plotting.plot_1d_comp_multinom(model, data_fs)
param_str = ", ".join([f"{name}={val:.4f}" for name, val in zip(param_names, best_params)])
plt.title(f"Growth model\nLL={best_ll:.2f}, AIC={aic:.2f}\n{param_str}")
plt.savefig(f"{output_dir}/growth_fit.png")
plt.close()

# Model 4: Bottlegrowth model
func = dadi.Demographics1D.bottlegrowth
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = [0.1, 2.0, 0.05]  # nuB, nuF, T
lower_bound = [0.01, 0.01, 0.001]
upper_bound = [1.0, 100.0, 10.0]
param_names = ["nuB", "nuF", "T"]
best_params = []
best_ll = -np.inf
for i in range(n_runs):
    p0_i = dadi.Misc.perturb_params(p0, fold=2, lower_bound=lower_bound, upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0_i, data_fs, func_ex, pts_l,
                                     lower_bound=lower_bound,
                                     upper_bound=upper_bound,
                                     maxiter=max_iter)
    model = func_ex(popt, data_fs.sample_sizes, pts_l)
    ll = dadi.Inference.ll_multinom(model, data_fs)
    if ll > best_ll:
        best_ll = ll
        best_params = popt
model = func_ex(best_params, data_fs.sample_sizes, pts_l)
theta = dadi.Inference.optimal_sfs_scaling(model, data_fs)
aic = 2*len(best_params) - 2*best_ll
results["bottlegrowth"] = {
    "ll": best_ll,
    "aic": aic,
    "theta": theta,
    "params": dict(zip(param_names, best_params))
}
# plot the model fit
plt.figure()
dadi.Plotting.plot_1d_comp_multinom(model, data_fs)
param_str = ", ".join([f"{name}={val:.4f}" for name, val in zip(param_names, best_params)])
plt.title(f"Bottlegrowth model\nLL={best_ll:.2f}, AIC={aic:.2f}\n{param_str}")
plt.savefig(f"{output_dir}/bottlegrowth_fit.png")
plt.close()

# Model 5: Three epoch model
func = dadi.Demographics1D.three_epoch
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = [0.5, 5.0, 0.1, 0.1]  # nuB, nuF, TB, TF
lower_bound = [0.01, 0.01, 0.001, 0.001]
upper_bound = [10.0, 100.0, 5.0, 5.0]
param_names = ["nuB", "nuF", "TB", "TF"]
best_params = []
best_ll = -np.inf
for i in range(n_runs):
    p0_i = dadi.Misc.perturb_params(p0, fold=2, lower_bound=lower_bound, upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0_i, data_fs, func_ex, pts_l,
                                     lower_bound=lower_bound,
                                     upper_bound=upper_bound,
                                     maxiter=max_iter)
    model = func_ex(popt, data_fs.sample_sizes, pts_l)
    ll = dadi.Inference.ll_multinom(model, data_fs)
    if ll > best_ll:
        best_ll = ll
        best_params = popt
model = func_ex(best_params, data_fs.sample_sizes, pts_l)
theta = dadi.Inference.optimal_sfs_scaling(model, data_fs)
aic = 2*len(best_params) - 2*best_ll
results["three_epoch"] = {
    "ll": best_ll,
    "aic": aic,
    "theta": theta,
    "params": dict(zip(param_names, best_params))
}
# plot the model fit
plt.figure()
dadi.Plotting.plot_1d_comp_multinom(model, data_fs)
param_str = ", ".join([f"{name}={val:.4f}" for name, val in zip(param_names, best_params)])
plt.title(f"Three epoch model\nLL={best_ll:.2f}, AIC={aic:.2f}\n{param_str}")
plt.savefig(f"{output_dir}/three_epoch_fit.png")
plt.close()

# check and run exponential_growth model if available 
# (this avoids crashes in previous versions of dadi)
if "exponential_growth" in dir(dadi.Demographics1D):
    func = dadi.Demographics1D.exponential_growth
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    p0 = [0.1, 0.1, 0.1]  # nu0, nuF, T
    lower_bound = [0.001, 0.001, 0.001]
    upper_bound = [10.0, 100.0, 10.0]
    param_names = ["nu0", "nuF", "T"]
    best_params = []
    best_ll = -np.inf
    for i in range(n_runs):
        p0_i = dadi.Misc.perturb_params(p0, fold=2, lower_bound=lower_bound, upper_bound=upper_bound)
        popt = dadi.Inference.optimize_log(p0_i, data_fs, func_ex, pts_l,
                                         lower_bound=lower_bound,
                                         upper_bound=upper_bound,
                                         maxiter=max_iter)
        model = func_ex(popt, data_fs.sample_sizes, pts_l)
        ll = dadi.Inference.ll_multinom(model, data_fs)
        if ll > best_ll:
            best_ll = ll
            best_params = popt
    model = func_ex(best_params, data_fs.sample_sizes, pts_l)
    theta = dadi.Inference.optimal_sfs_scaling(model, data_fs)
    aic = 2*len(best_params) - 2*best_ll
    results["exponential_growth"] = {
        "ll": best_ll,
        "aic": aic,
        "theta": theta,
        "params": dict(zip(param_names, best_params))
    }
    # plot the model fit
    plt.figure()
    dadi.Plotting.plot_1d_comp_multinom(model, data_fs)
    param_str = ", ".join([f"{name}={val:.4f}" for name, val in zip(param_names, best_params)])
    plt.title(f"Exponential growth model\nLL={best_ll:.2f}, AIC={aic:.2f}\n{param_str}")
    plt.savefig(f"{output_dir}/exponential_growth_fit.png")
    plt.close()
for model_name, model_results in results.items():
    with open(f"{output_dir}/{model_name}_results.txt", 'w') as f:
        f.write(f"Model: {model_name}\n")
        f.write(f"Log-likelihood: {model_results['ll']:.4f}\n")
        f.write(f"AIC: {model_results['aic']:.4f}\n")
        f.write(f"Optimal theta: {model_results['theta']:.4f}\n")
        f.write("Parameters:\n")
        for param, value in model_results['params'].items():
            f.write(f"{param}: {value:.4f}\n")

### ---------- comparing 1d models 
# Compare models by AIC
sorted_results = sorted(results.items(), key=lambda x: x[1]['aic'])
# Generate summary table
with open(f"{output_dir}/model_comparison.txt", 'w') as f:
    f.write("Model\tLog-likelihood\tAIC\tTheta\tParameters\n")
    for model_name, model_results in sorted_results:
        param_str = ", ".join([f"{k}={v:.4f}" for k, v in model_results['params'].items()])
        f.write(f"{model_name}\t{model_results['ll']:.4f}\t{model_results['aic']:.4f}\t{model_results['theta']:.4f}\t{param_str}\n")
# Plot model comparison
plt.figure(figsize=(10, 6))
models = [x[0] for x in sorted_results]
aics = [x[1]['aic'] for x in sorted_results]
plt.bar(models, aics)
plt.xticks(rotation=45, ha="right")
plt.ylabel("AIC")
plt.title("Model Comparison by AIC")
plt.tight_layout()
plt.savefig(f"{output_dir}/model_comparison.png")
print("Analysis complete. Results saved to:", output_dir)