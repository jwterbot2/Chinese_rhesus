##### ---------- ##### ---------- ##### ---------- ##### ---------- #####
# Dadi pipeline 
##### ---------- ##### ---------- ##### ---------- ##### ---------- #####
#!/usr/bin/env python3
"""
Script to run optimizations for the models
"""
import dadi
import numpy as np
import matplotlib.pyplot as plt
import os
output_dir = 'resize_models_results'
os.makedirs(output_dir, exist_ok=True)
data_fs = "[/path/to/dadi-formatted/SFS-file-name].fs"
fs = dadi.Spectrum.from_file(data_fs)
if not fs.folded:
    fs = fs.fold()
pts_l = [40, 50, 60] #grid points for extrapolation

# define the models to test
# Model 0: No size change (constant population size)
def no_change(params, ns, pts):
    """
    No size change model.
    params = [nu]
    ns = sample size
    """
    nu = params[0]
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, 5.0, nu)
    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs
    
# Model 1: One size change
def one_change(params, ns, pts):
    """
    One size change model.
    params = [nu, T]
    ns = sample size
    """
    nu, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T, nu)
    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs

# Model 2: Two size changes
def two_changes(params, ns, pts):
    """
    Two size changes model.
    params = [nu1, nu2, T1, T2]
    ns = sample size
    """
    nu1, nu2, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)
    phi = dadi.Integration.one_pop(phi, xx, T2, nu2)
    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs

# Model 3: Three size changes
def three_changes(params, ns, pts):
    """
    Three size changes model.
    params = [nu1, nu2, nu3, T1, T2, T3]
    ns = sample size
    """
    nu1, nu2, nu3, T1, T2, T3 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)
    phi = dadi.Integration.one_pop(phi, xx, T2, nu2)
    phi = dadi.Integration.one_pop(phi, xx, T3, nu3)
    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs

# Model 4: Four size changes
def four_changes(params, ns, pts):
    """
    Four size changes model.
    params = [nu1, nu2, nu3, nu4, T1, T2, T3, T4]
    ns = sample size
    """
    nu1, nu2, nu3, nu4, T1, T2, T3, T4 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)
    phi = dadi.Integration.one_pop(phi, xx, T2, nu2)
    phi = dadi.Integration.one_pop(phi, xx, T3, nu3)
    phi = dadi.Integration.one_pop(phi, xx, T4, nu4)
    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs

# model configurations
models = [
    {
        'name': 'no_change',
        'func': no_change,
        'param_names': ['nu'],
        'lower_bound': [0.01],
        'upper_bound': [100],
        'initial_values': [1.0]
    },
    {
        'name': 'one_change',
        'func': one_change,
        'param_names': ['nu', 'T'],
        'lower_bound': [0.01, 0.01],
        'upper_bound': [100, 10],
        'initial_values': [2.0, 0.5]
    },
    {
        'name': 'two_changes',
        'func': two_changes,
        'param_names': ['nu1', 'nu2', 'T1', 'T2'],
        'lower_bound': [0.01, 0.01, 0.01, 0.01],
        'upper_bound': [100, 100, 10, 10],
        'initial_values': [2.0, 0.5, 0.5, 0.1]
    },
    {
        'name': 'three_changes',
        'func': three_changes,
        'param_names': ['nu1', 'nu2', 'nu3', 'T1', 'T2', 'T3'],
        'lower_bound': [0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        'upper_bound': [100, 100, 100, 10, 10, 10],
        'initial_values': [2.0, 1.0, 0.5, 0.5, 0.2, 0.1]
    },
    {
        'name': 'four_changes',
        'func': four_changes,
        'param_names': ['nu1', 'nu2', 'nu3', 'nu4', 'T1', 'T2', 'T3', 'T4'],
        'lower_bound': [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        'upper_bound': [100, 100, 100, 100, 10, 10, 10, 10],
        'initial_values': [2.0, 1.5, 1.0, 0.5, 0.5, 0.3, 0.2, 0.1]
    }
]

# optimization for each model (five)
results = []
for model in models:
    print(f"\nOptimizing {model['name']} model...")
    func_ex = dadi.Numerics.make_extrap_log_func(model['func'])
    best_p = None
    best_ll = -np.inf
    for i in range(5):
        print(f"  Optimization run {i+1}/5")
        if i == 0:
            p0 = model['initial_values']
        else:
            p0 = [v * np.random.uniform(0.8, 1.2) for v in model['initial_values']]
            p0 = [max(min(p0[j], model['upper_bound'][j]), model['lower_bound'][j]) 
                  for j in range(len(p0))]
        popt = dadi.Inference.optimize_log(
            p0, fs, func_ex, pts_l,
            lower_bound=model['lower_bound'],
            upper_bound=model['upper_bound'],
            verbose=False, maxiter=100
        )
        model_fs = func_ex(popt, fs.sample_sizes, pts_l)
        ll = dadi.Inference.ll_multinom(model_fs, fs)
        print(f"    log-likelihood: {ll}")
        if ll > best_ll:
            best_ll = ll
            best_p = popt
    k = len(model['param_names'])
    aic = 2 * k - 2 * best_ll
    best_model = func_ex(best_p, fs.sample_sizes, pts_l)
    results.append({
        'name': model['name'],
        'parameters': dict(zip(model['param_names'], best_p)),
        'log_likelihood': best_ll,
        'AIC': aic,
        'model_fs': best_model
    })
    with open(f"{output_dir}/{model['name']}_results.txt", 'w') as f:
        f.write(f"Model: {model['name']}\n")
        f.write(f"Parameters: {dict(zip(model['param_names'], best_p))}\n")
        f.write(f"Log-likelihood: {best_ll}\n")
        f.write(f"AIC: {aic}\n")
# model vs data plot 
    fig = plt.figure(figsize=(8, 6))
    dadi.Plotting.plot_1d_comp_multinom(best_model, fs)
    plt.title(f"{model['name']} (AIC = {aic:.2f})")
    plt.savefig(f"{output_dir}/{model['name']}_fit.png", dpi=300)
    plt.close()

# Compare models by AIC
results.sort(key=lambda x: x['AIC'])
best_model = results[0]
print("\n\nModel Comparison by AIC:")
for result in results:
    print(f"{result['name']}: AIC = {result['AIC']:.2f}, Δ AIC = {result['AIC'] - best_model['AIC']:.2f}")
print(f"\nBest model: {best_model['name']} (AIC = {best_model['AIC']:.2f})")
print("Parameters:", best_model['parameters'])
print("Log-likelihood:", best_model['log_likelihood'])

# comparison plot of all models
plt.figure(figsize=(10, 8))
for i, result in enumerate(results):
    plt.subplot(len(results), 1, i+1)
    dadi.Plotting.plot_1d_comp_multinom(result['model_fs'], fs, residual=True)
    plt.title(f"{result['name']} (AIC = {result['AIC']:.2f})")

plt.tight_layout()
plt.savefig(f"{output_dir}/all_models_comparison.png", dpi=300)
plt.close()
# summary table
with open(f"{output_dir}/model_comparison.txt", 'w') as f:
    f.write("Model Comparison by AIC:\n")
    f.write("-------------------------\n")
    for result in results:
        f.write(f"{result['name']}: AIC = {result['AIC']:.2f}, Δ AIC = {result['AIC'] - best_model['AIC']:.2f}\n")
    f.write(f"\nBest model: {best_model['name']} (AIC = {best_model['AIC']:.2f})\n")
    f.write(f"Parameters: {best_model['parameters']}\n")
    f.write(f"Log-likelihood: {best_model['log_likelihood']}\n")
