#Part4. PCC analysis
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from drift_simulation_3600s import PARAMS, drift_rhs, IDX
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10
rcParams['mathtext.default'] = 'regular'
os.makedirs('figure1', exist_ok=True)
 
def perturb_params(base, rng, noise=0.2):
    p = {k: v.copy() if isinstance(v, dict) else v for (k, v) in base.items()}
 
    def g(mu):
        return mu * np.clip(rng.normal(1.0, noise), 0.4, 1.6)
    p['gamma_0'] = g(p['gamma_0'])
    p['theta_0'] = g(p['theta_0'])
    p['theta_eq'] = g(p['theta_eq'])
    p['k_wetting'] = g(p['k_wetting'])
    p['R'] = g(p['R'])
    p['alpha_SNARE'] = g(p['alpha_SNARE'])
    p['beta_ESCRT'] = g(p['beta_ESCRT'])
    p['k_fus'] = g(p['k_fus'])
    p['k_endo'] = g(p['k_endo'])
    p['k_traffic'] = g(p['k_traffic'])
    p['k_lyso'] = g(p['k_lyso'])
    p['E_c'] = g(p['E_c'])
    for k in ('fus', 'rep', 'rec', 'endo', 'osm'):
        p['k_act'][k] = g(p['k_act'][k])
        p['k_dec'][k] = g(p['k_dec'][k])
    return p
 
def run_one(p, t_end=3600.0, n_pts=120):
    y0 = [p['theta_0'], p['gamma_0'], 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, p['GP_0'], 1.0, 0.0, 0.0, 0.0, 0.0]
    t_eval = np.linspace(0.0, t_end, n_pts)
    sol = solve_ivp(drift_rhs, [0.0, t_end], y0, t_eval=t_eval, args=(p,), method='RK45', rtol=1e-06, atol=1e-09)
    if not sol.success:
        return None
    last = sol.y[:, -1]
    N = last[[IDX['N_PM'], IDX['N_EE'], IDX['N_LE'], IDX['N_Lyso'], IDX['N_Cyto']]]
    M = last[[IDX['M_fus'], IDX['M_rep'], IDX['M_rec'], IDX['M_endo'], IDX['M_osm']]]
    M = M / M.sum()
    return np.concatenate([N, M])
 
def monte_carlo(label, n=100, seed=42):
    base = PARAMS[label]
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        p = perturb_params(base, rng)
        v = run_one(p)
        if v is not None and np.all(np.isfinite(v)):
            rows.append(v)
    arr = np.vstack(rows)
    print(f'  [{label}] {len(rows)}/{n} simulations succeeded')
    return arr
 
def pearson_matrix(arr):
    R = np.zeros((5, 5))
    P = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            (x, y) = (arr[:, i], arr[:, 5 + j])
            if np.std(x) < 1e-12 or np.std(y) < 1e-12:
                (R[i, j], P[i, j]) = (0.0, 1.0)
            else:
                (r, pv) = pearsonr(x, y)
                R[i, j] = r
                P[i, j] = pv
    return (R, P)
 
def plot_heatmap(R, P, title, save_path):
    (fig, ax) = plt.subplots(figsize=(7, 6))
    im = ax.imshow(R, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label('Pearson r', rotation=270, labelpad=15)
    dist_lbl = ['N_PM', 'N_EE', 'N_LE', 'N_Lyso', 'N_Cyto']
    mod_lbl = ['fusion', 'repair', 'receptor', 'endocyt.', 'osmotic']
    ax.set_xticks(range(5))
    ax.set_xticklabels(mod_lbl, rotation=30, ha='right')
    ax.set_yticks(range(5))
    ax.set_yticklabels(dist_lbl)
    ax.set_xlabel('Module weights ($w_k$)')
    ax.set_ylabel('Subcellular distribution')
    ax.set_title(title, fontsize=12, fontweight='bold')
    for i in range(5):
        for j in range(5):
            r = R[i, j]
            pv = P[i, j]
            star = '***' if pv < 0.001 else '**' if pv < 0.01 else '*' if pv < 0.05 else ''
            ax.text(j, i, f'{r:+.2f}\n{star}', ha='center', va='center', color='white' if abs(r) > 0.6 else 'black', fontsize=9, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  [OK] saved: {save_path}')
 
def save_xlsx(R, P, label):
    dist_lbl = ['N_PM', 'N_EE', 'N_LE', 'N_Lyso', 'N_Cyto']
    mod_lbl = ['w_fusion', 'w_repair', 'w_receptor', 'w_endocytosis', 'w_osmotic']
    df_r = pd.DataFrame(R, index=dist_lbl, columns=mod_lbl)
    df_p = pd.DataFrame(P, index=dist_lbl, columns=mod_lbl)
    fname = f'pearson_correlation_{label}.xlsx'
    with pd.ExcelWriter(fname) as wb:
        df_r.to_excel(wb, sheet_name='Pearson_r')
        df_p.to_excel(wb, sheet_name='p_value')
    print(f'  [OK] saved: {fname}')
 
def main():
    print('=' * 80)
    print('Fig. 4 j-k - Pearson correlation analysis (Monte Carlo n=100)')
    print('=' * 80)
    for label in ('micro', 'nano'):
        print(f'\n>>> {label}-BMC')
        arr = monte_carlo(label, n=100, seed=42)
        (R, P) = pearson_matrix(arr)
        print(f'  Pearson r matrix [5 distribution rows x 5 module cols]:')
        dist_lbl = ['N_PM    ', 'N_EE    ', 'N_LE    ', 'N_Lyso  ', 'N_Cyto  ']
        mod_lbl = ['fusion', 'repair', 'receptor', 'endocyt.', 'osmotic']
        print('              ' + '  '.join([f'{m:>9s}' for m in mod_lbl]))
        for (i, dl) in enumerate(dist_lbl):
            row = '    ' + dl + '  ' + '  '.join([f'{R[i, j]:+9.3f}' for j in range(5)])
            print(row)
        title = f"Fig. 4{('j' if label == 'micro' else 'k')} | "
        title += f"{('Micro' if label == 'micro' else 'Nano')}-BMC Pearson correlation"
        title += ' (n = 100 MC)'
        save = f"figure1/Fig4{('j' if label == 'micro' else 'k')}_{label}_pearson.png"
        plot_heatmap(R, P, title, save)
        save_xlsx(R, P, label)
    print('\nDone.')
if __name__ == '__main__':
    main()
