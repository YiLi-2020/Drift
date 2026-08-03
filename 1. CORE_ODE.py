#Part 1_4 module core ODE model (300s)
import os
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
os.makedirs('outputs', exist_ok=True)
MODULE_PARAMS = {'snare': {'k_act': 0.1157, 'k_dec': 0.009644, 'tau': 19.43, 'sig': 10.0}, 'escrt': {'k_act': 0.1025, 'k_dec': 0.0114, 'tau': 13.11, 'sig': 4.49}, 'glyco': {'k_act': 0.1426, 'k_dec': 0.01426, 'tau': 14.3, 'sig': 10.01}, 'pik3c': {'k_act': 0.07804, 'k_dec': 0.01171, 'tau': 24.56, 'sig': 10.0}}
SIZE_PARAMS = {'micro': {'gamma0': 0.025, 'theta0': 120.0, 'R': 1.5e-06, 'k_wet': 0.098, 'k_dmg': 0.445335254, 'E_c': 1.56614817e-12, 'beta_escrt': 0.000322578294, 'k_disorder': 0.00309425746, 'k_order': 0.0175571822, 'GP0': -0.28, 'GPinf': -0.269563632, 'k_fus': 0.0848689187, 'k_endo': 0.00627003059, 'k_traf': 0.0223608667, 'k_lyso': 0.0107779832}, 'nano': {'gamma0': 0.015, 'theta0': 60.0, 'R': 1e-07, 'k_wet': 0.098, 'k_dmg': 0.00117725953, 'E_c': 4.97508682e-19, 'beta_escrt': 0.394091784, 'k_disorder': 0.00713892771, 'k_order': 0.0115884103, 'GP0': -0.28, 'GPinf': -0.116451303, 'k_fus': 0.016376311, 'k_endo': 0.0792104649, 'k_traf': 0.0273488643, 'k_lyso': 0.0200161786}}
ALPHA_SNARE = 0.5
BETA_GLYCO = 0.4
 
def _sigmoid(t, tau, sig):
    return 1.0 / (1.0 + np.exp(-(t - tau) / sig))
 
def bmc_rhs(t, y, p, mod):
    (snare, escrt, glyco, pik3c, theta, I, GP, Nmem, Nee, Nle, Nly, Ncy) = y
    dsnare = mod['snare']['k_act'] * _sigmoid(t, mod['snare']['tau'], mod['snare']['sig']) * (1 - snare) - mod['snare']['k_dec'] * snare
    descrt = mod['escrt']['k_act'] * _sigmoid(t, mod['escrt']['tau'], mod['escrt']['sig']) * (1 - escrt) - mod['escrt']['k_dec'] * escrt
    dglyco = mod['glyco']['k_act'] * _sigmoid(t, mod['glyco']['tau'], mod['glyco']['sig']) * (1 - glyco) - mod['glyco']['k_dec'] * glyco
    dpik3c = mod['pik3c']['k_act'] * _sigmoid(t, mod['pik3c']['tau'], mod['pik3c']['sig']) * (1 - pik3c) - mod['pik3c']['k_dec'] * pik3c
    gamma = p['gamma0'] * (1.0 - ALPHA_SNARE * snare)
    theta_eq = p['theta0'] * (1.0 - BETA_GLYCO * glyco)
    dtheta = -p['k_wet'] * (theta - theta_eq)
    area = np.pi * p['R'] ** 2 * np.sin(np.deg2rad(np.clip(theta, 0.0, 180.0))) ** 2
    energy = max(gamma, 0.0) * area
    dI = -p['k_dmg'] * (energy / p['E_c']) * I + p['beta_escrt'] * escrt * (1.0 - I)
    dGP = -p['k_disorder'] * snare * max(0.0, 1.0 - I) + p['k_order'] * (p['GPinf'] - GP)
    fus = p['k_fus'] * snare * (0.15 + 0.85 * max(0.0, 1.0 - I)) * Nmem
    endo = p['k_endo'] * glyco * I / (1.0 + pik3c) * Nmem
    dNmem = -(fus + endo)
    dNee = endo - p['k_traf'] * Nee
    dNle = p['k_traf'] * Nee - p['k_lyso'] * Nle
    dNly = p['k_lyso'] * Nle
    dNcy = fus
    return [dsnare, descrt, dglyco, dpik3c, dtheta, dI, dGP, dNmem, dNee, dNle, dNly, dNcy]
 
def simulate(label, t_end=299.5, n_pts=600, size_overrides=None, module_overrides=None):
    p = dict(SIZE_PARAMS[label])
    if size_overrides:
        p.update(size_overrides)
    mod = {k: dict(v) for (k, v) in MODULE_PARAMS.items()}
    if module_overrides:
        for (k, v) in module_overrides.items():
            mod[k].update(v)
    y0 = [0.1, 0.1, 0.1, 0.1, p['theta0'], 1.0, p['GP0'], 1.0, 0.0, 0.0, 0.0, 0.0]
    t_eval = np.linspace(0.0, t_end, n_pts)
    sol = solve_ivp(bmc_rhs, [0.0, t_end], y0, t_eval=t_eval, args=(p, mod), method='RK45', rtol=1e-07, atol=1e-10)
    (snare, escrt, glyco, pik3c, theta, I, GP, Nmem, Nee, Nle, Nly, Ncy) = sol.y
    gamma = p['gamma0'] * (1.0 - ALPHA_SNARE * snare)
    area = np.pi * p['R'] ** 2 * np.sin(np.deg2rad(np.clip(theta, 0.0, 180.0))) ** 2
    energy = gamma * area
    theta_eq_final = p['theta0'] * (1.0 - BETA_GLYCO * glyco[-1])
    wetting = np.clip((p['theta0'] - theta) / max(p['theta0'] - theta_eq_final, 1e-12), 0.0, 1.0)
    dmg = np.maximum(-np.gradient(I, sol.t), 0.0)
    repair = np.maximum(np.gradient(I, sol.t), 0.0)
    df = pd.DataFrame({'time_s': sol.t, 'SNARE_activity': snare, 'ESCRT_activity': escrt, 'Glycoprotein_activity': glyco, 'PIK3C_activity': pik3c, 'Wetting_degree': wetting, 'Contact_angle_deg': theta, 'Contact_area_m2': area, 'Surface_tension_N_m': gamma, 'Interface_energy_J': energy, 'Membrane_integrity': I, 'GP_value': GP, 'BMC_membrane': Nmem, 'BMC_early_endosome': Nee, 'BMC_late_endosome': Nle, 'BMC_lysosome': Nly, 'BMC_cytoplasm': Ncy, 'Cumulative_damage': np.cumsum(dmg) * (sol.t[1] - sol.t[0] if len(sol.t) > 1 else 1.0), 'Cumulative_repair': np.cumsum(repair) * (sol.t[1] - sol.t[0] if len(sol.t) > 1 else 1.0), 'Cumulative_wetting': wetting, 'Cumulative_entry': 1.0 - Nmem})
    return df
 
def plot_core_figures(df_micro, df_nano):
    colors = {'micro': '#E63946', 'nano': '#2A9D8F'}
    (fig, axes) = plt.subplots(2, 3, figsize=(14, 8))
    ax = axes[0, 0]
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.plot(df['time_s'], df['Wetting_degree'], color=colors[lbl], lw=2, label=lbl)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Wetting')
    ax.set_title('Wetting kinetics')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[0, 1]
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.plot(df['time_s'], df['Contact_angle_deg'], color=colors[lbl], lw=2, label=lbl)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Contact angle (deg)')
    ax.set_title('Contact angle')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[0, 2]
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.semilogy(df['time_s'], df['Interface_energy_J'].clip(lower=1e-25), color=colors[lbl], lw=2, label=lbl)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Interfacial energy (J)')
    ax.set_title('Interface energy')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[1, 0]
    for (col, c) in zip(['SNARE_activity', 'ESCRT_activity', 'Glycoprotein_activity', 'PIK3C_activity'], ['#E63946', '#F4A261', '#2A9D8F', '#264653']):
        ax.plot(df_micro['time_s'], df_micro[col], color=c, lw=1.8, label=col.replace('_activity', ''))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Activity')
    ax.set_title('Pathway activity (micro)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1, 1]
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.plot(df['time_s'], df['Membrane_integrity'], color=colors[lbl], lw=2, label=lbl)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Integrity')
    ax.set_title('Membrane integrity')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = axes[1, 2]
    cols = ['BMC_membrane', 'BMC_early_endosome', 'BMC_late_endosome', 'BMC_lysosome', 'BMC_cytoplasm']
    ax.stackplot(df_micro['time_s'], *[df_micro[c] for c in cols], labels=[c.replace('BMC_', '') for c in cols], alpha=0.85)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Fraction')
    ax.set_title('BMC distribution (micro)')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig1_Droplet_Physics_Dynamics.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    (fig, axes) = plt.subplots(1, 2, figsize=(12, 4.5))
    for (ax, df, title) in [(axes[0], df_micro, 'micro'), (axes[1], df_nano, 'nano')]:
        for (col, c) in zip(['SNARE_activity', 'ESCRT_activity', 'Glycoprotein_activity', 'PIK3C_activity'], ['#E63946', '#F4A261', '#2A9D8F', '#264653']):
            ax.plot(df['time_s'], df[col], color=c, lw=2, label=col.replace('_activity', ''))
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Activity')
        ax.set_title(f'Four pathways - {title}')
        ax.legend()
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig2_Four_Pathway_Comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    (fig, axes) = plt.subplots(1, 2, figsize=(12, 4.5))
    cols = ['BMC_membrane', 'BMC_early_endosome', 'BMC_late_endosome', 'BMC_lysosome', 'BMC_cytoplasm']
    palette = ['#A5A5A5', '#F4A261', '#E76F51', '#9D4EDD', '#2A9D8F']
    for (ax, df, title) in [(axes[0], df_micro, 'micro'), (axes[1], df_nano, 'nano')]:
        ax.stackplot(df['time_s'], *[df[c] for c in cols], labels=[c.replace('BMC_', '') for c in cols], colors=palette, alpha=0.85)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Fraction')
        ax.set_title(f'BMC distribution - {title}')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig3_BMC_Distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    (fig, axes) = plt.subplots(1, 2, figsize=(12, 4.5))
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        axes[0].plot(df['time_s'], df['GP_value'], color=colors[lbl], lw=2, label=lbl)
        axes[1].plot(df['time_s'], df['Membrane_integrity'], color=colors[lbl], lw=2, label=lbl)
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('GP')
    axes[0].set_title('GP dynamics')
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Integrity')
    axes[1].set_title('Membrane integrity')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig7_Membrane_Effects.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
 
def main():
    df_micro = simulate('micro')
    df_nano = simulate('nano')
    df_micro.to_excel('outputs/simulation_micro_BMC.xlsx', index=False)
    df_nano.to_excel('outputs/simulation_nano_BMC.xlsx', index=False)
    plot_core_figures(df_micro, df_nano)
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        last = df.iloc[-1]
        print(lbl, 'cyto', float(last['BMC_cytoplasm']), 'I', float(last['Membrane_integrity']), 'theta', float(last['Contact_angle_deg']))
if __name__ == '__main__':
    main()
