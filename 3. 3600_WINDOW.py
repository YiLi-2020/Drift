#Part3. 3600s simulation window
import os
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10
rcParams['mathtext.default'] = 'regular'
os.makedirs('figure1', exist_ok=True)
PARAMS = {'micro': {'gamma_0': 0.025, 'theta_0': 120.0, 'theta_eq': 88.7, 'k_wetting': 0.05, 'R': 1.5e-06, 'k_act': {'fus': 0.06, 'rep': 0.05, 'rec': 0.015, 'endo': 0.012, 'osm': 0.04}, 'k_dec': {'fus': 0.012, 'rep': 0.012, 'rec': 0.025, 'endo': 0.025, 'osm': 0.015}, 'tau': {'fus': 10.0, 'rep': 25.0, 'rec': 120.0, 'endo': 120.0, 'osm': 60.0}, 'sigma': {'fus': 5.0, 'rep': 8.0, 'rec': 30.0, 'endo': 30.0, 'osm': 20.0}, 'alpha_SNARE': 0.5, 'beta_ESCRT': 0.3, 'k_gamma': 0.005, 'k_dmg': 0.025, 'E_c': 8e-14, 'k_disorder': 0.008, 'k_order': 0.0015, 'GP_0': -0.28, 'GP_inf': 0.2, 'k_fus': 0.02, 'k_endo': 0.0008, 'k_traffic': 0.005, 'k_lyso': 0.003}, 'nano': {'gamma_0': 0.015, 'theta_0': 60.0, 'theta_eq': 44.3, 'k_wetting': 0.1, 'R': 1e-07, 'k_act': {'fus': 0.025, 'rep': 0.005, 'rec': 0.05, 'endo': 0.06, 'osm': 0.005}, 'k_dec': {'fus': 0.02, 'rep': 0.025, 'rec': 0.012, 'endo': 0.012, 'osm': 0.025}, 'tau': {'fus': 30.0, 'rep': 120.0, 'rec': 10.0, 'endo': 10.0, 'osm': 120.0}, 'sigma': {'fus': 15.0, 'rep': 30.0, 'rec': 5.0, 'endo': 5.0, 'osm': 30.0}, 'alpha_SNARE': 0.3, 'beta_ESCRT': 0.2, 'k_gamma': 0.005, 'k_dmg': 0.025, 'E_c': 3e-15, 'k_disorder': 0.004, 'k_order': 0.0015, 'GP_0': -0.28, 'GP_inf': 0.17, 'k_fus': 0.003, 'k_endo': 0.012, 'k_traffic': 0.005, 'k_lyso': 0.003}}
IDX = {'theta': 0, 'gamma_eff': 1, 'M_fus': 2, 'M_rep': 3, 'M_rec': 4, 'M_endo': 5, 'M_osm': 6, 'I': 7, 'GP': 8, 'N_PM': 9, 'N_EE': 10, 'N_LE': 11, 'N_Lyso': 12, 'N_Cyto': 13}
 
def drift_rhs(t, y, p):
    (theta, g_eff) = (y[0], y[1])
    (Mfus, Mrep, Mrec, Mendo, Mosm) = (y[2], y[3], y[4], y[5], y[6])
    (I, GP) = (y[7], y[8])
    (N_PM, N_EE, N_LE, N_Lyso, N_Cyto) = (y[9], y[10], y[11], y[12], y[13])
    dtheta = -p['k_wetting'] * (theta - p['theta_eq'])
    A_contact = np.pi * p['R'] ** 2 * np.sin(np.deg2rad(theta)) ** 2
    E = g_eff * A_contact
    dg_eff = -p['alpha_SNARE'] * Mfus * g_eff + p['k_gamma'] * (p['gamma_0'] - g_eff)
 
    def S(k):
        return 1.0 / (1.0 + np.exp(-(t - p['tau'][k]) / p['sigma'][k]))
    dMfus = p['k_act']['fus'] * S('fus') * (1 - Mfus) - p['k_dec']['fus'] * Mfus
    dMrep = p['k_act']['rep'] * S('rep') * (1 - Mrep) - p['k_dec']['rep'] * Mrep
    dMrec = p['k_act']['rec'] * S('rec') * (1 - Mrec) - p['k_dec']['rec'] * Mrec
    dMendo = p['k_act']['endo'] * S('endo') * (1 - Mendo) - p['k_dec']['endo'] * Mendo
    dMosm = p['k_act']['osm'] * S('osm') * (1 - Mosm) - p['k_dec']['osm'] * Mosm
    dI = -p['k_dmg'] * (E / p['E_c']) * I + p['beta_ESCRT'] * Mrep * (1 - I)
    dGP = -p['k_disorder'] * Mfus * max(0.0, 1 - I) + p['k_order'] * (p['GP_inf'] - GP)
    flux_fus = p['k_fus'] * Mfus * N_PM
    flux_endo = p['k_endo'] * Mendo * N_PM
    dN_PM = -(flux_fus + flux_endo)
    dN_EE = flux_endo - p['k_traffic'] * N_EE
    dN_LE = p['k_traffic'] * N_EE - p['k_lyso'] * N_LE
    dN_Lyso = p['k_lyso'] * N_LE
    dN_Cyto = flux_fus
    return [dtheta, dg_eff, dMfus, dMrep, dMrec, dMendo, dMosm, dI, dGP, dN_PM, dN_EE, dN_LE, dN_Lyso, dN_Cyto]
 
def simulate(label: str, t_end: float=3600.0, n_pts: int=720) -> pd.DataFrame:
    p = PARAMS[label]
    y0 = [p['theta_0'], p['gamma_0'], 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, p['GP_0'], 1.0, 0.0, 0.0, 0.0, 0.0]
    t_eval = np.linspace(0.0, t_end, n_pts)
    sol = solve_ivp(drift_rhs, [0.0, t_end], y0, t_eval=t_eval, args=(p,), method='RK45', rtol=1e-07, atol=1e-10)
    df = pd.DataFrame({'time_s': sol.t, 'Contact_angle_deg': sol.y[IDX['theta']], 'Surface_tension_N_m': sol.y[IDX['gamma_eff']], 'M_fusion': sol.y[IDX['M_fus']], 'M_repair': sol.y[IDX['M_rep']], 'M_receptor': sol.y[IDX['M_rec']], 'M_endocytosis': sol.y[IDX['M_endo']], 'M_osmotic': sol.y[IDX['M_osm']], 'Membrane_integrity': sol.y[IDX['I']], 'GP_value': sol.y[IDX['GP']], 'BMC_plasma_membrane': sol.y[IDX['N_PM']], 'BMC_early_endosome': sol.y[IDX['N_EE']], 'BMC_late_endosome': sol.y[IDX['N_LE']], 'BMC_lysosome': sol.y[IDX['N_Lyso']], 'BMC_cytoplasm': sol.y[IDX['N_Cyto']]})
    df['Contact_area_m2'] = np.pi * p['R'] ** 2 * np.sin(np.deg2rad(df['Contact_angle_deg'])) ** 2
    df['Interface_energy_J'] = df['Surface_tension_N_m'] * df['Contact_area_m2']
    return df
 
def main():
    print('=' * 80)
    print('DRIFT 3600s endpoint simulation')
    print('=' * 80)
    df_micro = simulate('micro', t_end=3600.0, n_pts=720)
    df_nano = simulate('nano', t_end=3600.0, n_pts=720)
    df_micro.to_excel('pathway_simulation_micro_3600s.xlsx', index=False)
    df_nano.to_excel('pathway_simulation_nano_3600s.xlsx', index=False)
    print('[OK] saved: pathway_simulation_micro_3600s.xlsx')
    print('[OK] saved: pathway_simulation_nano_3600s.xlsx')
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        last = df.iloc[-1]
        print(f'\n--- {lbl} @ t=3600s ---')
        print(f"  theta            = {last['Contact_angle_deg']:.2f} deg")
        print(f"  E_peak           = {df['Interface_energy_J'].max():.2e} J")
        print(f"  GP (final)       = {last['GP_value']:+.3f}")
        print(f"  I (final)        = {last['Membrane_integrity']:.3f}")
        print(f'  Module weights (steady-state, normalized):')
        cols = ['M_fusion', 'M_repair', 'M_receptor', 'M_endocytosis', 'M_osmotic']
        wend = df.iloc[-1][cols].values
        wnorm = wend / wend.sum()
        for (c, w) in zip(cols, wnorm):
            print(f'    {c:18s} = {w:.3f}')
        print(f'  Subcellular distribution:')
        dist_cols = ['BMC_plasma_membrane', 'BMC_early_endosome', 'BMC_late_endosome', 'BMC_lysosome', 'BMC_cytoplasm']
        for c in dist_cols:
            print(f'    {c:24s} = {last[c]:.3f}')
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 4, hspace=0.5, wspace=0.4)
    colors = {'micro': '#E63946', 'nano': '#2A9D8F'}
    ax = fig.add_subplot(gs[0, 0])
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.plot(df['time_s'], df['Contact_angle_deg'], color=colors[lbl], lw=2, label=lbl)
    ax.set_xlim(0, 300)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$\\theta$ (deg)')
    ax.set_title('(a) Contact angle  [0-300 s]')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[0, 1])
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.semilogy(df['time_s'], df['Interface_energy_J'].clip(lower=1e-25), color=colors[lbl], lw=2, label=lbl)
    ax.set_xlim(0, 300)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$E$ (J)')
    ax.set_title('(b) Interfacial energy  [0-300 s]')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[0, 2])
    for (c, color) in zip(['M_fusion', 'M_repair', 'M_receptor', 'M_endocytosis', 'M_osmotic'], ['#E63946', '#F4A261', '#2A9D8F', '#264653', '#9D4EDD']):
        ax.plot(df_micro['time_s'], df_micro[c], color=color, lw=1.8, label=c.replace('M_', ''))
    ax.set_xlim(0, 300)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$M_k(t)$')
    ax.set_title('(c) Modules - micro [0-300 s]')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[0, 3])
    for (c, color) in zip(['M_fusion', 'M_repair', 'M_receptor', 'M_endocytosis', 'M_osmotic'], ['#E63946', '#F4A261', '#2A9D8F', '#264653', '#9D4EDD']):
        ax.plot(df_nano['time_s'], df_nano[c], color=color, lw=1.8, label=c.replace('M_', ''))
    ax.set_xlim(0, 300)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$M_k(t)$')
    ax.set_title('(d) Modules - nano  [0-300 s]')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 0])
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.plot(df['time_s'] / 60.0, df['Membrane_integrity'], color=colors[lbl], lw=2, label=lbl)
    ax.set_xlim(0, 60)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('$I(t)$')
    ax.set_title('(e) Integrity  [0-60 min]')
    ax.legend()
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 1])
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        ax.plot(df['time_s'] / 60.0, df['GP_value'], color=colors[lbl], lw=2, label=lbl)
    ax.scatter([60, 60], [0.19, 0.17], color=['#E63946', '#2A9D8F'], s=80, zorder=5, edgecolor='black', label='1 h C-Laurdan exp.')
    ax.set_xlim(0, 60)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('GP value')
    ax.set_title('(f) GP dynamics  [0-60 min]')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 2])
    cols = ['BMC_plasma_membrane', 'BMC_early_endosome', 'BMC_late_endosome', 'BMC_lysosome', 'BMC_cytoplasm']
    palette = ['#A5A5A5', '#F4A261', '#E76F51', '#9D4EDD', '#2A9D8F']
    ax.stackplot(df_micro['time_s'] / 60.0, *[df_micro[c] for c in cols], labels=[c.replace('BMC_', '') for c in cols], colors=palette, alpha=0.85)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Fraction')
    ax.set_title('(g) Compartments - micro')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 3])
    ax.stackplot(df_nano['time_s'] / 60.0, *[df_nano[c] for c in cols], labels=[c.replace('BMC_', '') for c in cols], colors=palette, alpha=0.85)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Fraction')
    ax.set_title('(h) Compartments - nano')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[2, 0])
    for (c, color) in zip(['M_fusion', 'M_repair', 'M_receptor', 'M_endocytosis', 'M_osmotic'], ['#E63946', '#F4A261', '#2A9D8F', '#264653', '#9D4EDD']):
        ax.plot(df_micro['time_s'] / 60.0, df_micro[c], color=color, lw=1.8, label=c.replace('M_', ''))
    ax.set_xlim(0, 60)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('$M_k(t)$')
    ax.set_title('(i) Modules - micro [0-60 min]')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[2, 1])
    for (c, color) in zip(['M_fusion', 'M_repair', 'M_receptor', 'M_endocytosis', 'M_osmotic'], ['#E63946', '#F4A261', '#2A9D8F', '#264653', '#9D4EDD']):
        ax.plot(df_nano['time_s'] / 60.0, df_nano[c], color=color, lw=1.8, label=c.replace('M_', ''))
    ax.set_xlim(0, 60)
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('$M_k(t)$')
    ax.set_title('(j) Modules - nano [0-60 min]')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[2, 2:], polar=True)
    cols = ['M_fusion', 'M_repair', 'M_receptor', 'M_endocytosis', 'M_osmotic']
    labels = ['fusion', 'repair', 'receptor', 'endocyt.', 'osmotic']
    ang = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    ang += ang[:1]
    for (lbl, df) in [('micro', df_micro), ('nano', df_nano)]:
        v = df.iloc[-1][cols].values
        v = v / v.sum()
        v = v.tolist() + [v[0]]
        ax.plot(ang, v, 'o-', color=colors[lbl], lw=2, label=lbl)
        ax.fill(ang, v, color=colors[lbl], alpha=0.2)
    ax.set_xticks(ang[:-1])
    ax.set_xticklabels(labels)
    ax.set_title('(k) Steady-state module weights')
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3)
    plt.suptitle('Fig. 4 (revised): DRIFT dual-time-window simulation (0-300 s early-phase + 0-60 min endpoint)', fontsize=14, fontweight='bold', y=1.02)
    plt.savefig('figure1/Fig4_DualTimeWindow.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print('\n[OK] saved: figure1/Fig4_DualTimeWindow.png')
    print('\nDone.')
if __name__ == '__main__':
    main()
