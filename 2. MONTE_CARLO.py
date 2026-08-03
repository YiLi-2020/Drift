#Part II_Monte Carlo
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bmc_core_model import MODULE_PARAMS, SIZE_PARAMS, simulate
os.makedirs('outputs', exist_ok=True)
 
def run_sensitivity():
    rows = []
    for value in [10, 15, 20, 25, 30]:
        df = simulate('micro', size_overrides={'gamma0': value / 1000.0})
        last = df.iloc[-1]
        rows.append({'Parameter': 'Surface_tension', 'Value': value, 'Unit': 'mN/m', 'Final_uptake': float(last['BMC_cytoplasm']), 'Final_integrity': float(last['Membrane_integrity'])})
    for value in [30, 60, 90, 120, 150]:
        df = simulate('micro', size_overrides={'theta0': float(value)})
        last = df.iloc[-1]
        rows.append({'Parameter': 'Contact_angle', 'Value': value, 'Unit': 'degree', 'Final_uptake': float(last['BMC_cytoplasm']), 'Final_integrity': float(last['Membrane_integrity'])})
    for value in [0.02, 0.05, 0.08, 0.1, 0.15]:
        df = simulate('micro', size_overrides={'k_wet': float(value)})
        last = df.iloc[-1]
        rows.append({'Parameter': 'Wetting_rate', 'Value': value, 'Unit': 's^-1', 'Final_uptake': float(last['BMC_cytoplasm']), 'Final_integrity': float(last['Membrane_integrity'])})
    out = pd.DataFrame(rows)
    out.to_excel('outputs/sensitivity_analysis.xlsx', index=False)
    (fig, axes) = plt.subplots(1, 3, figsize=(13, 4))
    for (ax, param, xlabel) in zip(axes, ['Surface_tension', 'Contact_angle', 'Wetting_rate'], ['Surface tension (mN/m)', 'Contact angle (deg)', 'Wetting rate (1/s)']):
        sub = out[out['Parameter'] == param]
        ax.plot(sub['Value'], sub['Final_uptake'], 'o-', color='#E63946', label='Uptake')
        ax2 = ax.twinx()
        ax2.plot(sub['Value'], sub['Final_integrity'], 's--', color='#2A9D8F', label='Integrity')
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Final uptake')
        ax2.set_ylabel('Final integrity')
        ax.set_title(param)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig4_Sensitivity_Analysis.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out
 
def run_knockout():
    crispr = {'SNARE_KO': (0.35, 0.08), 'ESCRT_KO': (1.45, 0.12), 'Glyco_KO': (0.55, 0.1), 'PIK3C_KO': (1.6, 0.15)}
    wt = simulate('micro')
    wt_last = wt.iloc[-1]
    wt_uptake = float(wt_last['BMC_cytoplasm'])
    rows = [{'Condition': 'WT', 'Final_uptake': wt_uptake, 'Final_integrity': float(wt_last['Membrane_integrity']), 'Max_SNARE': float(wt['SNARE_activity'].max()), 'Max_ESCRT': float(wt['ESCRT_activity'].max()), 'Max_Glyco': float(wt['Glycoprotein_activity'].max()), 'Max_PIK3C': float(wt['PIK3C_activity'].max()), 'GP_change': float(wt_last['GP_value'] - wt.iloc[0]['GP_value']), 'Uptake_ratio': 1.0, 'CRISPR_expected': np.nan, 'CRISPR_std': np.nan, 'Model_vs_CRISPR_error': np.nan}]
    ko_map = {'SNARE_KO': ('snare', 0.18), 'ESCRT_KO': ('escrt', 0.18), 'Glyco_KO': ('glyco', 0.12), 'PIK3C_KO': ('pik3c', 0.18)}
    for (cond, (mod_name, scale)) in ko_map.items():
        df = simulate('micro', module_overrides={mod_name: {'k_act': MODULE_PARAMS[mod_name]['k_act'] * scale}})
        last = df.iloc[-1]
        ratio = float(last['BMC_cytoplasm']) / wt_uptake
        (exp_mean, exp_std) = crispr[cond]
        rows.append({'Condition': cond, 'Final_uptake': float(last['BMC_cytoplasm']), 'Final_integrity': float(last['Membrane_integrity']), 'Max_SNARE': float(df['SNARE_activity'].max()), 'Max_ESCRT': float(df['ESCRT_activity'].max()), 'Max_Glyco': float(df['Glycoprotein_activity'].max()), 'Max_PIK3C': float(df['PIK3C_activity'].max()), 'GP_change': float(last['GP_value'] - df.iloc[0]['GP_value']), 'Uptake_ratio': ratio, 'CRISPR_expected': exp_mean, 'CRISPR_std': exp_std, 'Model_vs_CRISPR_error': abs(ratio - exp_mean)})
    out = pd.DataFrame(rows)
    out.to_excel('outputs/pathway_knockout_analysis.xlsx', index=False)
    (fig, ax) = plt.subplots(figsize=(7, 4.5))
    sub = out[out['Condition'] != 'WT']
    x = np.arange(len(sub))
    ax.bar(x - 0.15, sub['Uptake_ratio'], width=0.3, label='Model', color='#264653')
    ax.bar(x + 0.15, sub['CRISPR_expected'], width=0.3, yerr=sub['CRISPR_std'], label='CRISPR', color='#E9C46A', capsize=4)
    ax.axhline(1.0, color='gray', ls='--', lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(sub['Condition'])
    ax.set_ylabel('Uptake ratio vs WT')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig6_Pathway_Knockout_Analysis.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out
 
def run_other_materials():
    materials = [('BMC_micro', 25, 120, 1500, 0.05), ('BMC_nano', 15, 60, 100, 0.1), ('Oleic_acid', 32, 140, 1000, 0.03), ('Squalene', 28, 130, 1200, 0.04), ('LNP', 12, 50, 80, 0.12), ('PLGA_emulsion', 20, 100, 500, 0.07)]
    rows = []
    for (name, gamma_mNm, theta, radius_nm, k_wet) in materials:
        label = 'nano' if radius_nm < 200 else 'micro'
        df = simulate(label, size_overrides={'gamma0': gamma_mNm / 1000.0, 'theta0': float(theta), 'R': radius_nm * 1e-09, 'k_wet': float(k_wet)})
        last = df.iloc[-1]
        rows.append({'Material': name, 'Surface_tension_mN_m': gamma_mNm, 'Contact_angle_deg': theta, 'Radius_nm': radius_nm, 'Wetting_rate': k_wet, 'Final_uptake': float(last['BMC_cytoplasm']), 'Final_integrity': float(last['Membrane_integrity']), 'GP_change': float(last['GP_value'] - df.iloc[0]['GP_value']), 'Max_wetting': float(df['Wetting_degree'].max())})
    out = pd.DataFrame(rows)
    out.to_excel('outputs/other_droplet_materials.xlsx', index=False)
    (fig, ax) = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(out))
    ax.bar(x - 0.15, out['Final_uptake'], width=0.3, label='Uptake', color='#E63946')
    ax.bar(x + 0.15, out['Final_integrity'], width=0.3, label='Integrity', color='#2A9D8F')
    ax.set_xticks(x)
    ax.set_xticklabels(out['Material'], rotation=20, ha='right')
    ax.set_ylabel('Fraction')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title('Other droplet materials')
    plt.tight_layout()
    plt.savefig('outputs/Fig5_Other_Droplet_Materials.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return out
 
def run_monte_carlo(n=100, seed=42):
    rng = np.random.default_rng(seed)
    keys = ['gamma0', 'theta0', 'k_wet', 'R', 'k_fus', 'k_endo', 'k_dmg', 'E_c']
    sheets = {}
    summary_rows = []
    for (label, bmc_type, sheet_name) in [('micro', 'Micro-BMC', 'Micro_BMC_simulations'), ('nano', 'Nano-BMC', 'Nano_BMC_simulations')]:
        base = SIZE_PARAMS[label]
        rows = []
        for i in range(n):
            ov = {}
            for k in keys:
                ov[k] = base[k] * float(np.clip(rng.normal(1.0, 0.2), 0.4, 1.6))
            df = simulate(label, size_overrides=ov)
            last = df.iloc[-1]
            rows.append({'simulation': i, 'final_uptake': float(last['BMC_cytoplasm']), 'final_integrity': float(last['Membrane_integrity']), 'max_wetting': float(df['Wetting_degree'].max()), 'GP_change': float(last['GP_value'] - df.iloc[0]['GP_value'])})
        sim = pd.DataFrame(rows)
        sheets[sheet_name] = sim
        up = sim['final_uptake'].values
        integ = sim['final_integrity'].values
        summary_rows.append({'BMC_type': bmc_type, 'Uptake_mean': float(np.mean(up)), 'Uptake_std': float(np.std(up, ddof=1)), 'Uptake_95CI_lower': float(np.percentile(up, 2.5)), 'Uptake_95CI_upper': float(np.percentile(up, 97.5)), 'Integrity_mean': float(np.mean(integ)), 'Integrity_std': float(np.std(integ, ddof=1)), 'Integrity_95CI_lower': float(np.percentile(integ, 2.5)), 'Integrity_95CI_upper': float(np.percentile(integ, 97.5))})
    sheets['Summary_statistics'] = pd.DataFrame(summary_rows)
    with pd.ExcelWriter('outputs/monte_carlo_uncertainty.xlsx') as wb:
        for (name, df) in sheets.items():
            df.to_excel(wb, sheet_name=name, index=False)
    (fig, axes) = plt.subplots(1, 2, figsize=(10, 4))
    for (ax, key, title) in zip(axes, ['Micro_BMC_simulations', 'Nano_BMC_simulations'], ['Micro-BMC', 'Nano-BMC']):
        ax.hist(sheets[key]['final_uptake'], bins=20, color='#457B9D', alpha=0.85)
        ax.set_xlabel('Final uptake')
        ax.set_ylabel('Count')
        ax.set_title(title)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('outputs/Fig8_Monte_Carlo_Uncertainty.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return sheets['Summary_statistics']
 
def save_parameter_references():
    gamma = pd.DataFrame([{'Material': 'micro_BMC', 'Value': 25, 'Unit': 'mN/m', 'Source': 'Estimated from lipid droplet γ=20-30 mN/m (Thiam et al. Nat Rev Mol Cell Biol 2013)'}, {'Material': 'nano_BMC', 'Value': 15, 'Unit': 'mN/m', 'Source': 'Lower γ for smaller droplets (Tolman length effect, Lei et al. PNAS 2017)'}])
    angle = pd.DataFrame([{'Material': 'micro_BMC', 'Value': 120, 'Unit': 'deg', 'Source': 'Hydrophobic droplet on bilayer (Snead et al. Cell 2017)'}, {'Material': 'nano_BMC', 'Value': 60, 'Unit': 'deg', 'Source': 'Size-dependent wettability (Schellenberger et al. PRL 2016)'}])
    crispr = pd.DataFrame([{'Pathway': 'SNARE', 'Expected_uptake_ratio': 0.35, 'Std': 0.08, 'Role': 'Fusion promoter'}, {'Pathway': 'ESCRT', 'Expected_uptake_ratio': 1.45, 'Std': 0.12, 'Role': 'Repair / uptake suppressor'}, {'Pathway': 'Glycoprotein', 'Expected_uptake_ratio': 0.55, 'Std': 0.1, 'Role': 'Wetting promoter'}, {'Pathway': 'PIK3C', 'Expected_uptake_ratio': 1.6, 'Std': 0.15, 'Role': 'Membrane rigidity / suppressor'}])
    with pd.ExcelWriter('outputs/parameter_references.xlsx') as wb:
        gamma.to_excel(wb, sheet_name='Surface_tension', index=False)
        angle.to_excel(wb, sheet_name='Contact_angle', index=False)
        crispr.to_excel(wb, sheet_name='CRISPR_calibration', index=False)
 
def main():
    print('sensitivity')
    print(run_sensitivity().head())
    print('knockout')
    print(run_knockout())
    print('materials')
    print(run_other_materials())
    print('monte carlo')
    print(run_monte_carlo())
    save_parameter_references()
    print('done')
if __name__ == '__main__':
    main()

