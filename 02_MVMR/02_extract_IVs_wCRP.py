# ================================================================================
# MVMR Data Prep v3 - Including CRP (with capped instruments)
# Strategy: Limit CRP to top 200 most significant independent SNPs
# ================================================================================
import gzip, sys, os, csv
sys.stdout.reconfigure(encoding='utf-8')

PVAL_THRESH = 5e-8
CLUMP_KB = 10000
MAX_CRP_SNPS = 200  # Cap CRP instruments

DATA_4TH = r"F:\MR\Interleukin-6 levels and keloids\数据\第四次"
DATA_MVMR = r"F:\MR\Interleukin-6 levels and keloids\数据\Multivariable MR 数据"

exposures = {
    "IL6":   os.path.join(DATA_4TH, "ebi-a-GCST90012005.vcf.gz"),
    "IL1RA": os.path.join(DATA_MVMR, "ebi-a-GCST90012004.vcf.gz"),
    "TNFR1": os.path.join(DATA_MVMR, "ebi-a-GCST90012015.vcf.gz"),
    "CRP":   os.path.join(DATA_MVMR, "ebi-a-GCST90029070.vcf.gz")
}
outcome_vcf = os.path.join(DATA_4TH, "out-ebi-a-GCST90018874.vcf.gz")

def parse_vcf_line(line):
    parts = line.strip().split('\t')
    if len(parts) < 9:
        return None
    fmt_keys = parts[8].split(':')
    vals_list = parts[9].split(':')
    fmt_dict = dict(zip(fmt_keys, vals_list))
    return {
        'CHR': parts[0], 'POS': int(parts[1]), 'ID': parts[2],
        'REF': parts[3], 'ALT': parts[4],
        'ES': float(fmt_dict.get('ES', 0)),
        'SE': float(fmt_dict.get('SE', 1)),
        'LP': float(fmt_dict.get('LP', 0)),
        'AF': float(fmt_dict.get('AF', 0)) if 'AF' in fmt_dict else None,
    }

def extract_significant_snps(vcf_path, pval_thresh=5e-8, max_snps=None):
    snps = []
    with gzip.open(vcf_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = parse_vcf_line(line)
            if p is None:
                continue
            pval = 10 ** (-p['LP'])
            if pval < pval_thresh:
                snps.append({
                    'SNP': p['ID'], 'CHR': p['CHR'], 'BP': p['POS'],
                    'beta': p['ES'], 'se': p['SE'], 'pval': pval,
                    'effect_allele': p['ALT'].split(',')[0] if ',' in p['ALT'] else p['ALT'],
                    'other_allele': p['REF'],
                    'eaf': p['AF'] if p['AF'] is not None else 0,
                    'neg_log10_p': p['LP'],
                })
    return snps

def extract_snp_effects_fast(vcf_path, snp_set):
    results = {}
    with gzip.open(vcf_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = parse_vcf_line(line)
            if p is None or p['ID'] not in snp_set:
                continue
            pval = 10 ** (-p['LP'])
            results[p['ID']] = {
                'beta': p['ES'], 'se': p['SE'], 'pval': pval,
                'ea': p['ALT'].split(',')[0] if ',' in p['ALT'] else p['ALT'],
                'oa': p['REF'],
                'eaf': p['AF'] if p['AF'] is not None else 0,
            }
    return results

print("=" * 60)
print("MVMR Data Prep v3 - Including CRP (capped)")
print("=" * 60)

# Step 1: Extract & clump for each exposure
print("\n--- Step 1: Extract & clump instruments ---")
all_clumped = {}
all_ivs = set()

for exp_name, vcf_path in exposures.items():
    print(f"\nProcessing {exp_name}...")
    
    if exp_name == "CRP":
        # For CRP: use lower threshold to get manageable number
        # Try P < 5e-100 first to get ~top SNPs
        # Actually, let's just extract all significant and take top
        snps = extract_significant_snps(vcf_path, pval_thresh=5e-8)
        print(f"  All significant SNPs (P<5e-8): {len(snps):,}")
    else:
        snps = extract_significant_snps(vcf_path)
        print(f"  Significant SNPs: {len(snps)}")
    
    if len(snps) == 0:
        all_clumped[exp_name] = []
        continue
    
    snps.sort(key=lambda x: x['pval'])
    
    # Distance-based clumping
    clumped = []
    last_chr = None
    last_bp = None
    for s in snps:
        if (last_chr is None or s['CHR'] != last_chr or 
            abs(s['BP'] - last_bp) > CLUMP_KB * 1000):
            clumped.append(s)
            last_chr = s['CHR']
            last_bp = s['BP']
    
    # Cap CRP instruments
    if exp_name == "CRP":
        print(f"  After clumping: {len(clumped):,} SNPs")
        clumped = clumped[:MAX_CRP_SNPS]
        print(f"  After capping to top {MAX_CRP_SNPS}: {len(clumped)} SNPs")
    else:
        print(f"  After clumping: {len(clumped)} SNPs")
    
    chr_counts = {}
    for s in clumped:
        chr_counts[s['CHR']] = chr_counts.get(s['CHR'], 0) + 1
    print(f"  Chr dist: {chr_counts}")
    
    all_clumped[exp_name] = clumped
    for s in clumped:
        all_ivs.add(s['SNP'])

print(f"\nTotal unique IV SNPs: {len(all_ivs):,}")

# Step 2: Extract effects for union SNPs
print("\n--- Step 2: Extract SNP effects ---")
exposure_effects = {}
for exp_name, vcf_path in exposures.items():
    print(f"  {exp_name}...", end=' ', flush=True)
    exposure_effects[exp_name] = extract_snp_effects_fast(vcf_path, all_ivs)
    print(f"{len(exposure_effects[exp_name])}/{len(all_ivs)} matched")

print(f"  Keloid (outcome)...", end=' ', flush=True)
outcome_effects = extract_snp_effects_fast(outcome_vcf, all_ivs)
print(f"{len(outcome_effects)}/{len(all_ivs)} matched")

# Step 3: Build unified dataset
print("\n--- Step 3: Build unified dataset ---")
merged_rows = []
for snp in sorted(all_ivs):
    row = {'SNP': snp}
    complete = True
    for exp_name in exposures.keys():
        if snp in exposure_effects.get(exp_name, {}):
            e = exposure_effects[exp_name][snp]
            row[f'beta_{exp_name}'] = e['beta']
            row[f'se_{exp_name}'] = e['se']
            row[f'ea_{exp_name}'] = e['ea']
            row[f'oa_{exp_name}'] = e['oa']
            row[f'eaf_{exp_name}'] = e['eaf']
            row[f'pval_{exp_name}'] = e['pval']
        else:
            complete = False
            break
    if snp in outcome_effects:
        o = outcome_effects[snp]
        row['beta_keloid'] = o['beta']
        row['se_keloid'] = o['se']
        row['ea_keloid'] = o['ea']
        row['oa_keloid'] = o['oa']
        row['eaf_keloid'] = o['eaf']
        row['pval_keloid'] = o['pval']
    else:
        complete = False
    if complete:
        merged_rows.append(row)

print(f"\n  Complete cases: {len(merged_rows)}/{len(all_ivs)}")

# Also save IL6+CRP only (minimal MVMR)
print("\n--- Also saving IL6+CRP-only subset ---")
il6_crp_snps = set()
for s in all_clumped.get("IL6", []):
    il6_crp_snps.add(s['SNP'])
for s in all_clumped.get("CRP", []):
    il6_crp_snps.add(s['SNP'])

il6crp_rows = []
for snp in sorted(il6_crp_snps):
    row = {'SNP': snp}
    complete = True
    for exp_name in ['IL6', 'CRP']:
        if snp in exposure_effects.get(exp_name, {}):
            e = exposure_effects[exp_name][snp]
            row[f'beta_{exp_name}'] = e['beta']
            row[f'se_{exp_name}'] = e['se']
            row[f'ea_{exp_name}'] = e['ea']
            row[f'oa_{exp_name}'] = e['oa']
            row[f'eaf_{exp_name}'] = e['eaf']
            row[f'pval_{exp_name}'] = e['pval']
        else:
            complete = False
            break
    if snp in outcome_effects:
        o = outcome_effects[snp]
        row['beta_keloid'] = o['beta']
        row['se_keloid'] = o['se']
        row['ea_keloid'] = o['ea']
        row['oa_keloid'] = o['oa']
        row['eaf_keloid'] = o['eaf']
        row['pval_keloid'] = o['pval']
    else:
        complete = False
    if complete:
        il6crp_rows.append(row)
print(f"  IL6+CRP complete cases: {len(il6crp_rows)}/{len(il6_crp_snps)}")

# Step 4: Save CSVs
print("\n--- Step 4: Saving ---")

# Full model (4 exposures)
csv_full = os.path.join(DATA_MVMR, "MVMR_input_wCRP.csv")
fn_full = ['SNP']
for e in exposures.keys():
    for p in ['beta','se','ea','oa','eaf','pval']:
        fn_full.append(f'{p}_{e}')
for p in ['beta','se','ea','oa','eaf','pval']:
    fn_full.append(f'{p}_keloid')

with open(csv_full, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fn_full)
    w.writeheader()
    w.writerows(merged_rows)
print(f"  Full model (4 exposures): {len(merged_rows)} SNPs → MVMR_input_wCRP.csv")

# IL6+CRP only
csv_min = os.path.join(DATA_MVMR, "MVMR_input_IL6_CRP.csv")
fn_min = ['SNP']
for e in ['IL6','CRP']:
    for p in ['beta','se','ea','oa','eaf','pval']:
        fn_min.append(f'{p}_{e}')
for p in ['beta','se','ea','oa','eaf','pval']:
    fn_min.append(f'{p}_keloid')

with open(csv_min, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fn_min)
    w.writeheader()
    w.writerows(il6crp_rows)
print(f"  Minimal model (IL6+CRP): {len(il6crp_rows)} SNPs → MVMR_input_IL6_CRP.csv")

print("\nDone!")
