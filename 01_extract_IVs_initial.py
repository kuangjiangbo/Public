# ================================================================================
# MVMR Data Preparation - Extract & Clump Instruments from VCF files
# Uses Python for efficient VCF processing, then R for MVMR analysis
# ================================================================================
import gzip, sys, os, math, csv
sys.stdout.reconfigure(encoding='utf-8')

PVAL_THRESH = 5e-8
CLUMP_KB = 10000

DATA_4TH = r"F:\MR\Interleukin-6 levels and keloids\数据\第四次"
DATA_MVMR = r"F:\MR\Interleukin-6 levels and keloids\数据\Multivariable MR 数据"

exposures = {
    "IL6":   os.path.join(DATA_4TH, "ebi-a-GCST90012005.vcf.gz"),
    "IL1RA": os.path.join(DATA_MVMR, "ebi-a-GCST90012004.vcf.gz"),
    "TNFR1": os.path.join(DATA_MVMR, "ebi-a-GCST90012015.vcf.gz")
}
outcome_vcf = os.path.join(DATA_4TH, "out-ebi-a-GCST90018874.vcf.gz")

def extract_significant_snps(vcf_path, pval_thresh=5e-8):
    """Extract genome-wide significant SNPs from VCF file"""
    snps = []
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            snp_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            
            fmt = parts[8]
            vals = parts[9]
            
            fmt_keys = fmt.split(':')
            vals_list = vals.split(':')
            fmt_dict = dict(zip(fmt_keys, vals_list))
            
            es = float(fmt_dict.get('ES', 0))
            se = float(fmt_dict.get('SE', 1))
            lp = float(fmt_dict.get('LP', 0))
            af = float(fmt_dict.get('AF', 0))
            
            pval = 10 ** (-lp)
            
            if pval < pval_thresh:
                snps.append({
                    'SNP': snp_id, 'CHR': chrom, 'BP': pos,
                    'beta': es, 'se': se, 'pval': pval,
                    'effect_allele': alt.split(',')[0] if ',' in alt else alt,
                    'other_allele': ref, 'eaf': af
                })
    return snps

def extract_snp_effects(vcf_path, snp_set):
    """Extract effects for specific SNPs from VCF"""
    results = {}
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            snp_id = parts[2]
            if snp_id not in snp_set:
                continue
            
            ref = parts[3]
            alt = parts[4]
            fmt = parts[8]
            vals = parts[9]
            
            fmt_keys = fmt.split(':')
            vals_list = vals.split(':')
            fmt_dict = dict(zip(fmt_keys, vals_list))
            
            results[snp_id] = {
                'SNP': snp_id, 'CHR': parts[0], 'BP': int(parts[1]),
                'beta': float(fmt_dict.get('ES', 0)),
                'se': float(fmt_dict.get('SE', 1)),
                'pval': 10 ** (-float(fmt_dict.get('LP', 0))),
                'ea': alt.split(',')[0] if ',' in alt else alt,
                'oa': ref,
                'eaf': float(fmt_dict.get('AF', 0))
            }
    return results

print("=" * 60)
print("MVMR Data Preparation")
print("=" * 60)

# Step 1: Extract significant SNPs from each exposure
print("\n--- Step 1: Extracting significant SNPs ---")
all_clumped = {}
all_ivs = set()

for exp_name, vcf_path in exposures.items():
    print(f"\nProcessing {exp_name}...")
    snps = extract_significant_snps(vcf_path)
    print(f"  Significant SNPs: {len(snps)}")
    
    if len(snps) == 0:
        all_clumped[exp_name] = []
        continue
    
    snps.sort(key=lambda x: x['pval'])
    
    clumped = []
    last_chr = None
    last_bp = None
    
    for s in snps:
        if last_chr is None or s['CHR'] != last_chr or abs(s['BP'] - last_bp) > CLUMP_KB * 1000:
            clumped.append(s)
            last_chr = s['CHR']
            last_bp = s['BP']
    
    print(f"  After clumping: {len(clumped)} SNPs")
    
    chr_counts = {}
    for s in clumped:
        chr_counts[s['CHR']] = chr_counts.get(s['CHR'], 0) + 1
    print(f"  Chromosome distribution: {chr_counts}")
    
    all_clumped[exp_name] = clumped
    for s in clumped:
        all_ivs.add(s['SNP'])

print(f"\nTotal unique IV SNPs: {len(all_ivs)}")

# Step 2: Extract effects for all IV SNPs from each exposure and outcome
print("\n--- Step 2: Extracting effects for union IV SNPs ---")

exposure_effects = {}
for exp_name, vcf_path in exposures.items():
    print(f"  Extracting {exp_name}...")
    exposure_effects[exp_name] = extract_snp_effects(vcf_path, all_ivs)
    print(f"    Matched: {len(exposure_effects[exp_name])}/{len(all_ivs)}")

print(f"  Extracting outcome (Keloid)...")
outcome_effects = extract_snp_effects(outcome_vcf, all_ivs)
print(f"    Matched: {len(outcome_effects)}/{len(all_ivs)}")

# Step 3: Build unified dataframe
print("\n--- Step 3: Building unified dataset ---")
all_snps_sorted = sorted(all_ivs)
merged_rows = []

for snp in all_snps_sorted:
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

print(f"  Complete cases (union IVs with all data): {len(merged_rows)}/{len(all_snps_sorted)}")

# Step 4: Save CSV for R
output_csv = os.path.join(DATA_MVMR, "MVMR_input_data.csv")
fieldnames = ['SNP']
for e in exposures.keys():
    for p in ['beta','se','ea','oa','eaf','pval']:
        fieldnames.append(f'{p}_{e}')
for p in ['beta','se','ea','oa','eaf','pval']:
    fieldnames.append(f'{p}_keloid')

with open(output_csv, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(merged_rows)

print(f"\nSaved: {output_csv}")
print(f"  {len(merged_rows)} SNPs, {len(exposures)+1} traits")
print("\nData preparation complete!")
