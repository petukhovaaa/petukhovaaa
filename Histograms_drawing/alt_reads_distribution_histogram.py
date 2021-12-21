import pandas as pd
from io import StringIO
import sys
import matplotlib.pyplot as plt


path = sys.argv[1]
vcf_file = sys.argv[2]
output_name = sys.argv[3]


# Step 1: finding the header string in the input vcf-file and opening the file

with open(path + vcf_file) as f:
    h = 0
    for line in f:
        if not line.startswith('#CHROM'):
            h+=1
        else:
            break # we found header of the input vcf-file

read_vcf = pd.read_csv(path + vcf_file, sep='\t', header = h)


# Step 2: parsing the pool column and creating column with fraction of alt reads

read_vcf['AD'] = read_vcf.pool.apply(lambda p: p.split(':')[1])
read_vcf['DP'] = read_vcf.pool.apply(lambda p: p.split(':')[3])
read_vcf['AD_alt'] = read_vcf.AD.apply(lambda p: p.split(',')[1])
read_vcf.loc[read_vcf['DP'] == '0', 'DP'] = 1 # DP normalisation
read_vcf['AD_alt_AF'] = read_vcf.apply(lambda p: int(p['AD_alt'])/int(p['DP']), axis = 1)


# Step 3: picture drawing

fig, ax = plt.subplots()

fig.set_size_inches(18.5, 10.5)
ax.hist(read_vcf_exome.loc[read_vcf_exome['DP'] >= 200, 'AD_alt_AF'].sort_values(), len(set(read_vcf_exome['AD_alt_AF']))//100, linewidth=0.6, edgecolor="white")
ax.set_xlabel('Fraction of reads wilt alt alleles')
ax.set_ylabel('Quantity of the reads')
ax.set_title(r'Histogram of fraction of reads with alt alleles')


plt.savepig(path + ouptut_name + '.png')
