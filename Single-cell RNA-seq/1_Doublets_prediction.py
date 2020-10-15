import os
os.getcwd()
os.chdir('/path/to/0_doublets_detection_scrublet')
import sys
print(sys.path)

#srublet usage
#%matplotlib inline
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42


############input 10X matrix
input_dir = '/path/to/filtered_feature_bc_matrix/'
counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=1))
out_df = pd.read_csv(input_dir + 'barcodes.tsv.gz', header = None, index_col=None, names=['barcode'])

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))


############initiate Scrublet object
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
############calculating doublet score
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)

############histogram of double score
scrub.call_doublets(threshold=0.25)
scrub.plot_histogram()
plt.savefig("INTS11_sample_double_score_histogram.pdf")


############dimensionality reduction
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')
scrub.plot_embedding('UMAP', order_points=True)
plt.savefig("INTS11_sample_double_score_umap.pdf")
print (scrub.detected_doublet_rate_)

############output doublet result
out_df['doublet_scores'] = doublet_scores
out_df['predicted_doublets'] = predicted_doublets
out_df.to_csv(input_dir + '/sample_doublet.txt', index=False,header=True)
out_df.head()

