import sys
sys.modules[__name__].__dict__.clear()
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

samples = np.array(('5953STDY8551191', '5953STDY8551192', '5953STDY8551193', '5953STDY8551194',
                   '5953STDY8551195', '5953STDY8551196', '5953STDY8551197', '5953STDY8551198'))

for i in range(len(samples)):
    print(i)
    output_dir = '/home/jovyan/HB_ZIK/scrublet/'
    input_dir = '/home/jovyan/data/ZikaGlioblastomas/cellranger302_count_32771_' + samples[i] + '_GRCh38-3_0_0_premrna/filtered_feature_bc_matrix/filtered_feature_bc_matrix'
    counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
    genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))

    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)

    predicted_doublets = predicted_doublets*1
    predicted_doublets = predicted_doublets.astype(int)
    detected_doublets_rate = round(scrub.detected_doublet_rate_, 4)
    overall_doublets_rate = round(scrub.overall_doublet_rate_, 4)

    np.savetxt(output_dir + '/' + samples[i] + '_' + 'doublets_scores.txt', doublet_scores)   
    np.savetxt(output_dir +  '/' + samples[i] + '_' + 'predicted_doublets.txt', predicted_doublets)                              
    with open(output_dir +  '/' + samples[i] + '_' + 'detected_doublets_rate.txt', 'w') as f:
      f.write('%f' % detected_doublets_rate)  

    with open(output_dir + '/' + samples[i] + '_' +  'overall_doublets_rate.txt', 'w') as f:
      f.write('%f' % overall_doublets_rate)

    f = scrub.plot_histogram()
    f[0].savefig(output_dir + '/' + samples[i] + '_' + "doubletScore_histogram.pdf", bbox_inches='tight')
    
