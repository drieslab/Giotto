import scrublet as scr

def python_scrublet(counts_matrix, expected_doublet_rate, min_counts, min_cells, min_gene_variability_pctl, n_prin_comps):
  
  min_counts = int(min_counts)
  min_cells = int(min_cells)
  min_gene_variability_pctl = int(min_gene_variability_pctl)
  n_prin_comps = int(n_prin_comps)
  
  scrub = scr.Scrublet(counts_matrix=counts_matrix, expected_doublet_rate=expected_doublet_rate)
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=min_counts,min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps)
  
  return_list = []
  return_list.append(doublet_scores)
  return_list.append(predicted_doublets)
  
  return(return_list)
