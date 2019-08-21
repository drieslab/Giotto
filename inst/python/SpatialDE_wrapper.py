import numpy as np
import pandas as pd
import sys
import NaiveDE
import SpatialDE

def Spatial_DE(filterd_exprs, coordinates):
    if(filterd_exprs.shape[0] != coordinates.shape[0]):
        sys.exit("The number of cells in expression file and location file don't match\n")
    else:
        ## results and ms_results
        coordinates_cp = coordinates.copy()
        coordinates_cp['total_counts'] = filterd_exprs.sum(1)
        
        dfm = NaiveDE.stabilize(filterd_exprs.T).T
        res = NaiveDE.regress_out(coordinates_cp, dfm.T, 'np.log(total_counts)').T
        res['log_total_count'] = np.log(coordinates_cp['total_counts'])
        
        results = SpatialDE.run(coordinates, res)
        
        de_results = results[(results.qval < 0.05)].copy()
        if(de_results.shape[0] > 0):
            ms_results = SpatialDE.model_search(coordinates, res, de_results)
            result_dic = {"results":results, "ms_results":ms_results}
        
        else:
            print("No spatially variable genes found! \n")
            result_dic = {"results":results}
        
        return result_dic
    
    
def Spatial_DE_AEH(filterd_exprs,coordinates,results,pattern_num,l = 1.05, verbosity = 1):
    ## Automatic expression histology
        coordinates_cp =coordinates.copy()
        coordinates_cp['total_counts'] = filterd_exprs.sum(1)
        
        dfm = NaiveDE.stabilize(filterd_exprs.T).T
        res = NaiveDE.regress_out(coordinates_cp, dfm.T, 'np.log(total_counts)').T
        
        results['pval'] = results['pval'].clip(lower = results.query('pval > 0')['pval'].min() / 2)
        results['qval'] = results['qval'].clip(lower = results.query('qval > 0')['qval'].min() / 2)

        sres = results.query('qval < 0.05 & g != "log_total_count"').copy()
        
        X = coordinates.values

        histology_results, patterns = SpatialDE.spatial_patterns(X, res, sres, int(pattern_num), l = l,verbosity=verbosity)
        
        pattern_dic = {"histology_results":histology_results,"patterns":patterns}
        
        return pattern_dic
