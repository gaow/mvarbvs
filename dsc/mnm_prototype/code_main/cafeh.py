import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

from cafeh.cafeh_summary import fit_cafeh_summary, fit_cafeh_z
from cafeh.cafeh_genotype import fit_cafeh_genotype
    
from cafeh.model_queries import *

LD_df = pd.read_csv(LD,sep='\t',index_col=0)
beta_df = pd.read_csv(bhat,sep='\t',index_col=0)
stderr_df = pd.read_csv(shat,sep='\t',index_col=0)
n_df = pd.read_csv(n_file,sep='\t',index_col=0)
cafeh = fit_cafeh_summary(LD_df, beta_df, stderr_df, n=n_df)

cafeh.save(m_name + '.model', save_data=True, save_ld=True)
variant_report = summary_table(cafeh)
variant_report.to_csv(m_name + '.result', sep='\t')

credible_sets = cafeh.credible_sets
purity = cafeh.purity

# filtered_cs = {k:v for k,v in credible_sets.items() if purity[k] > 0.5}
# filtered_purity = {k:v for k,v in purity.items() if purity[k] > 0.5}
# with open(out_cs_name,'wt') as out_h:
#   json_object = json.dumps(filtered_cs)
#   out_h.write(json_object)
# with open(out_purity_name,'wt') as out_h:
#   json_object = json.dumps(filtered_purity)
#   out_h.write(json_object)
  
pip = cafeh.get_pip()

trait_pip = cafeh.get_study_pip()
