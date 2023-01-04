import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

from cafeh.cafeh_summary import fit_cafeh_summary, fit_cafeh_z
from cafeh.cafeh_genotype import fit_cafeh_genotype
    
from cafeh.model_queries import *
import pickle

cafehg_reload = pickle.load(open(cafeh_file + '.model', 'rb'))

thold = [0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.5,0.4,0.3,0.2,0.1]

sets = {}
purity = {}

for i in thold:
  sets[str(i)] = cafehg_reload.get_credible_sets(i)
  purity[str(i)] = cafehg_reload.get_purity(i)
  
