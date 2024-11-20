import sys
sys.path.insert(0, 'manual_load_enformer_scripts/')
import scoring
import numpy as np

class common_scoring_methods:
    
    '''
    Define the function for common scorings 
    '''
    
    def __init__(self, map_a, map_b):
        self.map_a = map_a
        self.map_b = map_b
    
    def mse(self): # mean squared error
        return scoring.mse(self.map_a, self.map_b)
    
    def corr(self): # spearman correlation
        return scoring.spearman(self.map_a, self.map_b)
    
    def ssi(self): # structural similarity index
        return scoring.ssim_map(self.map_a, self.map_b)
    
    def scc(self): # stratum adjusted correlation
        return scoring.scc(self.map_a, self.map_b)
 
