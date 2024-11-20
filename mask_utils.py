import numpy as np
import pandas as pd
import sys
sys.path.insert(0, 'manual_load_enformer_scripts/')
from bin_utils import get_bin

def mask_matrices(REF_pred, ALT_pred, SVTYPE, SVLEN, var_rel_pos, bin_size, offset, target_length_cropped):

    '''
    This applied to non-BND predicted matrices.
    
    Mask reference and alternate predicted matrices based on the type of variant.
    
    '''
    if SVTYPE in ['DEL', 'DUP', "INS"]:
        
        # Insertions: Add nans to reference matrix and crop ends, then mirror nans to alternate matrix
        
        if SVTYPE == 'DEL':
            
            # Deletions: Same but swapping reference and alternate values
            REF_pred, ALT_pred = ALT_pred, REF_pred
            var_rel_pos.reverse()
            
            
        # Get variant relative positions in the unmasked predicted maps
        var_start_REF = get_bin(var_rel_pos[0], bin_size, offset)

        var_start_ALT = get_bin(var_rel_pos[1], bin_size, offset)
        var_end_ALT = get_bin(var_rel_pos[1] + SVLEN, bin_size, offset)

        

        #### Mask REF map
        
        REF_pred_masked = REF_pred.copy()
        ALT_pred_masked = ALT_pred.copy()

        
        # Insert the empty bins in the reference matrix where the variant is in the alternate matrix
        to_add = var_end_ALT - var_start_ALT

        for j in range(var_start_REF, var_start_REF + to_add): # range only includes the first variable 
            REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 0)
            # REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 1)

            
        # Crop the outside of the reference matrix 
        to_remove_left = var_start_REF - var_start_ALT
        to_remove_right = len(REF_pred_masked) - target_length_cropped - to_remove_left

        if to_remove_left != 0:
            REF_pred_masked = REF_pred_masked[to_remove_left:]
        if to_remove_right != 0:
            REF_pred_masked = REF_pred_masked[:-to_remove_right]



        
        if SVTYPE == 'DEL':
            
            # Deletions: Swap reference and alternate values back
            REF_pred_masked, ALT_pred_masked = ALT_pred_masked, REF_pred_masked
        
    # Inversions: Mask REF and mirror nans to ALT
    elif SVTYPE == 'INV':
        
        var_start = get_bin(var_rel_pos[0], bin_size, offset)
        var_end = get_bin(var_rel_pos[0] + SVLEN, bin_size, offset)
 

        # Mask REF map: make variant bin(s) nan 
        
        REF_pred_masked = REF_pred.copy()

        REF_pred_masked[var_start:var_end + 1] = np.nan
        # REF_pred_masked[:, var_start:var_end + 1] = np.nan

        
        # Mask ALT map: make all nans in REF_pred also nan in ALT_pred
        
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        
        
    if len(REF_pred_masked) != target_length_cropped:
        # print(len(REF_pred_masked))
        # print(target_length_cropped)
        raise ValueError('Masked reference matrix is not the right size.')    
        
        
    if len(ALT_pred_masked) != target_length_cropped:
        raise ValueError('Masked alternate matrix is not the right size.')
    
    
    return REF_pred_masked, ALT_pred_masked

def get_masked_BND_maps(matrices, rel_pos_map, target_length_cropped):

    # Get REF and ALT vectors, excluding diagonal
    indexes_left = rel_pos_map
    indexes_right = target_length_cropped - rel_pos_map

    REF_L = get_left_BND_map(matrices[0], rel_pos_map)
    REF_R = get_right_BND_map(matrices[1], rel_pos_map)

    return (assemple_BND_maps(REF_L, REF_R, rel_pos_map, target_length_cropped),matrices[2])

def get_left_BND_map(pred_matrix, rel_pos_map):

    '''
    Take upper left quarter (or more or less if shifted) of the matrix.

    '''
    
    left_BND_map = pred_matrix[:int(rel_pos_map)]
    return left_BND_map
    
    
    
def get_right_BND_map(pred_matrix, rel_pos_map):
    
    '''
    Take lower right quarter (or more or less if shifted) of the matrix.

    ''' 
    right_BND_map = pred_matrix[int(rel_pos_map):]
    return right_BND_map

def assemple_BND_maps(vector_repr_L, vector_repr_R, BND_rel_pos_map, matrix_len, num_diags = 1):
    
    '''
    This applies to BND predcitions.

    Get predicted matrix from Akita predictions. 
    Output is a 448x448 array with the contact frequency at each 2048 bp bin corresponding to a 917,504 bp sequence (32 bins are cropped on each end from the prediction).
 
    '''
    # Create empty matrix of NAs to place values in
    z = np.zeros(matrix_len) #,matrix_len))
    z[:BND_rel_pos_map] = vector_repr_L
    z[BND_rel_pos_map:] = vector_repr_R
    return z
