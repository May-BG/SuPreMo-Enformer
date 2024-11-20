import math
import numpy as np
import pandas as pd
from typing import Any
import tensorflow.compat.v2 as tf
import tensorflow_hub as hub

from pathlib import Path
import sys
sys.path.insert(0, 'manual_load_enformer_scripts/')
from bin_utils import get_bin
from scoring_methods import common_scoring_methods
from mask_utils import mask_matrices, get_masked_BND_maps

#enformer_model = hub.load("https://kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model
enformer_model = hub.load("/pollard/home/xzhang/.cache/kagglehub/models/deepmind/enformer/tensorFlow2/enformer/1").model
repo_path = Path(__file__).parents[1]
human_targets_file = f'{repo_path}/Enformer_model/targets_human.txt'

human_targets = pd.read_csv(human_targets_file, sep='\t')
human_target_index_dict = dict(zip(human_targets["description"], human_targets["index"]))

def one_hot_encode(sequence: str,
                   alphabet: str = 'ACGT',
                   neutral_alphabet: str = 'N',
                   neutral_value: Any = 0,
                   dtype=np.float32) -> np.ndarray:
  """One-hot encode sequence."""
  def to_uint8(string):
    return np.frombuffer(string.encode('ascii'), dtype=np.uint8)
  hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
  hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
  hash_table[to_uint8(neutral_alphabet)] = neutral_value
  hash_table = hash_table.astype(dtype)
  return hash_table[to_uint8(sequence)]

def get_scores(POS, SVTYPE, SVLEN, sequences, scores, shift, revcomp, get_tracks: bool, cell_type, protocol, binding_factor, selected_tracks): 
    
    '''
    Get disruption scores, disruption tracks, and/or predicted maps from variants and the sequences generated from them.
    
    ''' 
    var_rel_pos = sequences[-1]

    seq_length = len(sequences[0])
    target_length = 896
    bin_size = 128
    offset = 1088
    target_length_cropped = target_length

    rel_pos_map = get_bin(var_rel_pos[0], bin_size, offset)
   
    # Error if variant position is too close to end of prediction window
    if any([int(x) <= bin_size*1088 or int(x) >= seq_length - bin_size*1088 for x in var_rel_pos]):
        raise ValueError('Variant outside prediction window after cropping.')
    
    
    
    # Make predictions
    
    sequences = [x for x in sequences if type(x) == str]
    
    inputs = [np.expand_dims(one_hot_encode(s), 0).astype(np.float32) for s in sequences]
    inputs = np.vstack(inputs)    
    predictions = enformer_model.predict_on_batch(inputs)
    if selected_tracks:
        selected_indices = [human_target_index_dict[track] for track in selected_tracks if track in human_target_index_dict]
        for idx in selected_indices:
            matrices.append(predictions['human'][:, :, idx].numpy())
    else:
        if protocol == "CHIP":
              
              #cell_type = "HFF-Myc originated from foreskin fibroblast"
              #binding_factor = "CTCF"
            idx = human_target_index_dict.get("{}:{}:{}".format(protocol, binding_factor,cell_type), None) ## change the order
        else:#
            idx = human_target_index_dict.get("{}:{}".format(protocol, cell_type), None) ## change the order
        if idx is None:
            idx = 1186
        matrices = predictions['human'][:, :, idx]
    if matrices.shape[0] == 2:
        matrices = [matrices[0].numpy(), matrices[1].numpy()]
    elif matrices.shape[0] == 3:
        matrices = [matrices[0].numpy(), matrices[1].numpy(), matrices[2].numpy()]
 
    if revcomp:
        matrices = [np.flipud(x) for x in matrices]
    
    
    
    # Mask matrices
    
    if SVTYPE != 'BND' and abs(int(SVLEN)) > bin_size/2:

        var_rel_pos2 = var_rel_pos.copy()
        matrices = mask_matrices(matrices[0], matrices[1], SVTYPE, abs(int(SVLEN)), var_rel_pos2, bin_size, offset, target_length_cropped)

        # If masking, the relative postion on the map depends on whether it's a duplication or deletion
        # If duplication, the relative position of variant in the ALT sequence should be used
        if SVTYPE == 'DUP':
            rel_pos_map = get_bin(var_rel_pos[1], bin_size, offset)
    
    if SVTYPE == "BND":
        matrices = get_masked_BND_maps(matrices, rel_pos_map, target_length_cropped)

    # Calculate scores 
    scores_results = {} 
    for score in scores:
        scores_results[score] = getattr(common_scoring_methods(matrices[0], matrices[1]), score)() 

        if get_tracks:
            scores_results[f"{score}_track"] = matrices

    return scores_results

