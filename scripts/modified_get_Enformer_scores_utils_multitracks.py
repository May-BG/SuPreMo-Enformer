import math
import numpy as np
import pandas as pd
from typing import Any
import tensorflow.compat.v2 as tf
import tensorflow_hub as hub

from pathlib import Path
import sys
from bin_utils import get_bin
from scoring_methods import common_scoring_methods
from mask_utils import mask_matrices, get_masked_BND_maps

# Load the Enformer model
enformer_model = hub.load("https://kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model
#enformer_model = hub.load("/pollard/home/xzhang/.cache/kagglehub/models/deepmind/enformer/tensorFlow2/enformer/1").model
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
    return hash_table[to_uint8(sequence)]

def get_scores(
    POS, SVTYPE, SVLEN, sequences, scores, shift, revcomp, get_tracks: bool, 
    cell_type, protocol, binding_factor, selected_tracks
):
    """
    Get disruption scores, including MSE and Correlation, for each track in selected_tracks.

    Args:
        POS: Variant position.
        SVTYPE: Structural variant type.
        SVLEN: Structural variant length.
        sequences: List of sequences to process.
        scores: List of scores to calculate.
        shift: Shift amount for sequences.
        revcomp: Boolean to indicate reverse complement processing.
        get_tracks: Boolean to include tracks in the output.
        cell_type: Cell type information.
        protocol: Protocol information.
        binding_factor: Binding factor information.
        selected_tracks: List of tracks to process.

    Returns:
        dict: MSE and correlation scores for each selected track.
    """
    var_rel_pos = sequences[-1]

    seq_length = len(sequences[0])
    target_length = 896
    bin_size = 128
    offset = 1088
    target_length_cropped = target_length

    rel_pos_map = get_bin(var_rel_pos[0], bin_size, offset)

    # Error if variant position is too close to the end of the prediction window
    if any([int(x) <= bin_size * 1088 or int(x) >= seq_length - bin_size * 1088 for x in var_rel_pos]):
        raise ValueError("Variant outside prediction window after cropping.")

    # Make predictions
    sequences = [x for x in sequences if isinstance(x, str)]
    inputs = [np.expand_dims(one_hot_encode(s), 0).astype(np.float32) for s in sequences]
    inputs = np.vstack(inputs)
    predictions = enformer_model.predict_on_batch(inputs)

    # Process selected tracks
    if selected_tracks:
        selected_indices = [
            human_target_index_dict[track] for track in selected_tracks if track in human_target_index_dict
        ]
    else:
        raise ValueError("No selected tracks provided.")

    scores_results = {}

    for track, idx in zip(selected_tracks, selected_indices):
        # Extract predictions for the track
        matrices = predictions["human"][:, :, idx].numpy()

        # Reverse complement if required
        if revcomp:
            matrices = np.flipud(matrices)

        # Mask matrices
        if SVTYPE != "BND" and abs(int(SVLEN)) > bin_size / 2:
            var_rel_pos2 = var_rel_pos.copy()
            matrices = mask_matrices(
                matrices[0], matrices[1], SVTYPE, abs(int(SVLEN)), var_rel_pos2, bin_size, offset, target_length_cropped
            )

        if SVTYPE == "BND":
            matrices = get_masked_BND_maps(matrices, rel_pos_map, target_length_cropped)

        # Calculate scores (MSE and Correlation for each track)
        mse = common_scoring_methods(matrices[0], matrices[1]).mse()
        corr = common_scoring_methods(matrices[0], matrices[1]).corr()

        # Add to results with track name in column
        scores_results[f"mse_{track}"] = mse
        scores_results[f"corr_{track}"] = corr

    return scores_results

# Example Usage
if __name__ == "__main__":
    # Example parameters
    POS = 123456
    SVTYPE = "DEL"
    SVLEN = 500
    sequences = ["ACGT" * 100, "TGCA" * 100, [150, 300]]
    scores_to_use = ["mse", "corr"]
    shift = 0
    revcomp = False
    get_tracks = True
    cell_type = "foreskin melanocyte male newborn"
    protocol = "CHIP"
    binding_factor = "H3K4me1"
    selected_tracks = ["DNASE:epidermal melanocyte", "CHIP:H3K4me1:foreskin melanocyte male newborn"]

    # Get scores
    results = get_scores(
        POS, SVTYPE, SVLEN, sequences, scores_to_use, shift, revcomp,
        get_tracks, cell_type, protocol, binding_factor, selected_tracks
    )

    # Print results
    for key, value in results.items():
        print(f"{key}: {value}")

