def get_scores(
    POS, SVTYPE, SVLEN, sequences, scores, shift, revcomp, get_tracks: bool,
    cell_type=None, protocol=None, binding_factor=None, selected_tracks=None
): 
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
    if any([int(x) <= bin_size * 1088 or int(x) >= seq_length - bin_size * 1088 for x in var_rel_pos]):
        raise ValueError('Variant outside prediction window after cropping.')

    # Make predictions
    sequences = [x for x in sequences if type(x) == str]
    inputs = [np.expand_dims(one_hot_encode(s), 0).astype(np.float32) for s in sequences]
    inputs = np.vstack(inputs)
    predictions = enformer_model.predict_on_batch(inputs)

    # Process selected tracks
    matrices = []
    if selected_tracks:
        selected_indices = [human_target_index_dict[track] for track in selected_tracks if track in human_target_index_dict]
        for idx in selected_indices:
            matrices.append(predictions['human'][:, :, idx].numpy())
    else:
        # Default behavior if no tracks are selected
        if cell_type is None:
            matrices = [predictions['human'][:, :, 1186].numpy()]  # Default to CHIP:CTCF:HFF
        else:
            if cell_type == "HFF" and ((protocol is None) or (binding_factor is None)):
                protocol = "CHIP"
                cell_type = "HFF-Myc originated from foreskin fibroblast"
                binding_factor = "CTCF"
            idx = human_target_index_dict.get(f"{protocol}:{cell_type}:{binding_factor}", 1186)
            matrices = [predictions['human'][:, :, idx].numpy()]

    if revcomp:
        matrices = [np.flipud(x) for x in matrices]

    # Mask matrices
    if SVTYPE != 'BND' and abs(int(SVLEN)) > bin_size / 2:
        var_rel_pos2 = var_rel_pos.copy()
        matrices = mask_matrices(
            matrices[0], matrices[1], SVTYPE, abs(int(SVLEN)), var_rel_pos2, bin_size, offset, target_length_cropped
        )

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

