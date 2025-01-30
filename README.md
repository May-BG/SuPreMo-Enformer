# SuPreMo-

# Run the Python script with the specified arguments
python /pollard/data/projects/xzhang/tcga/SuPreMo/enformer/SuPreMo-enformer/manual_load_enformer_scripts/streamlined_SuPreMo.py \
    "$data_dir$vcf_file" \
    --dir "$data_dir/results" \
    --file "skin_melanoma_DEL_part_$1.out" \
    --get_Enformer_scores \
--selected_tracks $(awk -F'\t' '{print $NF}' /pollard/data/projects/xzhang/tcga/skin_melanoma_DEL/scripts/target_new.txt | sort | uniq)
