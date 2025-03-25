# Filter A2A Loops for Active Genes

## Parameters:

- **`-uts`**: Unfiltered list of TSS/TES regions (required)
- **`-hf`**: List of regions of active genes (referred to as "histone file" for L36/H3K4me3) (required)
- **`-p2p`**: A2A peak file from FitHiChIP (required)
- **`-tff`**: Regions of interest, typically from a transcription factor (required)
- **`-g`**: List of gene names of interest (optional)
- **`-o`**: Path to the output directory - parent directories and the output directory will be created, if not already there (optional, default: current directory)
- **`-ms`**: Minimum score (number of contacts) to keep in the P2P (optional, default: 1)
- **`-kindTS`**: Specifies the type of transcription site region (optional, default: 'TSS', choices: 'both', 'TSS', 'TES')

## Results:

1. **`interaction_per_enhancer_region_iag.bed`**: A file summarizing the number of interactions per enhancer region.
2. **`interaction_per_iag.txt`**: A file summarizing the number of interactions per gene.
3. **`p2p_plus_numInteractions_plus_associated_iaGeneList.bedpe`**: A filtered P2P file with the number of interactions and associated genes in a list format.
4. **`parameters_used.json`**: A JSON file summarizing the parameters that were used.
5. **`TSS_TES_for_iag_only_first_appearance_per_gene.bed`**: A file containing the TSS/TES regions (from `-uts`) but only for active (and specified) genes, keeping only the first row per gene.

## Notes:

- Ensure that **BedTools** is installed before running the script.
