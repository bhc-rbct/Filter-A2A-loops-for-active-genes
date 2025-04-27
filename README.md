# Filter A2A Loops for Active Genes

**p2pForActiveGenes.py**

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

## Run on https://jupyter.rbct.de/

1. Click on File > New > Terminal 
2. In the Terminal enter:
    - git clone https://github.com/bhc-rbct/Filter-A2A-loops-for-active-genes.git
    - conda activate p2p

## Notes:

- Ensure that **BedTools** is installed before running the script.

# Get common p2p between two P2P files

**find_common_p2p.py**

## Parameters:

- **`-p2ps`**: Paths to the two P2P files (`p2p_plus_numInteractions_plus_associated_iaGeneList.bedpe` files) to find common p2ps. Provide two file paths separated by a space. (required)
- **`-o`**: Two file names separated by a space. The first one is used for the file with the number of interactions from the first given p2p file and the second from the second p2p file (required)

## Results:

1. Two files named according to the -o parameter, containing the common peaks and the corresponding number of interactions.
