# Filter A2A loops for active genes
**Parameters:** 

*-uts*: unfilered list of TSS/TES regions *[required]*           

*-hf*: list of regions of active genes *[required]* (”histon file”, for L36/H3K4me3 

*-p2p:* A2A peak file from FitHiChIP *[required]* 

*-tff*: regions of interest *[required]* (usually from transcription factor)

*-g*: list of genesnames one is interested in *[not required]* 

*-o*: path to the output directory *[not required]* (default: current directory)

*-ms*: lowest score (number of contacts) to keep P2P *[not required]* (default: 1)

*-kindTS*: kind of transcription site, which should be used as region *[not required]* (default: 'TSS', choices='both', 'TSS', 'TES')

**Results:** 

- with suffixes *_P2P.bedpe*: filtered peaks
- with suffix *_perGene.txt*: number of interactions summarized per gene
- with prefix *TSS_TES_cTSS*: filtered TSS for the active (and interesting) genes. 
One can use this file as input for the next runs of the script (together with –fts instead of –ts and –hf), to skip the first two steps of the script of filtering for active Genes

**Notes:** 

- BedTools has to be installed
