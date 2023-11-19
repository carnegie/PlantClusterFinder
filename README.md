# PlantClusterFinder (version 1.3)


PlantClusterFinder (PCF) detects metabolic gene clusters in a sequenced genome. It uses a gene location file provided by the user (see below) and a PGDB created with Pathway Tools as well as further information (see below) to identify enzyme-coding genes (metabolic genes) located together on a chromosome. Initially only continous stretches of metabolic genes lying directly next to each other are allowed. This condition is relaxed by iteratively increasing the intervening (non-metabolic) gene size by one. Several criteria to select for clusters are provided. In addition to this, clusters can be prevented from forming by a section of criteria.
Details of PCF (version 1.0) can be found in PMID:
28228535.  
 
The major differences between this version (1.3) and previous versions (1.0 and 1.2) are:
 
    1) Physical breaks of the genome or sequencing gaps of unknown size are typically encoded by stretches of Ns in the genome assembly fasta file. Previously we inserted 20 hypothetical genes at each break. This however diluted the background of low quality genomes with non-enzymes, and hence the likelyhood of a cluster to be classified as top x% of enzyme dense regions was better than in a genome that had good quality. In version 1.3 we identify these breaks (but no longer insert 20 hypo genes) and prevent formation of a cluster over these gaps. 
    2) Any sequencing information that is missing is typically hard masked with Ns. Previously, any intergenic region affected by at least one N was evaluated for its length, and hypothetical genes were inserted accordingly (See Schlapfer et al, PMID:28228535). This led sometimes to unrealistic prevention of detecting gene clusters. It is unlikely that missing information about a single nucleotide would (if it would be known) lead to the finding of multiple gene models. Thus here we changed the code to insert 2 hypothetical genes only if a strech of unknown sequence is larger than nth percentile of gene sequences (set to 5). We also provide the option to NOT insert any hypothetical genes all together. Instead by default we use MaxSeqGapSize set to 100000 and MaxInterGeneDistByMedian set to 50 resulting in similar cluster predictions as in PCF version 1.0.
    3) Large gene poor intergenic regions are present in genomes. In this version we provide the option of several parameters to prevent clusters from spanning such large gene poor regions.  


Authors: Pascal Schlapfer, December 2017
         Bo Xue, December 2017


 
 Mandatory inputs: all files need full path, order does not matter.
 
    -pgdb fullpath: Give full path of PGDB flat file folder, where the flat files of the species pgdb is stored.
    
    -rmdf fullpath: TSV file of the metabolic domains of reactions. Check if file is outdated (needs to contain all reactions that are present in the PGDB that have genes annotated to).
    
    -md list:  A list of Metabolic domains that should be analyzed. For example {'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'}
    
    -psf fullpath: Protein sequence file of the genome. Can only contain Protein IDs as header e.g. >protein1-ID
    
    -gtpf fullpath: TSV file with a header describing gene, transcript, and protein identifier and then for all the identfiers the listings (gene-ID, transcript-ID, protein-ID). This is used to map information regarding transcripts (currently not used) and proteins (used for MCL clustering, potentially also used in PGDBs) to the genes (e.g. the gene location information).
    
    -glof fullpath: TSV file with gene ID, start bp, end bp, chromosome / scaffold, strand (encoded as 1 or -1). This file can be generated with a biomart or out of a gff3 file.
    
    -dnaf fullpath: Fasta file of the HARDMASKED genome nucleotide sequences.
    
    -sitf fullpath: TSV file describing the classes of enzymes that should be classified as signature or tailoring enzymes to identify clusters containing such.
    
    -gout fullpath: Path and file name to the gene-Outputfile, make sure you have access to the file and folder. This is the first results output file that will be generated.
    
    -cout fullpath: Path and file name to the cluster-Outputfile, make sure you have access to the file and folder. This is the second Results output file that will be generated.
 
 **Note**: If MCL cluster file isn't provided, the pipeline would require [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [MCL](https://micans.org/mcl/) to be installed.

 Optional Inputs:
 
    'GenesDat_GeneInfo': needs to be followed by a list of Attributes in genes.dat. Defines where gene ID information should be searched for. Default is 'Unique-ID', 'Accession-1', 'Accession-2'
    
    'MCLClusterFile': needs to be followed by either 1 or 0. With a default of 0. 1 indicates that the input file is not a protein fasta file, but a precomputed MCL clustering file (precomputing can be usefull for speed). If such a file is used, then nothing but gene IDs and tabs and new lines are allowed in this file.
                     
    'HypoGenePercentile': needs to be followed by a number. Defines the length of sequencing gap that should be populated by two hypothetical genes. If this is 10, then two hypothetical genes are introduced as soon as the sequencing gap exceeds the (100-10) lower percentile of the background size of all genes of the genome (default is 5).
    
    'MaxSeqGapSize': needs to be followed by a number. Maximal masked nucleotide region (in bp) that a cluster is allowed to bridge (default is 100000). 
    
    'MaxInterGeneDist': needs to be followed by a number. Maximal intergenic distance ( in bp) that still can be crossed by a gene cluster. (default is -1, thus inactive).

    'MaxInterGeneDistByMedian': needs to be followed by a number. Maximal intergenic distance that still can be crossed by a gene cluster, here defined by this number times the median of the gene sizes. (default is 50. Set it to -1 if it should be inactive).

    'PercentileForMaxInterGeneDist': needs to be followed by a number. Calculates the maximal intergenic distance that still can be crossed by a gene cluster. This number determines the percentile to choose from all the distances. If this is 99.9, then a cluster is allowed to bridge if the largest intergenic distance in the cluster is below 99.9% of the background size of all intergenic distances of the genome (default is -1, thus inactive).
    
    'SeqGapSizesChromBreak': needs to be followed by a list of numbers. The masked nucleotides sequences of these specific lengths will prevent cluster formation (Default is an empty list: []).
    
    'PreScreenGaps': needs to be followed by either 1 or zero. Defines if sequencing gaps should be prescreened to delete them if they are anyway too small to matter. (default is 0)
    
    'OverwriteSeqGapInfo': needs to be followed by either 1 or 0. Forces masked nucleotide analysis for finding sequencing gaps (default is 0, no enforcement).
    
    'OverwriteMCLClustering': needs to be followed by either 1 or 0. Forces (if set to 1) that the MCL clustering is redone even if the result files are already present (default is 0).
    
    'Verbose': needs to be followed by a number (0 is default, 1 gives more screen-output, 2 gives exhaustive screen-output).
    
    'EnzMinForClusters': needs to be followed by a number: Defines the minimum number of enyzmes that need to be present in a cluster (default is 3).
    
    'KbRangeLD': needs to be followed by a number. Genes that are separated by more than this size (in bp) are not considered as local duplicates. (default is 100000)
    
    'GeneRangeLD': needs to be followed by a number. Number of intervening genes that are allowed to separate two local duplicated genes. Two genes that are separated by more than this many intervening genes are not called local duplicates. (default is 10)
    
    'TopPercentClusters': needs to be followed by a number. Clusters are labeled if they are within this top X% of enzyme dense clusters compared to the background. (default is 5)
    
    'MinStepSize': needs to be followed by a number. Minimal intervening gene size that should be computed (default is 5).
    
    'MaxStepSize': needs to be followed by a number. Maximal intervening gene size that should be computed (default is 20).
    
    'Criterion': needs to be followed by a number. Criterion after which stepsize (intervening gene size) is chosen (default is 2).
        1: as published in Schlapfer et al 2017 (PMID:28228535)
            sum(enzymes in clusters) - sum(non-enzymes in clusters) > 0.
        2: stepsize 3 (because most organisms take this one)
        3: sum(enzymes in clusters) - sum(non-enzymes in clusters, without hypo genes) -2*sum(hypogenes)
        4: sum(clusters in stepsize n) > sum(clusters in stepsize n+1).
        
    -help: Display help
    
    -para: needs to be followed by a number: how many cpus can be used by the algorithm (used for MCL clustering and Sequencing gap search). The default is 1.
    
    'UnmaskedDNA': needs to be followed by a number. Defines if plant Cluster finder should continue despite a non-masked genome sequence file was used (default is 0, not continue).
    
    'PGDBIdsToMap': needs to be followed by a string. Can contain the Letters G and/or T and/or P. In that case (the pgdb is mapped to Gene (G) and/or Transcript (T) and/or Protein (P) Ids of the Gene conversion file. The default value is 'G'.
    
    'RunAfterPart': Needs to be followed by a number. Sets the parts of the algorithm that should be recomputed. (default is 0)
    
    'Tempsaves': Defines if Temporary saves should be made. (default is 0)
    
    'TempsavesOverwrite': Defines if Temporary saves should be overwritten. (Default is 0).
    
    'RemoveNonProtLocations': Defines if the protein-fasta file should be used to clean out the gene position file from sequences that do not show up on the protein fasta file, and thus do not have had a chance to be predicted to be an enzyme. (Default 0)
    
    'InsertHypos': needs to be followed by 0 or 1 (can be extended with code). Defines if the genome should be populated by hypothetical genes in sequencing gaps that are not physical breaks. (Default 0)
    
    'HypoAmount': needs to be followed by a number. Defines by how many hypothetical genes a sequencing gap that is not a physical break should be populated by artificial genes to punish such sequencing gaps from becoming clusters. (Default 0)

    'OutputType': Defines the format of the Outputfiles: 'old' is using the old format to report all clusters of all stepsizes, 'verbose' defines that top percent clusters of all stepsizes are reported. 'simple' just reports the top percent clusters of the chosen stepsize (default simple).


 Outputs:
    None (Two files are generated: vGeneOutputFile, containing all the information about genes and vClusterOutputFile, containing all the information about clusters)

 USAGE:
 
    MATLAB (paths with windows, please be aware that certain parts need to be run in linux!):
        PlantClusterFinder('-pgdb', vPGDB_FlatFileFolder, '-rmdf', vMD_reactions_File, '-md', vMD_to_annotate, '-psf', vProtein_sequence_FastaFile, '-gtpf', vGeneTranscriptProtein_mapping_File, '-glof', vGeneLocation_File, '-dnaf', vDNA_FastaFile, '-sitf', vSignatureTailorFile, '-gout', vGeneOutputFile, '-cout', vClusterOutputFile, varargin)
    
    Example: (use absolute paths!)
        PlantClusterFinder('-pgdb', '[PlantClusterFinder]\csubellipsoidea\pgdb\csubellipsoideacyc\1.0\data\', ...
        '-rmdf', '[PlantClusterFinder]\Inputs\ReactionMetabolicDomainClassification.txt', ...
        '-md', {'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'; 'Carbohydrates Metabolism'; 'Cofactors Metabolism'; 'Detoxification Metabolism'; 'Energy Metabolism'; 'Fatty Acids and Lipids Metabolism'; 'Hormones Metabolism'; 'Inorganic Nutrients Metabolism'; 'Nitrogen-Containing Compounds'; 'Nucleotides Metabolism'; 'Phenylpropanoid Derivatives'; 'Polyketides'; 'Primary-Specialized Interface Metabolism'; 'Redox Metabolism'; 'Specialized Metabolism'; 'Sugar Derivatives'; 'Terpenoids'}, ...
        '-psf', '[PlantClusterFinder]\csubellipsoidea\CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa', ...
        '-gtpf', '[PlantClusterFinder]\csubellipsoidea\gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt', ...
        '-glof', '[PlantClusterFinder]\csubellipsoidea\glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt', ...
        '-dnaf', '[PlantClusterFinder]\csubellipsoidea\CsubellipsoideaC_169_227_v2.0.hardmasked.fa', ...
        '-sitf', '[PlantClusterFinder]\Inputs\scaffold-tailoring-reactions-05082016.tab', ...
        '-gout', '[PlantClusterFinder]\csubellipsoidea\csubellipsoidea_Gene_v1_3.txt', ...
        '-cout', '[PlantClusterFinder]\csubellipsoidea\csubellipsoidea_Clust_v1_3.txt', ...
        'SeqGapSizesChromBreak', [10000], 'PGDBIdsToMap', 'GTP');
         % Note, PGDBIdsToMap is only needed here because the pgdb
         % contains protein Ids in place of gene identifiers.
        
    Shell, standalone (linux, but you can compile a windows version, see how to compile new version file):
        To use standalone, download matlab runtime from matlab website. Make sure, you download the v91. Then replace "/share/apps/MATLAB/MATLAB_Runtime/v91" int the example with your path to your runtime.
        
        USE ABSOLUTE PATHS!
        
        ./run_PlantClusterFinder.sh /share/apps/MATLAB/MATLAB_Runtime/v91 -pgdb "/PATH/TO/csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/" -rmdf "/PATH/TO/Inputs/ReactionMetabolicDomainClassification.txt" -md "{'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'; 'Carbohydrates Metabolism'; 'Cofactors Metabolism'; 'Detoxification Metabolism'; 'Energy Metabolism'; 'Fatty Acids and Lipids Metabolism'; 'Hormones Metabolism'; 'Inorganic Nutrients Metabolism'; 'Nitrogen-Containing Compounds'; 'Nucleotides Metabolism'; 'Phenylpropanoid Derivatives'; 'Polyketides'; 'Primary-Specialized Interface Metabolism'; 'Redox Metabolism'; 'Specialized Metabolism'; 'Sugar Derivatives'; 'Terpenoids'}" -psf "/PATH/TO/csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa" -gtpf "/PATH/TO/csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt" -glof "/PATH/TO/csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt" -dnaf "/PATH/TO/csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa" -sitf "/PATH/TO/Inputs/scaffold-tailoring-reactions-05082016.tab" -gout "/PATH/TO/csubellipsoidea/ARAGene1_3_memex.txt" -cout "/PATH/TO/csubellipsoidea/ARAClust1_3_memex.txt" SeqGapSizesChromBreak '[10000]' PGDBIdsToMap GTP
        % Note, PGDBIdsToMap is only needed here because the pgdb
        % contains protein Ids in place of gene identifiers.
