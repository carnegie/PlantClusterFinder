%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PlantClusterFinder 1.3
% PlantClusterFinder (PCF) detects metabolic gene clusters in a sequenced
% genome. It uses a gene location file provided by the user (see below) and
% a PGDB created with Pathway Tools as well as further information (see
% below) to identify enzyme-coding genes (metabolic genes) located together
% on a chromosome. Initially only continous stretches of metabolic genes 
% lying directly next to each other are allowed. This condition is relaxed
% by iteratively increasing the intervening (non-metabolic) gene size by
% one. Several criteria to select for clusters are provided. In addition to
% this, clusters can be prevented from forming by a section of criteria.
% Details of PCF (version 1.0) can be found in PMID: 28228535.  
% 
% The major differences between this version (1.3) and previous versions 
% (1.0 and 1.2) are:
% 
%   1) Physical breaks of the genome or sequencing gaps of unknown size are
%      typically encoded by stretches of Ns in the genome assembly fasta
%      file. Previously we inserted N hypothetical genes. This however 
%      diluted the background of low quality genomes with non-enzymes, and 
%      hence the likelyhood of a cluster to be classified as top x% of 
%      enzyme dense regions was better than in a genome that had good 
%      quality. In version 1.3 we identify these breaks and prevent 
%      formation of a cluster over these gaps. 
%   2) Any sequencing information that is missing is typically hard masked 
%      with Ns. Previously, any intergenic region affected by at least one 
%      N was evaluated for its length, and hypothetical genes were inserted 
%      accordingly (See Schlapfer et al, PMID: 28228535). This led 
%      sometimes to unrealistic prevention of detecting gene clusters. It 
%      is unlikely that missing information about a single nucleotide would 
%      (if it would be known) lead to the finding of multiple gene models. 
%      Thus here we changed the code to insert 2 hypothetical genes only if 
%      a strech of unknown sequence is larger than nth percentile of gene 
%      sequences (set to 5). We also provide the option to NOT insert any 
%      hypothetical genes all together. Instead we by default we use 
%      MaxSeqGapSize set to 100000 and MaxInterGeneDistByMedian set to 50 
%      resulting in similar clusters as in PCF version 1.0.
%   3) Large gene poor intergenic regions are present in genomes. In this 
%      version we provide the option of several parameters to prevent 
%      clusters from spanning such large gene poor regions.  
%
% Authors: Pascal Schlapfer, December 2017
%          Bo Xue, December 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Mandatory inputs: all files need full path, order does not matter.
%  
%     -pgdb fullpath: Give full path of PGDB flat file folder, where the
%                     flat files of the species pgdb is stored.
%     -rmdf fullpath: TSV file of the metabolic domains of reactions. Check
%                     if file is outdated (needs to contain all reactions
%                     that are present in the PGDB that have genes
%                     annotated to).
%     -md list: A list of Metabolic domains that should be analyzed. For
%               example {'Amines and Polyamines Metabolism'; 'Amino Acids
%               Metabolism'}
%     -psf fullpath: Protein sequence file of the genome. Can only contain
%                    Protein IDs as header e.g. >protein1-ID
%     -gtpf fullpath: TSV file with a header describing gene, transcript,
%                     and protein identifier and then for all the
%                     identfiers the listings (gene-ID, transcript-ID,
%                     protein-ID). This is used to map information
%                     regarding transcripts (currently not used) and
%                     proteins (used for MCL clustering, potentially also
%                     used in PGDBs) to the genes (e.g. the gene location
%                     information).
%     -glof fullpath: TSV file with gene ID, start bp, end bp, chromosome /
%                     scaffold, strand (encoded as 1 or -1). This file can
%                     be generated with a biomart or out of a gff3 file.
%     -dnaf fullpath: Fasta file of the HARDMASKED genome nucleotide
%                     sequences.
%     -sitf fullpath: TSV file describing the classes of enzymes that
%                     should be classified as signature or tailoring
%                     enzymes to identify clusters containing such.
%     -gout fullpath: Path and file name to the gene-Outputfile, make sure
%                     you have access to the file and folder. This is the
%                     first results output file that will be generated.
%     -cout fullpath: Path and file name to the cluster-Outputfile, make
%                     sure you have access to the file and folder. This is
%                     the second Results output file that will be
%                     generated.
%  
% Optional Inputs:
%     'GenesDat_GeneInfo': needs to be followed by a list of Attributes in
%                          genes.dat. Defines where gene ID information
%                          should be searched for. Default is 'Unique-ID',
%                          'Accession-1', 'Accession-2'
%     'MCLClusterFile': needs to be followed by either 1 or 0. With a default
%                       of 0. 1 indicates that the input file is not a
%                       protein fasta file, but a precomputed MCL
%                       clustering file (precomputing can be usefull for
%                       speed). If such a file is used, then nothing but
%                       gene IDs and tabs and new lines are allowed in this
%                       file.
%     'HypoGenePercentile': needs to be followed by a number. Defines the
%                           length of sequencing gap that should be
%                           populated by two hypothetical genes. If this is
%                           10, then two hypothetical genes are introduced
%                           as soon as the sequencing gap exceeds the
%                           (100-10) lower percentile of the background
%                           size of all genes of the genome (default is 5).
%     'MaxSeqGapSize': needs to be followed by a number. Maximal masked
%                      nucleotide region (in bp) that a cluster is allowed
%                      to bridge (default is 100000). 
%     'MaxInterGeneDist': needs to be followed by a number. Maximal
%                         intergenic distance ( in bp) that still can be
%                         crossed by a gene cluster. (default is -1, thus
%                         inactive).
%     'MaxInterGeneDistByMedian': needs to be followed by a number. Maximal
%                                 intergenic distance that still can be
%                                 crossed by a gene cluster, here defined
%                                 by this number times the median of the
%                                 gene sizes. (default is 50. Set it to -1
%                                 if it should be inactive).
%     'PercentileForMaxInterGeneDist': needs to be followed by a number.
%                                      Calculates the maximal intergenic
%                                      distance that still can be crossed
%                                      by a gene cluster. This number
%                                      determines the percentile to choose
%                                      from all the distances. If this is
%                                      99.9, then a cluster is allowed to
%                                      bridge if the largest intergenic
%                                      distance in the cluster is below
%                                      99.9% of the background size of all
%                                      intergenic distances of the genome
%                                      (default is -1, thus inactive).
%     'SeqGapSizesChromBreak': needs to be followed by a list of numbers.
%                              The masked nucleotides sequences of these
%                              specific lengths will prevent cluster
%                              formation (Default is an empty list: []).
%     'PreScreenGaps': needs to be followed by either 1 or zero. Defines if
%                      sequencing gaps should be prescreened to delete them
%                      if they are anyway too small to matter. (default is
%                      0)
%     'OverwriteSeqGapInfo': needs to be followed by either 1 or 0. Forces
%                            masked nucleotide analysis for finding
%                            sequencing gaps (default is 0, no
%                            enforcement).
%     'OverwriteMCLClustering': needs to be followed by either 1 or 0.
%                               Forces (if set to 1) that the MCL
%                               clustering is redone even if the result
%                               files are already present (default is 0).
%     'Verbose': needs to be followed by a number (0 is default, 1 gives
%                more screen-output, 2 gives exhaustive screen-output).
%     'EnzMinForClusters': needs to be followed by a number: Defines the
%                          minimum number of enyzmes that need to be
%                          present in a cluster (default is 3).
%     'KbRangeLD': needs to be followed by a number. Genes that are
%                  separated by more than this size (in bp) are not
%                  considered as local duplicates. (default is 100000)
%     'GeneRangeLD': needs to be followed by a number. Number of
%                    intervening genes that are allowed to separate two
%                    local duplicated genes. Two genes that are separated
%                    by more than this many intervening genes are not
%                    called local duplicates. (default is 10)
%     'TopPercentClusters': needs to be followed by a number. Clusters are
%                           labeled if they are within this top X% of
%                           enzyme dense clusters compared to the
%                           background. (default is 5)
%     'MinStepSize': needs to be followed by a number. Minimal intervening
%                    gene size that should be computed (default is 5).
%     'MaxStepSize': needs to be followed by a number. Maximal intervening
%                    gene size that should be computed (default is 20).
%     'Criterion': needs to be followed by a number. Criterion after which
%                  stepsize (intervening gene size) is chosen (default is
%                  2).
%         1: as published in Schlapfer et al 2017 (PMID:28228535)
%             sum(enzymes in clusters) - sum(non-enzymes in clusters) > 0.
%         2: stepsize 3 (because most organisms take this one)
%         3: sum(enzymes in clusters) - sum(non-enzymes in clusters, without hypo genes) -2*sum(hypogenes)
%         4: sum(clusters in stepsize n) > sum(clusters in stepsize n+1).
%     -help: Display help
%     -para: needs to be followed by a number: how many cpus can be used by
%            the algorithm (used for MCL clustering and Sequencing gap
%            search). The default is 1.
%     'UnmaskedDNA': needs to be followed by a number. Defines if plant
%                    Cluster finder should continue despite a non-masked
%                    genome sequence file was used (default is 0, not
%                    continue).
%     'PGDBIdsToMap': needs to be followed by a string. Can contain the
%                     Letters G and/or T and/or P. In that case (the pgdb
%                     is mapped to Gene (G) and/or Transcript (T) and/or
%                     Protein (P) Ids of the Gene conversion file. The
%                     default value is 'G'.
%     'RunAfterPart': Needs to be followed by a number. Sets the parts of
%                     the algorithm that should be recomputed. (default is
%                     0)
%     'Tempsaves': Defines if Temporary saves should be made. (default is
%                  0)
%     'TempsavesOverwrite': Defines if Temporary saves should be
%                           overwritten. (Default is 0).
%     'RemoveNonProtLocations': Defines if the protein-fasta file should be
%                               used to clean out the gene position file
%                               from sequences that do not show up on the
%                               protein fasta file, and thus do not have
%                               had a chance to be predicted to be an
%                               enzyme. (Default 0)
%     'InsertHypos': needs to be followed by 0 or 1 (can be extended with
%                    code). Defines if the genome should be populated by
%                    hypothetical genes in sequencing gaps that are not
%                    physical breaks. (Default 0)
%     'HypoAmount': needs to be followed by a number. Defines by how many
%                   hypothetical genes a sequencing gap that is not a
%                   physical break should be populated by artificial genes
%                   to punish such sequencing gaps from becoming clusters.
%                   (Default 0)
%     'OutputType': Defines the format of the Outputfiles: 'old' is using
%                   the old format to report all clusters of all stepsizes,
%                   'verbose' defines that top percent clusters of all
%                   stepsizes are reported. 'simple' just reports the top
%                   percent clusters of the chosen stepsize (default
%                   simple).
% 
%  Outputs:
%     None (Two files are generated: vGeneOutputFile, containing all the
%          information about genes and vClusterOutputFile, containing all
%          the information about clusters)
% 
%  USAGE:
%  
%     MATLAB (paths with windows, please be aware that certain parts need
%            to be run in linux!):
%         PlantClusterFinder('-pgdb', vPGDB_FlatFileFolder, '-rmdf', ...
%                            vMD_reactions_File, '-md', ...
%                            vMD_to_annotate, '-psf', ...
%                            vProtein_sequence_FastaFile, '-gtpf', ...
%                            vGeneTranscriptProtein_mapping_File, ...
%                            '-glof', vGeneLocation_File, '-dnaf', ...
%                            vDNA_FastaFile, '-sitf', ...
%                            vSignatureTailorFile, '-gout', ...
%                            vGeneOutputFile, '-cout', ...
%                            vClusterOutputFile, varargin)
%     
%     Example: (use absolute paths!)
%         PlantClusterFinder('-pgdb', '[PlantClusterFinder]\csubellipsoidea\pgdb\csubellipsoideacyc\1.0\data\', ...
%         '-rmdf', '[PlantClusterFinder]\Inputs\ReactionMetabolicDomainClassification.txt', ...
%         '-md', {'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'; 'Carbohydrates Metabolism'; 'Cofactors Metabolism'; 'Detoxification Metabolism'; 'Energy Metabolism'; 'Fatty Acids and Lipids Metabolism'; 'Hormones Metabolism'; 'Inorganic Nutrients Metabolism'; 'Nitrogen-Containing Compounds'; 'Nucleotides Metabolism'; 'Phenylpropanoid Derivatives'; 'Polyketides'; 'Primary-Specialized Interface Metabolism'; 'Redox Metabolism'; 'Specialized Metabolism'; 'Sugar Derivatives'; 'Terpenoids'}, ...
%         '-psf', '[PlantClusterFinder]\csubellipsoidea\CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa', ...
%         '-gtpf', '[PlantClusterFinder]\csubellipsoidea\gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt', ...
%         '-glof', '[PlantClusterFinder]\csubellipsoidea\glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt', ...
%         '-dnaf', '[PlantClusterFinder]\csubellipsoidea\CsubellipsoideaC_169_227_v2.0.hardmasked.fa', ...
%         '-sitf', '[PlantClusterFinder]\Inputs\scaffold-tailoring-reactions-05082016.tab', ...
%         '-gout', '[PlantClusterFinder]\csubellipsoidea\csubellipsoidea_Gene_v1_3.txt', ...
%         '-cout', '[PlantClusterFinder]\csubellipsoidea\csubellipsoidea_Clust_v1_3.txt', ...
%         'SeqGapSizesChromBreak', [10000], 'PGDBIdsToMap', 'GTP');
%         % Note, PGDBIdsToMap is only needed here because the pgdb
%         % contains protein Ids in place of gene identifiers.
%         
%     Shell, standalone (linux, but you can compile a windows version, see
%                       how to compile new version file):
%         To use standalone, download matlab runtime from matlab website.
%         Make sure, you download the v91. Then replace
%         "/share/apps/MATLAB/MATLAB_Runtime/v91" int the example with your
%         path to your runtime.
%         
%         USE ABSOLUTE PATHS!
%
%         ./run_PlantClusterFinder.sh /share/apps/MATLAB/MATLAB_Runtime/v91 -pgdb "./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/" -rmdf "./Inputs/ReactionMetabolicDomainClassification.txt" -md "{'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'; 'Carbohydrates Metabolism'; 'Cofactors Metabolism'; 'Detoxification Metabolism'; 'Energy Metabolism'; 'Fatty Acids and Lipids Metabolism'; 'Hormones Metabolism'; 'Inorganic Nutrients Metabolism'; 'Nitrogen-Containing Compounds'; 'Nucleotides Metabolism'; 'Phenylpropanoid Derivatives'; 'Polyketides'; 'Primary-Specialized Interface Metabolism'; 'Redox Metabolism'; 'Specialized Metabolism'; 'Sugar Derivatives'; 'Terpenoids'}" -psf "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa" -gtpf "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt" -glof "./csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt" -dnaf "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa" -sitf "./Inputs/scaffold-tailoring-reactions-05082016.tab" -gout "./csubellipsoidea/csubellipsoidea_Gene1_3.txt" -cout "./csubellipsoidea/csubellipsoidea_Clust1_3.txt" SeqGapSizesChromBreak '[10000]' PGDBIdsToMap GTP
%         % Note, PGDBIdsToMap is only needed here because the pgdb
%         % contains protein Ids in place of gene identifiers.
% 
function PlantClusterFinder(varargin)
    %% Display help if no arguments
    if nargin == 0
        f_display_help()
        return;
    end
    
    %% Display parameters:
    if ispc()
        vc = 1;
        for vi = 1:size(varargin,2)
            if vc == 1
                fprintf('Parameter %i: %s\n', (vi-1)/2+1, varargin{1,vi})
                vc = vc + 1;
            else
                if ischar(varargin{1,vi})
                    fprintf('Value: %s\n', varargin{1,vi});
                    vc = 1;
                elseif isnumeric(varargin{1,vi})
                    fprintf('Value: %i\n', varargin{1,vi});
                    vc = 1;
                else
                    fprintf('Value: can not be displayed\n');
                    vc = 1;
                end
            end
        end
    else
        vc = 1;
        for vi = 1:size(varargin,2)
            if vc == 1
                fprintf('Parameter %i: %s\n', (vi-1)/2+1, varargin{1,vi})
                vc = vc + 1;
            else
                fprintf('Value: %s\n', varargin{1,vi})
                vc = 1;
            end
        end
    end
    
    %% Initializing parameters
    fprintf('Initializing parameters and running PlantClusterFinder\n');
    vVerbose = 0;
    vCurrArgin = 0;
    vGenesDat_GeneInfo = {'Unique-ID', 'Accession-1', 'Accession-2'};
    vMCLClusterFile = 0;
    vHypoGenePercentile = 5;
    vMaxSeqGapSize = 100000;
    vMaxInterGeneDist = -1;
    vMaxInterGeneDistByMedian = 50;
    vPercentileForMaxInterGeneDist = -1;
    vSeqGapSizesChromBreak = [];
    vPreScreenGaps = 0;
    vMinStepSize = 5;
    vMaxStepSize = 20;
    vCriterion = 2;
    vTopPercentClusters = 5;
    vOverwriteSeqGapInfo = 0;
    vOverwriteMCLClustering = 0;
    vKbRangeLD = 100000;
    vGeneRangeLD = 10;
    vEnzMinForClusters = 3;
    vPara = 1;
    vUnmaskedDNA = 0;
    vPGDBIdsToMap = 'G';
    vPartFinished = -1;
    vRunAfterPart = 0;
    vTempsaves = 0;
    vTempsavesOverwrite = 0;
    vRemoveNonProtLocations = 0;
    vInsertHypos = 0;
    vOutputType = 'simple';
    if vPartFinished < 0
        if nargin > vCurrArgin
            vi = 1;
            while vi <= nargin - vCurrArgin
                if ischar(varargin{1,vi})

                    if strcmp(varargin{1,vi},'-pgdb')
                        if ischar(varargin{1,vi + 1})
                            vPGDB_FlatFileFolder = varargin{1,vi + 1};
                            if ~ispc()
                                if isempty(regexp(vPGDB_FlatFileFolder,'/$','once'))
                                    vPGDB_FlatFileFolder = [vPGDB_FlatFileFolder '/']; %#ok<*AGROW>
                                end
                            else
                                if isempty(regexp(vPGDB_FlatFileFolder,'\\$','once'))
                                    vPGDB_FlatFileFolder = [vPGDB_FlatFileFolder '\'];
                                end
                            end
                        else
                            error('PGDB Flatfilefolder argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-rmdf')
                        if ischar(varargin{1,vi + 1})
                            vMD_reactions_File = varargin{1,vi + 1};
                        else
                            error('MD file for reactions argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-md')
                        if ischar(varargin{1,vi + 1})
                            eval(['vMD_to_annotate =' ' ' varargin{1,vi + 1} ';']);
                        elseif iscell(varargin{1,vi + 1})
                            vMD_to_annotate = varargin{1,vi + 1};
                        else
                            error('Metabolic domains to annotate list argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-psf')
                        if ischar(varargin{1,vi + 1})
                            vProtein_sequence_FastaFile = varargin{1,vi + 1};
                        else
                            error('Protein sequence file argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-gtpf')
                        if ischar(varargin{1,vi + 1})
                            vGeneTranscriptProtein_mapping_File = varargin{1,vi + 1};
                        else
                            error('Gene transcript protein id conversion file argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-glof')
                        if ischar(varargin{1,vi + 1})
                            vGeneLocation_File = varargin{1,vi + 1};
                        else
                            error('Gene location (postion) file argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-dnaf')
                        if ischar(varargin{1,vi + 1})
                            vDNA_FastaFile = varargin{1,vi + 1};
                        else
                            error('Masked dna fasta file argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-sitf')
                        if ischar(varargin{1,vi + 1})
                            vSignatureTailorFile = varargin{1,vi + 1};
                        else
                            error('Signature and tailoring reacitons file argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-gout')
                        if ischar(varargin{1,vi + 1})
                            vGeneOutputFile = varargin{1,vi + 1};
                        else
                            error('Outputfile for gene summary argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-cout')
                        if ischar(varargin{1,vi + 1})
                            vClusterOutputFile = varargin{1,vi + 1};
                        else
                            error('Outputfile for cluster summary argument is not of type string\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'GenesDat_GeneInfo')
                        if ischar(varargin{1,vi + 1})
                            eval(['vGenesDat_GeneInfo = ' varargin{1,vi + 1} ';']);
                        elseif iscell(varargin{1,vi + 1})
                            vGenesDat_GeneInfo = varargin{1,vi + 1};
                        else
                            error('GenesDat_GeneInfo does not follow format {''Attribute1'', ''Atribute2''}\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'MCLClusterFile')
                        if ischar(varargin{1,vi + 1})
                            vMCLClusterFile = str2num(varargin{1,vi + 1}); %#ok<*ST2NM>
                        elseif isnumeric(varargin{1,vi + 1})
                            vMCLClusterFile = varargin{1,vi + 1};
                        else
                            error('MCLClusterFile is not followed by a number (either 1 or 0)\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'HypoGenePercentile')
                        if isnumeric(varargin{1,vi + 1})
                            vHypoGenePercentile = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vHypoGenePercentile = str2num(varargin{1,vi + 1});
                        else
                            error('HypoGenePercentile is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'MaxSeqGapSize')
                        if isnumeric(varargin{1,vi + 1})
                            vMaxSeqGapSize = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vMaxSeqGapSize = str2num(varargin{1,vi + 1});
                        else
                            error('MaxSeqGapSize is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'MaxInterGeneDist')
                        if isnumeric(varargin{1,vi + 1})
                            vMaxInterGeneDist = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vMaxInterGeneDist = str2num(varargin{1,vi + 1});
                        else
                            error('MaxInterGeneDist is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'MaxInterGeneDistByMedian')
                        if isnumeric(varargin{1,vi + 1})
                            vMaxInterGeneDistByMedian = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vMaxInterGeneDistByMedian = str2num(varargin{1,vi + 1});
                        else
                            error('MaxInterGeneDist is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'PercentileForMaxInterGeneDist')
                        if isnumeric(varargin{1,vi + 1})
                            vPercentileForMaxInterGeneDist = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vPercentileForMaxInterGeneDist = str2num(varargin{1,vi + 1});
                        else
                            error('PercentileForMaxInterGeneDist is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'SeqGapSizesChromBreak')
                        if isnumeric(varargin{1,vi + 1})
                            vSeqGapSizesChromBreak = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            if strcmp(varargin{1,vi + 1},'')
                                eval('vSeqGapSizesChromBreak = [];');
                            else
                                eval(['vSeqGapSizesChromBreak = ' varargin{1,vi + 1} ';']);
                            end
                        else
                            error('SeqGapSizesChromBreak does not follow format [Number1, Number2]\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'PreScreenGaps')
                        if isnumeric(varargin{1,vi + 1})
                            vPreScreenGaps = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vPreScreenGaps = str2num(varargin{1,vi + 1});
                        else
                            error('PreScreenGaps is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'OverwriteSeqGapInfo')
                        if isnumeric(varargin{1,vi + 1})
                            vOverwriteSeqGapInfo = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vOverwriteSeqGapInfo = str2num(varargin{1,vi + 1});
                        else
                            error('OverwriteSeqGapInfo is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'OverwriteMCLClustering')
                        if isnumeric(varargin{1,vi + 1})
                            vOverwriteMCLClustering = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vOverwriteMCLClustering = str2num(varargin{1,vi + 1});
                        else
                            error('OverwriteMCLClustering is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'Verbose')
                        if isnumeric(varargin{1,vi + 1})
                            vVerbose = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vVerbose = str2num(varargin{1,vi + 1});
                        else
                            error('Verbose is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'EnzMinForClusters')
                        if isnumeric(varargin{1,vi + 1})
                            vEnzMinForClusters = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vEnzMinForClusters = str2num(varargin{1,vi + 1});
                        else
                            error('EnzMinForClusters is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'KbRangeLD')
                        if isnumeric(varargin{1,vi + 1})
                            vKbRangeLD = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vKbRangeLD = str2num(varargin{1,vi + 1});
                        else
                            error('KbRangeLD is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'GeneRangeLD')
                        if isnumeric(varargin{1,vi + 1})
                            vGeneRangeLD = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vGeneRangeLD = str2num(varargin{1,vi + 1});
                        else
                            error('GeneRangeLD is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'TopPercentClusters')
                        if isnumeric(varargin{1,vi + 1})
                            vTopPercentClusters = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vTopPercentClusters = str2num(varargin{1,vi + 1});
                        else
                            error('TopPercentClusters is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'MinStepSize')
                        if isnumeric(varargin{1,vi + 1})
                            vMinStepSize = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vMinStepSize = str2num(varargin{1,vi + 1});
                        else
                            error('MinStepSize is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'MaxStepSize')
                        if isnumeric(varargin{1,vi + 1})
                            vMaxStepSize = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vMaxStepSize = str2num(varargin{1,vi + 1});
                        else
                            error('MaxStepSize is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'Criterion')
                        if isnumeric(varargin{1,vi + 1})
                            vCriterion = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vCriterion = str2num(varargin{1,vi + 1});
                        else
                            error('Criterion is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'-para')
                        if isnumeric(varargin{1,vi + 1})
                            vPara = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vPara = str2num(varargin{1,vi + 1});
                        else
                            error('-para is not followed by a string that can be converted into a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'help') || strcmp(varargin{1,vi},'h') || strcmp(varargin{1,vi},'man') || ...
                            strcmp(varargin{1,vi},'info') || strcmp(varargin{1,vi},'?') || strcmp(varargin{1,vi},'-help') ...
                            || strcmp(varargin{1,vi},'-h') || strcmp(varargin{1,vi},'-info') || strcmp(varargin{1,vi},'-man') || strcmp(varargin{1,vi},'-?')
                        f_display_help()
                        return;
                    elseif strcmp(varargin{1,vi},'UnmaskedDNA')
                        if isnumeric(varargin{1,vi + 1})
                            vUnmaskedDNA = varargin{1,vi + 1};
                        elseif ischar(varargin{1,vi + 1})
                            vUnmaskedDNA = str2num(varargin{1,vi + 1});
                        else
                            error('UnmaskedDNA is not followed by a number\n');
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'PGDBIdsToMap')
                        if ~ischar(varargin{1,vi + 1})
                            error('UnmaskedDNA is not followed by a string\n');
                        else
                            vPGDBIdsToMap = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'RunAfterPart')
                        if ~ischar(varargin{1,vi + 1}) && ~isnumeric(varargin{1,vi + 1})
                            error('RunAfterPart is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vRunAfterPart = str2num(varargin{1,vi + 1});
                        elseif isnumeric(varargin{1,vi + 1})
                            vRunAfterPart = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'Tempsaves')
                        if ~ischar(varargin{1,vi + 1}) && ~isnumeric(varargin{1,vi + 1})
                            error('Tempsaves is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vTempsaves = str2num(varargin{1,vi + 1});
                        elseif isnumeric(varargin{1,vi + 1})
                            vTempsaves = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'TempsavesOverwrite')
                        if ~ischar(varargin{1,vi + 1}) && ~isnumeric(varargin{1,vi + 1})
                            error('TempsavesOverwrite is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vTempsavesOverwrite = str2num(varargin{1,vi + 1});
                        elseif isnumeric(varargin{1,vi + 1})
                            vTempsavesOverwrite = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'RemoveNonProtLocations')
                        if ~ischar(varargin{1,vi + 1}) && ~isnumeric(varargin{1,vi + 1})
                            error('RemoveNonProtLocations is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vRemoveNonProtLocations = str2num(varargin{1,vi + 1});
                        elseif isnumeric(varargin{1,vi + 1})
                            vRemoveNonProtLocations = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'InsertHypos')
                        if ~ischar(varargin{1,vi + 1}) && ~isnumeric(varargin{1,vi + 1})
                            error('InsertHypos is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vInsertHypos = str2num(varargin{1,vi + 1});
                        elseif isnumeric(varargin{1,vi + 1})
                            vInsertHypos = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'HypoAmount')
                        if ~ischar(varargin{1,vi + 1}) && ~isnumeric(varargin{1,vi + 1})
                            error('HypoAmount is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vHypoAmount = str2num(varargin{1,vi + 1});
                        elseif isnumeric(varargin{1,vi + 1})
                            vHypoAmount = varargin{1,vi + 1};
                        end
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'OutputType')
                        if ~ischar(varargin{1,vi + 1})
                            error('OutputType is not followed by a string\n');
                        elseif ischar(varargin{1,vi + 1})
                            vOutputType = varargin{1,vi + 1};
                            if strcmp(vOutputType,'simple') || strcmp(vOutputType,'Simple')
                                vOutputType = 'simple';
                            elseif strcmp(vOutputType,'verbose') || strcmp(vOutputType,'Verbose')
                                vOutputType = 'verbose';
                            elseif strcmp(vOutputType,'old') || strcmp(vOutputType,'Old')
                                vOutputType = 'old';
                            else
                                error('OutputType is not one of to following: simple, verbose\n');
                            end
                        end
                        vi = vi + 1;
                    else
                        if ischar(varargin{1,vi})
                            error('Argument %s not supported\n',varargin{1,vi});
                        elseif vi > 1
                            error('Argument %i (after %s) not supported', vi/2,varargin{1,vi-2});
                        else
                            error('Argument %i not supported');
                        end
                    end
                end
                vi = vi + 1;
            end
        end

        vMissingInfo = 0;
        if ~exist('vPGDB_FlatFileFolder','var')
            fprintf('Please give the pgdb flat file folder with -pgdb fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vMD_reactions_File','var')
            fprintf('Please give the metabolic domains for reactions file with -rmdf fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vMD_to_annotate','var')
            fprintf('Please give a list of metabolic domains to annotate by -md ''{''metabolic domain 1'';''metabolic domain 2''}');
            vMissingInfo = 1;
        end
        if ~exist('vProtein_sequence_FastaFile','var')
            fprintf('Please give the protein fasta file with -psf fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vGeneTranscriptProtein_mapping_File','var')
            fprintf('Please give gene transcript protein id conversion file with -gtpf fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vGeneLocation_File','var')
            fprintf('Please give gene location (position) file with -glof fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vDNA_FastaFile','var')
            fprintf('Please give the masked dna fasta file with -dnaf fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vSignatureTailorFile','var')
            fprintf('Please give the signature tailoring reactions file with -sitf fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vGeneOutputFile','var')
            fprintf('Please give the path to where the gene summary should be stored by -gout fullpath');
            vMissingInfo = 1;
        end
        if ~exist('vClusterOutputFile','var')
            fprintf('Please give the path to where the cluster summary should be stored by -gout fullpath');
            vMissingInfo = 1;
        end

        if vMissingInfo == 1
            error('Not all mandatory arguments were accepted, run aborted.\n');
        end

        %% Check files if we can open the files:
        fprintf('Check files for reading and writing\n');
        f_check_files_for_read_write(vPGDB_FlatFileFolder, ...
                                     vMD_reactions_File, ...
                                     vProtein_sequence_FastaFile, ...
                                     vGeneTranscriptProtein_mapping_File, ...
                                     vGeneLocation_File, vDNA_FastaFile, ...
                                     vSignatureTailorFile, vGeneOutputFile, ...
                                     vClusterOutputFile)
        vPartFinished = 0;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    if exist([vGeneLocation_File '_tempsave.mat'],'file') && vRunAfterPart ~= 0 
        vPartFinished_old = vPartFinished;
        load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
        if vPartFinished >= vRunAfterPart
            load([vGeneLocation_File '_tempsave.mat']); %#ok<*LOAD>
            vPartFinished = vRunAfterPart;
        elseif vTempsavesOverwrite == 1
            vPartFinished = vPartFinished_old;
            clear vPartFinished_old
        else
            load([vGeneLocation_File '_tempsave.mat']);
            clear vPartFinished_old
        end
    end
    
    %% Read in Gene-Transcript-Protein File
    if vPartFinished < 1
        if vVerbose >= 1
            fprintf(['Read in file describing gene-ID, transcript-ID' ...
                ', and protein-ID conversion\n']);
        end
        [vMapG, vMapT, vMapP, vMapG_MapT, vMapG_MapP] = ...
            f_read_in_ConversionFile(vGeneTranscriptProtein_mapping_File, ...
            vVerbose);
        vPartFinished = 1;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Map PGDB Gene-ids to reactions
    if vPartFinished < 2
        if vVerbose >= 1
            fprintf('Read in pgdb flat files\n');
        end
        [vGenes_all, ~, vReactions_all, vReactions_all_Att, vG2R] = ...
            f_map_pgdb_reactions_to_genes(vPGDB_FlatFileFolder, ...
            vVerbose);
        
        vPartFinished = 2;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Map PGDB reactions to metabolic domains
    if vPartFinished < 3
        if vVerbose >= 1
            fprintf('Read in and map metabolic domains\n');
        end
        [vRxnMetaDom, vRxnMetaDom_Att, vR_MD] = ...
            f_get_metabolic_domains(vMD_reactions_File, vMD_to_annotate);
        %Match PGDBs to MDs
        [vRpgdb_Rmd] = f_matchPGDBs_to_MDs(vG2R, vR_MD, ...
            vReactions_all, vReactions_all_Att, vRxnMetaDom, vRxnMetaDom_Att);
        vPartFinished = 3;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Lable Reactions with signature and tailoring
    if vPartFinished < 4
        if vVerbose >= 1
            fprintf(['Label reactions with signature and tailor' ...
                'ing information\n']);
        end
        [vSTClasses, vSigOrTail, vR_STClasses, vR_SigOrTail] = ...
            f_lable_Reactions_withSigandTail(vReactions_all, ...
            vReactions_all_Att, vSignatureTailorFile);
        vPartFinished = 4;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Map PGDB genes to Maping genes
    if vPartFinished < 5
        if vVerbose >= 1
            fprintf('Map pgdb genes to conversion file gene-IDs\n');
        end
        vGenesDat_GeneInfo = upper(regexprep(vGenesDat_GeneInfo,'-','_'));
        vMapG_Gpgdb = f_map_PGDBs_to_conversion_file(vMapG, vMapT, vMapP, vMapG_MapT, vMapG_MapP, vGenes_all, ...
            vGenesDat_GeneInfo, vPGDBIdsToMap, vVerbose);
        vPartFinished = 5;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Create MCL clustering mapped to Maping genes
    if vPartFinished < 6
        if vMCLClusterFile == 0
            if vVerbose >= 1
                fprintf(['Create MCL clustering based on protein file a' ...
                    'nd map cluster IDs to conversion gene-IDs\n']);
            end
            [vG_MCL, vProtein_IDs] = f_get_MCL_clustering(vProtein_sequence_FastaFile, ...
                vMapG, vMapP, vMapG_MapP, vPara, vOverwriteMCLClustering, vRemoveNonProtLocations, vVerbose);
        else
            if vVerbose >= 1
                fprintf(['Read in MCL clustering and map cluster ID' ...
                    's to conversion gene-IDs\n']);
            end
            [vG_MCL, vProtein_IDs] = ...
                f_read_MCL_clustering_File(vProtein_sequence_FastaFile, ...
                vMapG, vMapP, vMapG_MapP, vRemoveNonProtLocations, vVerbose);
        end
        vPartFinished = 6;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Read in Gene Position File
    if vPartFinished < 7
        if vVerbose >= 1
            fprintf('Read in gene position file\n');
        end
        [vGeneLocation, vGeneLocation_Att] = ...
            f_extract_results_with_header(vGeneLocation_File);
        if vVerbose >= 2
            fprintf('Gene Location File Attributes:\n')
            for vi = 1:size(vGeneLocation_Att,2)
                if ~iscell(vGeneLocation_Att(1,vi))
                    fprintf('Attribute %i is broken.\n', vi);
                else
                    fprintf('Attribute %i: %s\n', vi, char(vGeneLocation_Att(1,vi)));
                end
            end
        end
        %Remove protein or transcript information
        if vVerbose >= 2
            fprintf('Remove transcript information\n');
        end
        [vGeneLocation] = ...
            f_remove_protein_transcript_info_from_genepositionfile(...
            vGeneLocation, vGeneLocation_Att, vMapG, vMapT, vMapP, vMapG_MapT, ...
            vMapG_MapP, vVerbose);
        %Consolidate gene information (one entry per gene)
        if vVerbose >= 2
            fprintf('Consolidate gene information\n');
        end
        [vUniqueChrom, vGeneLocation2] = ...
            f_consolidate_gene_info(vGeneLocation, vGeneLocation_Att, vVerbose);
        %Sort genes according to their physical positon on the chromosomes
        if vVerbose >= 2
            fprintf('Sort genes on the chromosomes\n');
        end
        [vGeneLocation2] = f_sort_genes(vGeneLocation2, ...
            vGeneLocation_Att, vUniqueChrom, vVerbose);
        vPartFinished = 7;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Map Biomartfile to Gene Map
    if vPartFinished < 8
        if vVerbose >= 1
            fprintf('Map gene location to conversion gene-IDs\n');
        end
        [vBiomG_G] = f_map_GeneLocation_to_conversion_file(vGeneLocation2, ...
            vGeneLocation_Att, vMapG, vVerbose);
        
        %If chosen, remove gene position information about genes that are
        %not present in the protein fasta file.
        if vRemoveNonProtLocations == 1
            if vVerbose >= 2
                fprintf('Remove none protein coding genes in gene locations.\n');
            end
            [vGeneLocation2, vBiomG_G] = f_Remove_non_protein_fasta_gene_locations(vGeneLocation2, ...
                vGeneLocation_Att, vMapP, vMapG_MapP, vBiomG_G, vProtein_IDs);
        end
        
        vPartFinished = 8;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Get places on genome that do not have genes encoded:
    if vPartFinished < 9
        if vVerbose >= 1
            fprintf('Find intergenic regions\n');
        end
        [vInterChrom] = ...
            f_find_intergenic_regions(vGeneLocation2, vGeneLocation_Att, ...
            vUniqueChrom, vVerbose);
        vPartFinished = 9;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Calculate percentile of gene size that should yield hypo genes
    if vPartFinished < 10
        if vVerbose >= 1
            fprintf(['Calculate size of sequence gap that should be pop' ...
                'ulated by hypothetical genes\n']);
        end
        vGeneLengthHypoMin = ...
            ceil(prctile(str2num(char(vGeneLocation2.(char(...
            vGeneLocation_Att(1,3)))(:,1)))-str2num(char(vGeneLocation2.(...
            char(vGeneLocation_Att(1,2)))(:,1))),vHypoGenePercentile));
        vPartFinished = 10;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Get sequencing gaps
    if vPartFinished < 11    
        if vVerbose >= 1
            fprintf('Identify sequencing gaps (bases encoded by N)\n');
        end
        if vGeneLengthHypoMin > 0 && vPreScreenGaps == 1
            vSeqGaps = f_annotate_Sequencing_Gaps(vDNA_FastaFile, vInterChrom,...
                vUniqueChrom, vOverwriteSeqGapInfo, vPara, vUnmaskedDNA, vVerbose, 'GeneLengthHypoMin', vGeneLengthHypoMin, 'SeqGapSizesChromBreak',vSeqGapSizesChromBreak);
        else
            vSeqGaps = f_annotate_Sequencing_Gaps(vDNA_FastaFile, vInterChrom,...
                vUniqueChrom, vOverwriteSeqGapInfo, vPara, vUnmaskedDNA, vVerbose);
        end
        vPartFinished = 11;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Introduce Hardbreaks (split chromosomes)
    if vPartFinished < 12
        if vVerbose >= 1
            if vMaxInterGeneDist ~= -1 || vMaxSeqGapSize ~= -1 || size(vSeqGapSizesChromBreak,2) ~= 0 || vPercentileForMaxInterGeneDist ~= -1 ||  vMaxInterGeneDistByMedian~= -1
                fprintf('Prevent clusters:\n');
                if vMaxSeqGapSize ~= -1
                    fprintf('Prevent clusters spanning sequencing gaps longer or equal than %ibp\n', vMaxSeqGapSize);
                end
                if size(vSeqGapSizesChromBreak,2) ~= 0
                    vNumbers = num2str(vSeqGapSizesChromBreak(1,1));
                    for vi = 2:size(vSeqGapSizesChromBreak,2)
                        vNumbers = [vNumbers ', ' ...
                            num2str(vSeqGapSizesChromBreak(1,vi))];
                    end
                    if vVerbose >= 1
                        fprintf('Prevent clusters spanning sequencing gaps of length %sbp\n', vNumbers);
                    end
                end
                if vMaxInterGeneDist ~= -1
                    fprintf('Prevent clusters spanning intergenic distances longer or equal than %ibp\n', vMaxInterGeneDist);
                end
                if vPercentileForMaxInterGeneDist ~= -1
                    fprintf('Prevent clusters spanning the top %i%% intergenic distances\n', vPercentileForMaxInterGeneDist);
                end
                if vMaxInterGeneDistByMedian ~= -1
                    fprintf('Prevent clusters spanning intergenic distances longer or equal than %i * the median of gene sizes\n', vMaxInterGeneDistByMedian);
                end
            end
        end
        [vGeneLocation2, vGeneLocation_Att, vSeqGaps, vInterChrom, vBiomG_G] = ...
        f_introduce_HardB(vGeneLocation2, vGeneLocation_Att, ...
        vUniqueChrom, vInterChrom, vSeqGaps, vSeqGapSizesChromBreak, ...
        vMaxSeqGapSize, vMaxInterGeneDist, vPercentileForMaxInterGeneDist, vMaxInterGeneDistByMedian, vBiomG_G);
        
        vPartFinished = 12;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Populate Hypothetical genes
    if vPartFinished < 13
        if vVerbose >= 1
            fprintf(['Introduce Hypothetical genes for gap' ...
                's longer than %i bp\n'], vGeneLengthHypoMin);
        end
        if vInsertHypos > 0
            if vInsertHypos == 1
                if ~exist('vHypoAmount','var')
                    error('vHypoAmount does not exist. Needs to be set as an Argument.');
                end
                [vGeneLocation2, vGeneLocation_Att, ~, vBiomG_G] = ...
                f_introduceHypoGenes(vGeneLocation2, vGeneLocation_Att, ...
                vInterChrom, vSeqGaps, vGeneLengthHypoMin, vBiomG_G, vInsertHypos, vVerbose, 'HypoAmount', vHypoAmount);
            end
        end
        
        vPartFinished = 13;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Annotate Biomartfile with Enzymes (use Gene-protein mapping)
    if vPartFinished < 14
        if vVerbose >= 1
            fprintf('Annotate gene location file with enzyme information\n');
        end
        [vBiomG_E] = f_annotate_enyme_info_to_gene_location(vGeneLocation2, ...
            vGeneLocation_Att, vBiomG_G, vMapG_Gpgdb, vG2R, vVerbose);
        if vVerbose >= 2
            fprintf('Number of genes in location file identified as enzymes: %i\n', sum(vBiomG_E));
        end
        vPartFinished = 14;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Compute LD information
    if vPartFinished < 15
        if vVerbose >= 1
            fprintf('Compute local gene duplication per gene location\n');
        end
        [vG_LD, vG_LD_ClustIDs] = f_compute_LocalDuplication(vGeneLocation2,...
            vGeneLocation_Att, vG_MCL, vBiomG_G, vKbRangeLD, vGeneRangeLD, ...
            vVerbose);
        vPartFinished = 15;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end

    %% Perform clustering
    if vPartFinished < 16
        if vVerbose >= 1
            fprintf('Perform clustering\n');
        end
        [vClusters, vClusterBackground, vFinalStepsize, vStepSize_Max] = ...
            f_perform_clustering(vMinStepSize, vMaxStepSize, vGeneLocation2, ...
            vGeneLocation_Att, vBiomG_E, vEnzMinForClusters, vKbRangeLD, ...
            vGeneRangeLD, vG_MCL, vBiomG_G, vG2R, vMapG_Gpgdb, vCriterion, vVerbose);
        vPartFinished = 16;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Perform cutoff
    if vPartFinished < 17
        if vVerbose >= 1
            fprintf('Perform cuttof at %i%%\n', vTopPercentClusters);
        end
        [vClusters_Top] = f_perform_cutoff(vClusters, vClusterBackground, vTopPercentClusters);
        vPartFinished = 17;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Produce Gene Outputfile
    if vPartFinished < 18
        if vVerbose >= 1
            fprintf('Produce GeneOutputFile\n');
        end
        f_write_Gene_Output(vGeneOutputFile, vGeneLocation2, vGeneLocation_Att, vFinalStepsize, vStepSize_Max, vTopPercentClusters, vMD_to_annotate, ...
            vG2R, vMapG_Gpgdb, vBiomG_G, vReactions_all, vReactions_all_Att, vBiomG_E, ...
            vG_LD_ClustIDs, vG_LD, vClusters, vClusters_Top, vR_MD, vRpgdb_Rmd, ...
            vSTClasses, vR_STClasses, vSigOrTail, vR_SigOrTail, vVerbose, vOutputType);
        vPartFinished = 18;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1;
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
    
    %% Print cluster File
    if vPartFinished < 19
        if vVerbose >= 1
            fprintf('Produce ClusterOutputFile\n');
        end
        f_write_Cluster_Output(vClusterOutputFile, vClusters, vClusters_Top, vGeneLocation2, vGeneLocation_Att, vVerbose, vFinalStepsize, vOutputType);
        vPartFinished = 19;
        if vTempsaves == 1
            if vTempsavesOverwrite == 1
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif ~exist([vGeneLocation_File '_tempsave.mat'],'file')
                save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
            elseif exist([vGeneLocation_File '_tempsave.mat'],'file')
                vPartFinished_old_1 = vPartFinished;
                load([vGeneLocation_File '_tempsave.mat'],'vPartFinished')
                if vPartFinished_old_1 > vPartFinished
                    vPartFinished = vPartFinished_old_1; %#ok<NASGU>
                    clear vPartFinished_old
                    save([vGeneLocation_File '_tempsave.mat'],'-regexp','^(?!vRunAfterPart$).','-v7.3');
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function identifying wether reactions found in the pgdb are of a certain
% type of enyzme and if they are either signature or tailoring enzymes.
% Note, there have to be signature enzymes identified as 'signature' in the
% File vSignatureTailorFile.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vReactions_all: Matlab structure containing all the information
%                   (Attributes)that was searched for (described in
%                   vReactionsAttributes, default is {'UNIQUE-ID', ...
%                   'ENZYMATIC-REACTION', 'EC-NUMBER'};) in reactions.dat
%                   of the pgdb.
%   vReactions_all_Att: Attributes of vReactions_all, can be displayed by
%                       vReactions_all.(char(vReactions_all_Att(1,n)))
%   vSignatureTailorFile: TSV file describing the classes of enzymes that
%                         should be classified as signature or tailoring
%                         enzymes to identify clusters containing such.
% 
% Outputs: 
%   vSTClasses: Is a list of signature and tailoring enzyme classes.
%   vSigOrTail: Is a list of types of enzymes, e.g. signature, tailoring
%               and so on.
%   vR_STClasses: Mapping of Reactions in vReactions_all to the enzyme
%                 classes of vSTClasses.
%   vR_SigOrTail: Mapping of Reactions in vReactions_all to the enzyme
%                 types of vSigOrTail.
function [vSTClasses, vSigOrTail, vR_STClasses, vR_SigOrTail] = f_lable_Reactions_withSigandTail(vReactions_all, vReactions_all_Att, vSignatureTailorFile)
    vFIO = fopen(vSignatureTailorFile);
    vLine = fgetl(vFIO);
    vHeader = regexp(vLine,'\t','split');
    vInit = 20000;
    vSTdata = cell(vInit,size(vHeader,2));
    vLine = fgetl(vFIO);
    vL = 1;
    while ischar(vLine)
        vSplit = regexp(vLine,'\t','split');
        vSTdata(vL,:) = vSplit;
        vLine = fgetl(vFIO);
        vL = vL +1;
    end
    fclose(vFIO);
    vSTdata = vSTdata(1:(vL-1),:);
    
    vSTClasses = unique(vSTdata(:,2));
    vSigOrTail = unique(vSTdata(:,3));
    vMatchSTClasses = false(size(vSTdata,1),size(vSTClasses,1));
    vMatchSigOrTail = false(size(vSTdata,1),size(vSigOrTail,1));
    for vi = 1:size(vSTdata,1)
        vMatchSTClasses(vi,:) = strcmp(vSTClasses,vSTdata(vi,2));
        vMatchSigOrTail(vi,:) = strcmp(vSigOrTail,vSTdata(vi,3));
    end
    
    vSTdata_star_IDs = ~cellfun('isempty',regexp(vSTdata(:,1),'\.\*','once'));
    vSTdata_star = vSTdata(vSTdata_star_IDs,:);
    vSTdata = vSTdata(~vSTdata_star_IDs,:);
    vMatchSTClasses_star = vMatchSTClasses(vSTdata_star_IDs,:);
    vMatchSTClasses = vMatchSTClasses(~vSTdata_star_IDs,:);
    vMatchSigOrTail_star = vMatchSigOrTail(vSTdata_star_IDs,:);
    vMatchSigOrTail = vMatchSigOrTail(~vSTdata_star_IDs,:);
    
    vR_STClasses = false(size(vReactions_all.(char(vReactions_all_Att(1,1))),1),size(vSTClasses,1));
    vR_SigOrTail = false(size(vReactions_all.(char(vReactions_all_Att(1,1))),1),size(vSigOrTail,1));
    for vi = 1:size(vReactions_all.(char(vReactions_all_Att(1,1))),1)
        vIDs = strcmp(vSTdata(:,1),vReactions_all.(char(vReactions_all_Att(1,1)))(vi,1));
        for vj = 1:size(vReactions_all.(char(vReactions_all_Att(1,3))),2)
            vIDs = vIDs | strcmp(vSTdata(:,1),regexprep(vReactions_all.(char(vReactions_all_Att(1,3)))(vi,vj),'^EC-',''));
        end
        vIDs_find = find(vIDs,1,'first');
        if size(vIDs_find,1)>0
            vR_STClasses(vi,:) = vMatchSTClasses(vIDs_find,:);
            vR_SigOrTail(vi,:) = vMatchSigOrTail(vIDs_find,:);
        end
    end
    vColSig = strcmpi(vSigOrTail,'signature');
    vTemp = regexp(vSTdata_star(:,1),'\.\*','split');
    for vi = 1:size(vTemp)
        vTemp2 = vTemp{vi,1};
        vSTdata_star(vi,1) = vTemp2(1,1);
    end
    for vi = 1:size(vSTdata_star,1)
        vIDs = strncmpi(regexprep(vReactions_all.(char(vReactions_all_Att(1,3)))(:,1),'^EC-',''),vSTdata_star(vi,1),size(char(vSTdata_star(vi,1)),2));
        for vj = 2:size(vReactions_all.(char(vReactions_all_Att(1,3))),2)
            vIDs = vIDs | strncmpi(regexprep(vReactions_all.(char(vReactions_all_Att(1,3)))(:,vj),'^EC-',''),vSTdata_star(vi,1),size(char(vSTdata_star(vi,1)),2));
        end
        vIDs = vIDs & ~vR_SigOrTail(vColSig,:);
        vR_STClasses(vIDs,vMatchSTClasses_star(vi,:)) = 1;
        vR_SigOrTail(vIDs,vMatchSigOrTail_star(vi,:)) = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that computes wether a genes is locally duplicated, and if so,
% which MCL cluster it is associated with. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vG_MCL: Maping between the genes in the conversion-ID file and the
%           MCL-clustering file.
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%   vKbRangeLD: Genes that are encoded more than this bp apart, are not
%               called local duplicates.
%   vGeneRangeLD: Number of genes that space two genes from each other.
%                 Genes that have more than this many intervening genes
%                 between themare not called local duplicates.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vG_LD: A vector of length of the number of Genes in vGeneLocation2,
%          being 1 if the gene was considered to be locally duplicated (a
%          gene was within vKbRangeLD bp distance and not spaced more than
%          with vGeneRangeLD genes from another gene with the same MCL
%          clustering given by vG_MCL.
%   vG_LD_ClustIDs: The MCL clusters of the gene (comma separated if more
%                   than 1 peptide is associated to a gene
function [vG_LD, vG_LD_ClustIDs] = f_compute_LocalDuplication(vGeneLocation2, vGeneLocation_Att, vG_MCL, vBiomG_G, vKbRangeLD, vGeneRangeLD, vVerbose)
    vNgenes = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
    vG_LD = false(vNgenes,1);
    vG_LD_ClustIDs = cell(vNgenes,1);
    vGStartbp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))));
    vGEndbp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))));
    for vi = 1:vNgenes
        if vVerbose >= 2
            fprintf('Checking Local duplication of Gene %i of %i\n', vi, vNgenes);
        end
        if ~strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1),'Hypothetical_Gene_',18)
            vStartG = vi - vGeneRangeLD - 1;
            vEndG = vi + vGeneRangeLD + 1;
            if vStartG < 1
                vStartG = 1;
            end
            if vEndG > vNgenes
                vEndG = vNgenes;
            end
            vIDs = false(vNgenes,1);
            vIDs(vStartG:vEndG,1) = 1;
            vIDs(vi,1) = 0;
            vIDs = vIDs & strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4))),vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1));

            vStartbp = vGEndbp(vi,1)-vKbRangeLD;
            vEndbp = vGStartbp(vi,1)+vKbRangeLD;
            vIDs = vIDs & vGStartbp >= vStartbp & vGEndbp <= vEndbp;
            vIDs_find = find(vIDs);
            
            vMCL = sum(vG_MCL(sum(vBiomG_G(vi,:),1)>0,:),1)>0;
            vMCL_find = find(vMCL);
            clear vTempMCLString_array
            if size(vMCL_find,2)~=0
                vTempMCLString_array = cellstr(['MCL_' num2str(vMCL_find(1,1))]);
            end
            for vj = 2:size(vMCL_find,2)
                vTempMCLString_array = [vTempMCLString_array; cellstr(['MCL_' num2str(vMCL_find(1,1))])];
            end
            if exist('vTempMCLString_array','var')
                vTempMCLString_array = unique(vTempMCLString_array);
            end
            vTempMCLString = '';
            if exist('vTempMCLString_array','var')
                if size(vTempMCLString_array,1)~=0
                    vTempMCLString = char(vTempMCLString_array);
                end
                for vj = 2:size(vTempMCLString_array,1)
                    vTempMCLString = [vTempMCLString ', ' char(vTempMCLString_array)];
                end
            end
            vG_LD_ClustIDs(vi,1) = cellstr(vTempMCLString);
            vN_MCL = size(vMCL,2);
            
            vStop = 0;
            vj = 1;
            while vj <= size(vIDs_find,1) && vStop == 0
                if sum(vMCL)==0
                    fprintf('Warning: %s has no entry in MCL clustering file.\n',char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1)))
                    vG_LD(vi,1) = -1;
                elseif sum(vMCL == (sum(vG_MCL(sum(vBiomG_G(vIDs_find(vj,1),:),1)>0,:),1)>0)) == vN_MCL
                    vG_LD(vi,1) = 1;
                    vStop = 1;
                end
                vj = vj + 1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that introduces Hypothetical genes if a sequencing gap reached a
% critical size vGeneLengthHypoMin.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vInterChrom: A matrix defining chromosome, start bp, end bp, gene
%                before and gene after a intergenic region.
%   vSeqGaps: A string array (cell) matrix containing the
%             chromosome/scaffold, the identification if it is coding (A)
%             or masked (N) DNA, start bp and length of a masked region.
%             In addition there is the information stored in with
%             intergenic region the gap is encoded and its start and end
%             bp, corrected for maximal size (= intergenic region).
%   vGeneLengthHypoMin: The minimal length of bp of masked DNA that should
%                       be populated by two hypothetical genes (one minus,
%                       and one plus strand)
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%   vInsertHypos: [0,1, can be extended with more code] If set to 1, then
%                 Hypothetical genes are introduced. Then vHypoAmount needs
%                 to be set.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)% 
% Outputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vInterChrom: A matrix defining chromosome, start bp, end bp, gene
%                before and gene after a intergenic region.
%   vHypoAmount: Number (0-n, default 2) of Hypothetical genes that should
%                be inserted. Only needed if vInsertHypos is set to 1.
function [vGeneLocation2, vGeneLocation_Att, vInterChrom, vBiomG_G] = ...
    f_introduceHypoGenes(vGeneLocation2, vGeneLocation_Att, vInterChrom, vSeqGaps, vGeneLengthHypoMin, ...
    vBiomG_G, vInsertHypos, vVerbose, varargin)
    vCurrArgin = 8;
    if nargin > vCurrArgin
        vi = 1;
        while vi <= nargin - vCurrArgin
            if strcmp(varargin{1,vi},'HypoAmount')
                if ischar(varargin{1,vi+1})
                    vHypoAmount = str2num(varargin{1,vi + 1});
                else
                    vHypoAmount = varargin{1,vi + 1};
                end
                vi = vi + 1;
            end
            vi = vi + 1;
        end
    end
    
    if vInsertHypos > 0
        %Check where and how many Hypothetical genes should be placed.
        if vInsertHypos == 1
            vSeqGaps_5 = str2num(char(vSeqGaps(:,5)));
            vSeqGaps_7 = str2num(char(vSeqGaps(:,7)));
            vHypoToPlace = zeros(size(vSeqGaps_5,1),1);
            vBreakAfter = zeros(size(vSeqGaps_5,1),1);
            if exist('vHypoAmount','var')
                vNi = size(vSeqGaps,1);
                for vi = vNi:-1:1
                    if vVerbose >= 2
                        fprintf('Estimate number of gaps: %i of %i\n', vNi - vi + 1, vNi);
                    end
                    vBreakAfter(vi,1) = vInterChrom(vSeqGaps_5(vi,1),4);

                    if vSeqGaps_7(vi,1) >= vGeneLengthHypoMin
                        vHypoToPlace(vi,1) = vHypoAmount;
                    end
                end
            else
                error('vHypoAmount not defined, needs to be defined for InsertHypos = 1.');
            end
        else
            error('Not other option (>1) programmed.');
        end
        
        %Sort the Breaks:
        [~, vBreakAfter_sort_ID] = sort(vBreakAfter);
    
        %Remove several sequencing gaps within one intergenic region.
        vBreakafterkill = false(size(vSeqGaps,1),1);
        vOld_Breakafter_sort_ID = vBreakAfter_sort_ID(1,1);
        for vi = 2:size(vSeqGaps,1)
            if vBreakAfter(vBreakAfter_sort_ID(vi,1),1) == vBreakAfter(vOld_Breakafter_sort_ID,1)
                vBreakafterkill(vBreakAfter_sort_ID(vi,1),1) = 1;
                %It is possible that a sequencing gap does not meet the
                %minimal Hypo length, but the next one does. We populate 2
                %genes in total once one of them meets the criterion.
                if vHypoToPlace(vOld_Breakafter_sort_ID,1) < vHypoToPlace(vBreakAfter_sort_ID(vi,1),1)
                    vHypoToPlace(vOld_Breakafter_sort_ID,1) = vHypoToPlace(vBreakAfter_sort_ID(vi,1),1);
                end
            else
                vOld_Breakafter_sort_ID = vBreakAfter_sort_ID(vi,1);
            end
        end
        vSeqGaps = vSeqGaps(~vBreakafterkill,:);
%         vSeqGaps_5 = vSeqGaps_5(~vBreakafterkill,:);
%         vSeqGaps_7 = vSeqGaps_7(~vBreakafterkill,:);
        vBreakAfter = vBreakAfter(~vBreakafterkill,:);
        vHypoToPlace = vHypoToPlace(~vBreakafterkill,:);
        
        [~, vBreakAfter_sort_ID] = sort(vBreakAfter);
        
        %Increase the size of the gene locations.
        for vk = 1:size(vGeneLocation_Att,2)
            vGeneLocation3.(char(vGeneLocation_Att(1,vk))) = cell(size(vGeneLocation2.(char(vGeneLocation_Att(1,vk))),1) + sum(vHypoToPlace), ...
                size(vGeneLocation2.(char(vGeneLocation_Att(1,vk))),2));
        end
    
        %Increase the size of all variables that are connected to the gene
        %locations.
        vBiomG_G_3 = logical(spalloc(size(vBiomG_G,1) + sum(vHypoToPlace),size(vBiomG_G,2),size(vBiomG_G,2)*2));
        
        %Populate the new gene location variable with the hypothetical genes.
        vEnd_GeneLocation2 = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
        vEnd_GeneLocation3 = size(vGeneLocation3.(char(vGeneLocation_Att(1,1))),1);
        for vi = size(vSeqGaps,1):-1:1
            if vVerbose >= 2
                fprintf('Populate hypothetical genes: %i of %i\n', vNi - vi + 1, vNi);
            end
            if vHypoToPlace(vBreakAfter_sort_ID(vi,1),1) ~= 0
                vStart_GeneLocation2 = vBreakAfter(vBreakAfter_sort_ID(vi,1),1) + 1;
                vStart_GeneLocation3 = vEnd_GeneLocation3 - (vEnd_GeneLocation2 - vStart_GeneLocation2);
                for vl = 1:vHypoToPlace(vBreakAfter_sort_ID(vi,1),1)
                    vPlaceLastGeneEnd = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vStart_GeneLocation2-1,1)));
                    vPlaceBetweenGenes = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vStart_GeneLocation2,1))) - vPlaceLastGeneEnd - 1;
                    vPlaceToInsert = vPlaceLastGeneEnd + floor(vPlaceBetweenGenes/(vHypoToPlace(vBreakAfter_sort_ID(vi,1),1)+1)*vl);
                    for vk = 1:size(vGeneLocation_Att,2)
                        vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3:vEnd_GeneLocation3,:) = ...
                        vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation2:vEnd_GeneLocation2,:);
                        if vk == 1
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = cellstr(['Hypothetical_Gene_' num2str(sum(vHypoToPlace(1:vBreakAfter_sort_ID(vi,1)-1,1))+vl)]);
                        elseif vk == 2
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = cellstr(num2str(vPlaceToInsert));
                        elseif vk == 3
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = cellstr(num2str(vPlaceToInsert));
                        elseif vk == 4
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vBreakAfter(vBreakAfter_sort_ID(vi,1),1),1);
                        elseif vk == 5
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vBreakAfter(vBreakAfter_sort_ID(vi,1),1),1);
                        elseif vk == 6
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vBreakAfter(vBreakAfter_sort_ID(vi,1),1),1);
                        end
                    end
                end
                if vVerbose >= 2
                    fprintf('Inserting Hypothetical genes: Copy old %i to %i to new %i to %i\n', vStart_GeneLocation2, vEnd_GeneLocation2, vStart_GeneLocation3, vEnd_GeneLocation3);
                end
                vBiomG_G_3(vStart_GeneLocation3:vEnd_GeneLocation3,:) = vBiomG_G(vStart_GeneLocation2:vEnd_GeneLocation2,:);
                vIDs_to_change = vInterChrom(:,5)>vBreakAfter(vBreakAfter_sort_ID(vi,1),1);
                vInterChrom(vIDs_to_change,5) = vInterChrom(vIDs_to_change,5) + vHypoToPlace(vBreakAfter_sort_ID(vi,1),1);
                vIDs_to_change = vInterChrom(:,4)>vBreakAfter(vBreakAfter_sort_ID(vi,1),1);
                vInterChrom(vIDs_to_change,4) = vInterChrom(vIDs_to_change,4) + vHypoToPlace(vBreakAfter_sort_ID(vi,1),1);
                vEnd_GeneLocation2 = vStart_GeneLocation2-1;
                vEnd_GeneLocation3 = vStart_GeneLocation3-(1+vHypoToPlace(vBreakAfter_sort_ID(vi,1),1));
            end
        end
        vStart_GeneLocation3 = 1;
        vStart_GeneLocation2 = 1;
        vBiomG_G_3(vStart_GeneLocation3:vEnd_GeneLocation3,:) = vBiomG_G(vStart_GeneLocation2:vEnd_GeneLocation2,:);
        for vk = 1:size(vGeneLocation_Att,2)
            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(1:vEnd_GeneLocation3,:) = ...
            vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(1:vEnd_GeneLocation2,:);
        end
        vGeneLocation2 = vGeneLocation3;
        vBiomG_G = vBiomG_G_3;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that prevents formation of clusters Sequencing gaps that
% exceed vMaxSeqGapSize or that are defined not to be crossed by
% vHardBreakGaps. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vUniqueChrom: A list with unique entries of all chromosome and scaffold
%                 names.
%   vInterChrom: A matrix defining chromosome, start bp, end bp, gene
%                before and gene after a intergenic region.
%   vSeqGaps: A string array (cell) matrix containing the
%             chromosome/scaffold, the identification if it is coding (A)
%             or masked (N) DNA, start bp and length of a masked region.
%             In addition there is the information stored in with
%             intergenic region the gap is encoded and its start and end
%             bp, corrected for maximal size (= intergenic region).
%   vHardBreakGaps: A list of numbers. The masked nucleotides sequences of
%                   these specific lengths are preventing cluster formation
%                   (Default is an empty list []).
%   vMaxSeqGapSize: Maximal masked nucleotide region (in bp) that a cluster
%                   is allowed to bridge (default is -1, thus inactive).
%   vMaxInterGeneDist: Maximal intergenic distance that still can be 
%                      gene cluster.
%   vPercentileForMaxInterGeneDist: Prevents clusters to be spanned over 
%                                   the top n% of intergenic distances.
%   vMaxInterGeneDistByMedian: Maximal intergenic distance that still can 
%                              be spanned by a gene cluster: this number
%                              times the median of the gene sizes.
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%
% Outputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vSeqGaps: A string array (cell) matrix containing the
%             chromosome/scaffold, the identification if it is coding (A)
%             or masked (N) DNA, start bp and length of a masked region.
%             In addition there is the information stored in with
%             intergenic region the gap is encoded and its start and end
%             bp, corrected for maximal size (= intergenic region).
%   vInterChrom: A matrix defining chromosome, start bp, end bp, gene
%                before and gene after a intergenic region.
function [vGeneLocation2, vGeneLocation_Att, vSeqGaps, vInterChrom, vBiomG_G] = f_introduce_HardB(vGeneLocation2, vGeneLocation_Att, vUniqueChrom, vInterChrom, vSeqGaps, vHardBreakGaps, vMaxSeqGapSize, vMaxInterGeneDist, vPercentileForMaxInterGeneDist, vMaxInterGeneDistByMedian, vBiomG_G, varargin)
    %for every occurence of bad gap or gap tat exeeds the max size, reform
    %the chromosome strand such that it is nomore the same as before (split
    %chromosome
    vCurrArgin = 11;
    if nargin > vCurrArgin
        vi = 1;
        while vi <= nargin - vCurrArgin
            
            vi = vi +1;
        end
    end
    
    vSeqSize = str2num(char(vSeqGaps(:,end)));
    vSeqHard = false(size(vSeqSize,1),1);
    for vi = 1:size(vHardBreakGaps,2)
        vSeqHard = vSeqHard | vSeqSize == vHardBreakGaps(1,vi);
    end
    if vMaxSeqGapSize ~= -1
        vSeqHard = vSeqHard | vSeqSize >= vMaxSeqGapSize;
    end
    vSeqGaps_Hard = vSeqGaps(vSeqHard,:);
    vSeqGaps = vSeqGaps(~vSeqHard,:);
    
    %Check maximal intergenic distance via percent
    if vPercentileForMaxInterGeneDist ~= -1
        %Change From Median to Percentile
        %vMultiplierInterGeneDist = median(vInterChrom(:,3)-vInterChrom(:,2)+1);
        %vMaxInterGeneDist_temp = vMultiplierInterGeneDist*vPercentileForMaxInterGeneDist;
        %Use 95% Percentile Cutoff Instead
        vMaxInterGeneDist_temp = prctile(vInterChrom(:,3)-vInterChrom(:,2)+1, vPercentileForMaxInterGeneDist);
        fprintf('Percentile used for intergenic distance: %ibp\nClusters broken when spanning distance larger: %ibp.\n', vPercentileForMaxInterGeneDist, vMaxInterGeneDist_temp);
        if vMaxInterGeneDist == -1
            vMaxInterGeneDist = vMaxInterGeneDist_temp;
        else
            vMaxInterGeneDist = min([vMaxInterGeneDist,vMaxInterGeneDist_temp]);
        end
    end
    
    %Check maximal intergenic distance via median gene size
    if vMaxInterGeneDistByMedian ~= -1
        vMaxInterGeneDist_temp = prctile(str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(:,1))) - ...
                                         str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(:,1))), 50);
        if vMaxInterGeneDist == -1
            vMaxInterGeneDist = vMaxInterGeneDist_temp*vMaxInterGeneDistByMedian;
        else
            vMaxInterGeneDist = min([vMaxInterGeneDist,vMaxInterGeneDist_temp*vMaxInterGeneDistByMedian]);
        end
    end
    
    %Check intergenic distances
    if vMaxInterGeneDist ~= -1
        vG_start = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(:,1)));
        vG_end = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(:,1)));
        vAct_Crom = '';
        vNG = size(vGeneLocation2.(char(vGeneLocation_Att(1,4))),1);
        for vi = vNG:-1:2
            if strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi-1,1),vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1))
                if vG_end(vi,1) - vG_start(vi-1,1) >= vMaxInterGeneDist
                    if ~strcmp(vAct_Crom,vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1))
                        vChrom_counter = 1;
                        vAct_Crom = vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1);
                    else
                        vChrom_counter = vChrom_counter + 1;
                    end
                    vIDs = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1));
                    vIDs(1:vi-1,1) = 0;
                    vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vIDs,1) = cellstr([char(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1)) '_' num2str(vChrom_counter)]);
                end
            end
        end
    end
    
    %Check physical hard breaks
    vGeneLocation2.([char(vGeneLocation_Att(1,4)) '_old']) = vGeneLocation2.(char(vGeneLocation_Att(1,4)));
    vGeneLocation_Att = [vGeneLocation_Att, [char(vGeneLocation_Att(1,4)) '_old']];
    for vi = size(vUniqueChrom,1):-1:1
        vTempHardGaps = vSeqGaps_Hard(strcmp(vSeqGaps_Hard(:,1),vUniqueChrom(vi,1)),:);
        if size(vTempHardGaps,1) > 0
            vTempChrom = cell(size(vTempHardGaps,1)+1,1);
            for vj = 1:size(vTempChrom,1)
                vTempChrom(vj,1) = cellstr([char(vUniqueChrom(vi,1)) '_' num2str(vj)]);
            end
            for vj = size(vTempHardGaps,1):-1:1
                vBreakAfter = vInterChrom(str2num(char(vTempHardGaps(vj,5))),4);
                vTemp_IDs = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vUniqueChrom(vi,1));
                vTemp_IDs(1:vBreakAfter) = 0;
                vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vTemp_IDs,1) = vTempChrom(vj + 1,1);
            end
            vTemp_IDs = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vUniqueChrom(vi,1));
            vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vTemp_IDs,1) = vTempChrom(1,1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Functions that uses the output of f_get_Sequencing_Gaps and reformulates 
% it such that it is clear in which intergenic region the sequening gap
% lies.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vFile: Fasta file of the masked genome nucleotide sequences.
%   vInterGeneLocation: A matrix defining chromosome, start bp, end bp,
%                       gene before and gene after a intergenic region.
%   vUniqueChrom: A list with unique entries of all chromosome and scaffold
%                 names.
%   vOverwriteSeqGapInfo: If set to 1, the code enforces masked nucleotide
%                         analysis for finding sequencing gaps (default is
%                         0, no enforcement).
%   vUnmaskedDNA: Option to ignore the error if you use unmasked DNA files.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% Optional Inputs:
%   vSeqGapSizesChromBreak: A list of numbers. The masked nucleotides 
%                           sequences of these specific lengths are 
%                           preventing cluster formation (Default is an 
%                           empty list []). If this is given, then also
%                           vGeneLengthHypoMin has to be set. Will lead to
%                           prescreening of the sequencing gaps.
%   vGeneLengthHypoMin: A size that defines the minimal length of a
%                       sequencing gap that should still be processed. If
%                       this is bigger than 0, then all sequencing gaps
%                       smaller than this are deleted, unless they are
%                       among the list of vSeqGapSizesChromBreak, which
%                       also needs to be set. Default is -1, not active.
% 
% Outputs: 
%   vSeqGaps: A string array (cell) matrix containing the
%             chromosome/scaffold, the identification if it is coding (A)
%             or masked (N) DNA, start bp and length of a masked region.
%             In addition there is the information stored in with
%             intergenic region the gap is encoded and its start and end
%             bp, corrected for maximal size (= intergenic region).
function vSeqGaps = f_annotate_Sequencing_Gaps(vFile, vInterGeneLocation, vUniqueChrom, vOverwriteSeqGapInfo, vPara, vUnmaskedDNA, vVerbose, varargin)
    vCurrNargin = 7;
    vGeneLengthHypoMin = -1;
    vSeqGapSizesChromBreak = [];
    vSet = 0;
    if nargin > vCurrNargin
        vi = 1;
        while vi <= nargin - vCurrNargin
            if ischar(varargin{1,vi})
                if strcmp(varargin{1,vi},'SeqGapSizesChromBreak')
                    vSeqGapSizesChromBreak = varargin{1,vi+1};
                    vSet = vSet + 1;
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'GeneLengthHypoMin')
                    vGeneLengthHypoMin = varargin{1,vi+1};
                    vSet = vSet + 1;
                    vi = vi + 1;
                end
            end
            vi = vi + 1;
        end
    end
    if vGeneLengthHypoMin ~= -1 && vSet ~= 2
        error('f_annotate_Sequencing_Gaps needs both SeqGapSizesChromBreak and GeneLengthHypoMin to be set or to be unset.');
    end

    vMaxInit = 7000000;
    if ispc()
        if ~exist([vFile '_GAPOutput'], 'file') || vOverwriteSeqGapInfo
            f_get_Sequencing_Gaps(vFile, vPara, vVerbose, vOverwriteSeqGapInfo, vUnmaskedDNA)
        end
    else
        [~,vOutput] = system(['ls' ' ' '''' vFile '_GAPOutput''']);
        if strncmp(vOutput,'ls: cannot access', 17) || vOverwriteSeqGapInfo
            f_get_Sequencing_Gaps(vFile, vPara, vVerbose, vOverwriteSeqGapInfo, vUnmaskedDNA)
        end
    end
    
    vFIO = fopen([vFile '_GAPOutput']);
    vSeqGaps = cell(vMaxInit,4);
    vLine = fgetl(vFIO);
    vL = 1;
    while ischar(vLine)
        %Split the line into four columns using the space-character as
        %delimiter.
        vSplit = regexp(vLine,' ','split');
        if size(vSplit,2) ~= 4
            error('%s is corrupt at line %i. It should contain 4 columns split by spaces, it contains %i.\n\n', [vFile '_GAPOutput'],vL, size(vSplit,2));
        end
        
        if vGeneLengthHypoMin == -1 || sum(vSeqGapSizesChromBreak == str2num(char(vSplit(1,4)))) > 0 || (vGeneLengthHypoMin ~= -1 && vGeneLengthHypoMin <= str2num(char(vSplit(1,4))))
            if vL > vMaxInit
                vSeqGaps = [vSeqGaps; cell(1,4)];
            end
            vSeqGaps(vL,:) = vSplit;
            vL = vL + 1;
        end
        vLine = fgetl(vFIO);
    end
    fclose(vFIO);
    vSeqGaps = vSeqGaps(1:(vL-1),:);
    
    %find the intergenic region of the gap, if the gap is present within a
    %gene model, get rid of the gap.
    vSeqGaps = [vSeqGaps, cell(size(vSeqGaps,1),3)];
    vSeqGaps_3 = str2num(char(vSeqGaps(:,3)));
    vSeqGaps_4 = str2num(char(vSeqGaps(:,4)));
    vSeqGaps_4 = vSeqGaps_3 + vSeqGaps_4 - 1;
    
    vnvi = size(vSeqGaps,1);
    vKillSeqGaps = false(vnvi,1); % Check which sequences should not be processed
    vKillSeqGaps = vKillSeqGaps | ~strcmp(vSeqGaps(:,2),'N'); % Do not process all non N sequences.
    % Do not process sequences that do not lie on a chromosome we want.
    vHasGeneLocChrom = false(vnvi,1); 
    for vi = 1:size(vUniqueChrom,1)
        vHasGeneLocChrom = vHasGeneLocChrom | strcmp(vSeqGaps(:,1),vUniqueChrom(vi,1));
    end
    vBadUniqueChrom = unique(vSeqGaps(~vHasGeneLocChrom,1));
    for vi = 1:size(vBadUniqueChrom,1)
        warning('DNA Fasta file ''%s'' contains chromosomes or scaffolds not covered with the gene position file: %s\n', vFile, char(vBadUniqueChrom(vi,1)));
    end
    vKillSeqGaps = vKillSeqGaps | ~vHasGeneLocChrom;
    vSeqGaps = vSeqGaps(~vKillSeqGaps,:);
    vSeqGaps_3 = vSeqGaps_3(~vKillSeqGaps,:);
    vSeqGaps_4 = vSeqGaps_4(~vKillSeqGaps,:);
    vnvi = size(vSeqGaps,1);
    vKillSeqGaps = false(vnvi,1); % Check which sequences should not be processed
    for vi =1:vnvi
        vChromID = find(strcmp(vUniqueChrom,vSeqGaps(vi,1)));
        vInterGeneLocationIDs = ...
            vInterGeneLocation(:,1) == vChromID & ...
            vInterGeneLocation(:,2) <= vSeqGaps_3(vi,1) & ...
            vInterGeneLocation(:,3) >= vSeqGaps_3(vi,1) & ...
            vInterGeneLocation(:,2) <= vSeqGaps_4(vi,1) & ...
            vInterGeneLocation(:,3) >= vSeqGaps_4(vi,1);
        if sum(vInterGeneLocationIDs) == 0
            vKillSeqGaps(vi,1) = 1;
        else
            vInterGeneLocationIDs_find = find(vInterGeneLocationIDs);
            vSeqGaps(vi,5) = cellstr(num2str(vInterGeneLocationIDs_find));
            vSeqGapsStart = max([vInterGeneLocation(vInterGeneLocationIDs_find,2),str2num(char(vSeqGaps(vi,3)))]);
            vSeqGapsEnd = min([vInterGeneLocation(vInterGeneLocationIDs_find,3),str2num(char(vSeqGaps(vi,3))) + str2num(char(vSeqGaps(vi,4))) - 1]);
            vSeqGaps(vi,6) = cellstr(num2str(vSeqGapsStart));
            vSeqGaps(vi,7) = cellstr(num2str(vSeqGapsEnd-vSeqGapsStart+1));
        end
    end
    vSeqGaps = vSeqGaps(~vKillSeqGaps,:);
%     vSeqGaps_3 = vSeqGaps_3(~vKillSeqGaps,:);
%     vSeqGaps_4 = vSeqGaps_4(~vKillSeqGaps,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that screens a DNA fasta file for masked streches of DNA (masked
% by N).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vFile: Fasta file of the masked genome nucleotide sequences.
%   vPara: Defines if multiple cpus should be used for computation of parts
%          of algorithm. Active if larger than 1.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
%   vOverwriteSeqGapInfo: If set to 1, the code enforces masked nucleotide
%                         analysis for finding sequencing gaps (default is
%                         0, no enforcement).
%   vUnmaskedDNA: Option to ignore the error if you use unmasked DNA files.
% 
% Outputs: 
%   None (File is generated with the path of the first input argument, 
%        modifying the gene name)
function f_get_Sequencing_Gaps(vFile, vPara, vVerbose, vOverwriteSeqGapInfo, vUnmaskedDNA)
    if ispc()
        if vOverwriteSeqGapInfo == 1 || ~exist([vFile '_GAPOutput'], 'file')
            fprintf('Finding sequencing gaps is only set up for linux. Please precompute the the sequencing gaps (masked stretches of DNA) (e.g. on linux server).\n');
            fprintf('Please use the follwing commands on Linux (use which and programname to find your paths):\n');
            vAWK_install_path = 'awk_install_path';
            fprintf('%s -f enter_new_line_characters_in_fasta_file.awk ''%s''>''%s_temp1''\n',vAWK_install_path, vFile, vFile);
            fprintf('tr -d ''\r'' < ''%s_temp1''>''%s_temp2''\n', vFile, vFile);
            fprintf('tr -d ''\n'' < ''%s_temp2''>''%s_temp3''\n', vFile, vFile);
            fprintf('%s -f get_new_line_in.awk ''%s_temp3''>''%s_temp4''\n', vAWK_install_path, vFile, vFile);
            fprintf('%s -f get_positions_of_gap.awk ''%s_temp4''>''%s_GAPOutput''\n', vAWK_install_path, vFile, vFile);
            fprintf('Once you have done this, please store the resulting %s_GAPOutput file in the accoring folder of the DNA fasta file and rerun Plantcluster finder.\n', vFile); 
            error('Run aborted because PC is used in combination with sequencing gap finding\n');
        end
    else
        if vVerbose >= 1
            fprintf('Check awk installation.\n');
        end
        [~, vAWK_install_path] = system('which awk');
        if strncmp(vAWK_install_path,'/usr/bin/which: no awk',22)
            error('No awk installed. Abort.\n');
        else
            vAWK_install_path = regexprep(vAWK_install_path,'\n','');
        end

        % Linearize the fasta files of the genomes
        if vVerbose >= 1
            fprintf('Get information gaps in genomes.\n');
        end
        if vPara <= 1
            f_run_awk_skripts(vFile, vAWK_install_path, vVerbose, vOverwriteSeqGapInfo);
        else
            %Split Files into scaffolds
            vFile_ID = 1;
            vFIO = fopen(vFile);
            vFIO_temp = fopen([vFile '_' num2str(vFile_ID)],'w');
            vLine = fgetl(vFIO);
            if strncmp(vLine,'>',1)
                fprintf(vFIO_temp,[vLine '\n']);
            else
                error('First line in DNA Fasta file is not a header\n');
            end
            vLine = fgetl(vFIO);
            while ischar(vLine)
                if strncmp(vLine,'>',1)
                    fclose(vFIO_temp);
                    vFile_ID = vFile_ID + 1;
                    vFIO_temp = fopen([vFile '_' num2str(vFile_ID)],'w');
                    fprintf(vFIO_temp,[vLine '\n']);
                else
                    fprintf(vFIO_temp,[vLine '\n']);
                end
                vLine = fgetl(vFIO);
            end
            fclose(vFIO_temp);
            fclose(vFIO);
            
            %open matlabpool
            vPool = parpool('local',vPara-1);
            
            %Run the commands individually with parfor
            parfor vi = 1:size(vFile_ID)
                f_run_awk_skripts([vFile '_' num2str(vi)], vAWK_install_path, vVerbose, vOverwriteSeqGapInfo);
            end
            
            %close matlabpool
            delete(vPool);
            
            %Build result files together
            vFIO = fopen([vFile '_GAPOutput_count.txt'],'w');
            for vi = 1:size(vFile_ID)
                vFIO_temp = fopen([vFile '_' num2str(vi) '_GAPOutput_count.txt']);
                vLine = fgetl(vFIO_temp);
                while ischar(vLine)
                    fprintf(vFIO, [vLine '\n']);
                    vLine = fgetl(vFIO_temp);
                end
                fclose(vFIO_temp);
            end
            fclose(vFIO);
        end
    end

    %Produce the Gap analysis file
    [~,vOutput] = system(['ls' ' ' '''' vFile '_GAPOutput_count.txt''']);
    if vOverwriteSeqGapInfo == 1 || strncmp(vOutput,'ls: cannot access', 17)
        f_analyze_PlantClusterGapFile([vFile '_GAPOutput'], vUnmaskedDNA, vVerbose);
    else
        if vVerbose >= 1
            fprintf('File %s exists, not overwritten.\n', ['''' vFile '_GAPOutput_count.txt''']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that uses awk to process the masked DNA fasta file to search for
% sequencing gaps.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vFile: Full path to a dna Fasta file. Note, the folder the DNA fasta
%          file needs to be writable for the temporary and result files
%          generated.
%   vAWK_install_path: Path to the installation of awk.
%   vVerbose: Tells the code on how much of output should be printed out 
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
%   vOverwriteSeqGapInfo: If set to 1, the code enforces masked nucleotide
%                         analysis for finding sequencing gaps (default is
%                         0, no enforcement).
% 
% Outputs: 
%   None (Temporary (_temp1, _temp2, _temp3, _temp4) and Gapouptu
%        (_GAPOutput) files generated.
function f_run_awk_skripts(vFile, vAWK_install_path, vVerbose, vOverwriteSeqGapInfo)
    [~,vOutput] = system(['ls' ' ' '''' vFile '_temp1''']);
    if vOverwriteSeqGapInfo == 1 || strncmp(vOutput,'ls: cannot access', 17)
        if vVerbose >= 2
            fprintf('Next command running:\n');
            fprintf([vAWK_install_path ' ' '-f enter_new_line_characters_in_fasta_file.awk' ' ' '''' vFile '''>''' vFile '_temp1''']);
            fprintf('\n');
        end
        [vStatus, vOutput] = system([vAWK_install_path ' ' '-f enter_new_line_characters_in_fasta_file.awk' ' ' '''' vFile '''>''' vFile '_temp1''']);
        if vVerbose >= 2
            fprintf('Status was %i\n', vStatus);
            fprintf('Output was:\n');
            fprintf('%s\n', vOutput);
        end
    else
        if vVerbose >= 1
            fprintf('File %s exists, not overwritten.\n', ['''' vFile '_temp1''']);
        end
    end
    [~,vOutput] = system(['ls' ' ' '''' vFile '_temp2''']);
    if vOverwriteSeqGapInfo == 1 || strncmp(vOutput,'ls: cannot access', 17)
        if vVerbose >= 2
            fprintf('Next command running:\n');
            fprintf(['tr -d ''\r'' <' ' ' '''' vFile '_temp1''>''' vFile '_temp2''']);
            fprintf('\n');
        end
        [vStatus, vOutput] = system(['tr -d ''\r'' <' ' ' '''' vFile '_temp1''>''' vFile '_temp2''']);
        if vVerbose >= 2
            fprintf('Status was %i\n', vStatus);
            fprintf('Output was:\n');
            fprintf('%s\n', vOutput);
        end
    else
        if vVerbose >= 1
            fprintf('File %s exists, not overwritten.\n', ['''' vFile '_temp2''']);
        end
    end
    [~,vOutput] = system(['ls' ' ' '''' vFile '_temp3''']);
    if vOverwriteSeqGapInfo == 1 || strncmp(vOutput,'ls: cannot access', 17)
        if vVerbose >= 2
            fprintf('Next command running:\n');
            fprintf(['tr -d ''\n'' <' ' ' '''' vFile '_temp2''>''' vFile '_temp3''']);
            fprintf('\n');
        end
        [vStatus, vOutput] = system(['tr -d ''\n'' <' ' ' '''' vFile '_temp2''>''' vFile '_temp3''']);
        if vVerbose >= 2
            fprintf('Status was %i\n', vStatus);
            fprintf('Output was:\n');
            fprintf('%s\n', vOutput);
        end
    else
        if vVerbose >= 1
            fprintf('File %s exists, not overwritten.\n', ['''' vFile '_temp3''']);
        end
    end

    [~,vOutput] = system(['ls' ' ' '''' vFile '_temp4''']);
    if vOverwriteSeqGapInfo == 1 || strncmp(vOutput,'ls: cannot access', 17)
        if vVerbose >= 2
            fprintf('Next command running:\n');
            fprintf([vAWK_install_path ' ' '-f get_new_line_in.awk' ' ' '''' vFile '_temp3''>''' vFile '_temp4''']);
            fprintf('\n');
        end
        [vStatus, vOutput] = system([vAWK_install_path ' ' '-f get_new_line_in.awk' ' ' '''' vFile '_temp3''>''' vFile '_temp4''']);
        if vVerbose >= 2
            fprintf('Status was %i\n', vStatus);
            fprintf('Output was:\n');
            fprintf('%s\n', vOutput);
        end
    else
        if vVerbose >= 1
            fprintf('File %s exists, not overwritten.\n', ['''' vFile '_temp4''']);
        end
    end

    % Search for sequencing gaps (searches for continuous stretches of
    % ATCGatcg* and Nn s, and lists them in a new file
    [~,vOutput] = system(['ls' ' ' '''' vFile '_GAPOutput''']);
    if vOverwriteSeqGapInfo == 1 || strncmp(vOutput,'ls: cannot access', 17)
        if vVerbose >= 2
            fprintf('Next command running:\n');
            fprintf([vAWK_install_path ' ' '-f get_positions_of_gap.awk' ' ' '''' vFile '_temp4''>''' vFile '_GAPOutput''']);
            fprintf('\n');
        end
        [vStatus, vOutput] = system([vAWK_install_path ' ' '-f get_positions_of_gap.awk' ' ' '''' vFile '_temp4''>''' vFile '_GAPOutput''']);
        if vVerbose >= 2
            fprintf('Status was %i\n', vStatus);
            fprintf('Output was:\n');
            fprintf('%s\n', vOutput);
        end
    else
        if vVerbose >= 1
            fprintf('File %s exists, not overwritten.\n', ['''' vFile '_GAPOutput''']);
        end
    end
    [status,msg,~] = fileattrib([vFile '_GAPOutput'],'+w +x','a');
    if status~= 1
        fprintf('Could not change file attribute of %s\n', ['''' vFile '_GAPOutput''']);
        error('Message was %s', msg);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that reads in the metabolic domains of reactions
% and gives a mapping between reactions in these files and the metabolic
% domains that were asked to be annotated (usually a subset of all) by
% vMD_to_annotate.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vGenesMetabolicDomains_File: TSV of the metabolic domains of reactions.
%                                Check if file is outdated (needs to
%                                contain all reactions that are present in
%                                the PGDB that have genes annotated to).
%   vMD_to_annotate: A list of Metabolic domains that should be analyzed.
%                    For example {'Amines and Polyamines Metabolism'; ...
%                    'Amino Acids Metabolism'
% 
% Outputs: 
%   vRxnMetaDom: Matlab structure of all the items stored in the Metabolic
%                domain file
%   vRxnMetaDom_Att: Attributes of vRxnMetaDom (lets you access vRxnMetaDom
%                    with vRxnMetaDom.(char(vRxnMetaDom_Att(1,n))))
%   vR_MD: Maping between all Reactions stored in the Metabolic domain file
%          and the metabolic domains that were identified to be of interest
%          in vMD_to_annotate.
function [vRxnMetaDom, vRxnMetaDom_Att, vR_MD] = f_get_metabolic_domains(vGenesMetabolicDomains_File, vMD_to_annotate)
    [vRxnMetaDom, vRxnMetaDom_Att] = f_extract_results_with_header(vGenesMetabolicDomains_File); 
    for vi = 1:size(vRxnMetaDom_Att,2)
        for vj = 1:size(vRxnMetaDom.(char(vRxnMetaDom_Att(1,vi))),2)
            vRxnMetaDom.(char(vRxnMetaDom_Att(1,vi)))(:,vj) = strtrim(vRxnMetaDom.(char(vRxnMetaDom_Att(1,vi)))(:,vj));
        end
    end
    
    vUnique_MD = unique(vRxnMetaDom.(char(vRxnMetaDom_Att(1,3))));
    for vi = 4:size(vRxnMetaDom_Att,2)
        vUnique_MD = unique([vUnique_MD;vRxnMetaDom.(char(vRxnMetaDom_Att(1,vi)))]);
    end
    
    for vi = 1:size(vUnique_MD,1)
        if sum(strcmp(vMD_to_annotate, vUnique_MD(vi,1))) == 0
            warning('Discarding metabolic domain information of %s\n', char(vUnique_MD(vi,1)));
        end
    end
    
    vR_MD = false(size(vRxnMetaDom.(char(vRxnMetaDom_Att(1,1))),1),size(vMD_to_annotate,1));
    for vi = 1:size(vMD_to_annotate,1)
        for vj = 3:size(vRxnMetaDom_Att,2)
            vR_MD(strcmp(vRxnMetaDom.(char(vRxnMetaDom_Att(1,vj)))(:,1),vMD_to_annotate(vi,1)),vi) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that reads out the MCL clustering information and gives every
% gene one or several MCL clusters. This is then used for local
% duplication.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vFile: MCL-clustering output file (I20).
%   vMapG: Unique list of gene IDs stored in the conversion file.
%   vMapP: Unique list of protein IDs stored in the conversion file.
%   vMapG_MapP: Maping between vMapG and vMapP
%   vVerbose: Tells the code on how much of output should be printed out 
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vG_MCL: Maping between the genes in the conversion-ID file and the
%           MCL-clustering file.
%   vProtein_IDs: A cell array of all protein-IDs that are stored in the
%                 protein fasta file. only active if
%                 vRemoveNonProtLocations is 1. Else the array has no rows.
function [vG_MCL,vProtein_IDs] = f_read_MCL_clustering_File(vFile, vMapG, vMapP, vMapG_MapP, vRemoveNonProtLocations, vVerbose)
    if vRemoveNonProtLocations == 1
        vProtinit = 300000;
        vProtein_IDs = cell(vProtinit,1);
        vFIO = fopen(regexprep(vProtein_sequence_FastaFile, regexptranslate('escape', '.MCLTAB.dump.I20'),''));
        vLine = fgetl(vFIO);
        vPID = 1;
        while ischar(vLine)
            if strncmp(vLine,'>',1)
                vTempLine = regexp(vLine, ' ', 'split');
                if vPID <= vProtinit
                    vProtein_IDs(vPID,1) = cellstr(regexprep(vTempLine, '^>', ''));
                else
                    vProtein_IDs = [vProtein_IDs; cellstr(regexprep(vTempLine, '^>', ''))];
                end
                vPID = vPID + 1;
            end
            vLine = fgetl(vFIO);
        end
        fclose(vFIO);
        vProtein_IDs = vProtein_IDs(1:(vPID-1),:);
    else
        vProtein_IDs = cell(0,1);
    end
    
    vInitMax = 100000;
    vFile_Data = cell(vInitMax,1);
    vFIO = fopen(vFile);
    vLine = fgetl(vFIO);
    vL = 1;
    if ischar(vLine)
        vFile_Data(vL,1) = cellstr(vLine);
    end
    vLine = fgetl(vFIO);
    vL = 2;
    while ischar(vLine)
        if vVerbose >= 2
            fprintf('MCL_Line %i\n', vL);
        end
        if vL <= vInitMax
            vFile_Data(vL,1) = cellstr(vLine);
        else
            vFile_Data = [vFile_Data; cellstr(vLine)];
        end
        vLine = fgetl(vFIO);
        vL = vL + 1;
    end
    fclose(vFIO);
    vFile_Data = vFile_Data(1:vL-1,1);
    
    vG_MCL = false(size(vMapG_MapP,1),vL-1);
    nvMapP = size(vMapP,1);
    for vi = 1:nvMapP
        if vVerbose >= 2
            fprintf('MCL: process protein %i of %i\n', vi, nvMapP);
        end
        vIDs = ~cellfun('isempty',regexp(vFile_Data, regexptranslate('escape', vMapP(vi,1)), 'once'));
        vG_MCL(vMapG_MapP(:,vi)==1,vIDs) = 1;
    end
    
    for vi = 1:size(vMapG,1)
        if sum(vG_MCL(vi,:)) == 0
            fprintf('Warning: %s seems not to be part of the MCL clustering.\n',char(vMapG(vi,1)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that uses a protein fasta file of a genome annotation to create
% a blastdb, and then does an all against all blast and uses the result to
% do a MCL clustering. This is used for local duplication information.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vProtein_sequence_FastaFile: Protein sequence file of the genome.
%                                Should only contain Protein IDs as
%                                header e.g. >protein1-ID
%   vMapP: Unique list of protein IDs stored in the conversion file.
%   vMapG_MapP: Maping between vMapG and vMapP
%   vPara: Defines if multiple cpus should be used for computation of parts
%          of algorithm. Active if larger than 1.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vG_MCL: Maping between the genes in the conversion-ID file and the
%           MCL-clustering file.
%   vProtein_IDs: A cell array of all protein-IDs that are stored in the
%                 protein fasta file. only active if
%                 vRemoveNonProtLocations is 1. Else the array has no rows.
function [vG_MCL, vProtein_IDs] = f_get_MCL_clustering(vProtein_sequence_FastaFile, vMapG, vMapP, vMapG_MapP, vPara, vOverwriteMCLClustering, vRemoveNonProtLocations, vVerbose)
    if ~exist([vProtein_sequence_FastaFile '.MCLTAB.dump.I20'],'file') || vOverwriteMCLClustering == 1
        if ispc()
            fprintf('MCL clustering is only set up for linux. Please precompute the MCL clustering file (e.g. on linux server).\n');
            fprintf('Please use the follwing commands on Linux (use which and programname to find your paths):\n');
            vMakeblastdb_install_path = 'makeblastdb_install_path';
            vBlastp_install_path = 'blastp_install_path';
            vMcxload_install_path = 'mcxload_install_path';
            vMcl_install_path = 'mcl_install_path';
            vMcxdump_install_path = 'mcxdump_install_path';
            fprintf('%s -in %s -dbtype prot -out %s.blastDB -hash_index\n',vMakeblastdb_install_path,vProtein_sequence_FastaFile,vProtein_sequence_FastaFile);
            fprintf('%s -query %s -db %s.blastDB -evalue 1 -out %s.blastOUT -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send"\n', vBlastp_install_path, vProtein_sequence_FastaFile, vProtein_sequence_FastaFile, vProtein_sequence_FastaFile);
            fprintf('cut -f 1,2,3 %s.blastOUT > %s.MCLIN\n', vProtein_sequence_FastaFile, vProtein_sequence_FastaFile);
            fprintf('%s -abc %s.MCLIN --stream-mirror --stream-neg-log10 -stream-tf ''ceil(200)'' -o %s.MCLMATRIX -write-tab %s.MCLTAB\n',vMcxload_install_path, vProtein_sequence_FastaFile, vProtein_sequence_FastaFile, vProtein_sequence_FastaFile);
            fprintf('%s %s.MCLMATRIX -I 2\n',vMcl_install_path, vProtein_sequence_FastaFile);
            fprintf('%s -icl %s.MCLMATRIX.I20 -tabr %s.MCLTAB -o %s.MCLTAB.dump.I20\n', vMcxdump_install_path, vProtein_sequence_FastaFile, vProtein_sequence_FastaFile, vProtein_sequence_FastaFile);
            fprintf('Please use a Protein fasta file that contains only protein Identifiers as headers.\n')
            fprintf('Once you have done this, please give the %s.MCLTAB.dump.I20 as Filename for the proteinfasta file and use option ''MCLClusterFile'' 1 to start PlantClusterFinder.\n',vProtein_sequence_FastaFile); 
            error('Run aborted because PC is used in combination with MCL clustering');
        end

        [~, vMakeblastdb_install_path] = system('which makeblastdb');
        if strncmp(vMakeblastdb_install_path,'/usr/bin/which: no makeblastdb',size('/usr/bin/which: no makeblastdb',2))
            error('No makeblastdb installed. Abort.\n');
        else
            vMakeblastdb_install_path = regexprep(vMakeblastdb_install_path,'\n','');
        end

        [~, vBlastp_install_path] = system('which blastp');
        if strncmp(vBlastp_install_path,'/usr/bin/which: no blastp',size('/usr/bin/which: no blastp',2))
            error('No blastp installed. Abort.\n');
        else
            vBlastp_install_path = regexprep(vBlastp_install_path,'\n','');
        end

        [~, vMcxload_install_path] = system('which mcxload');
        if strncmp(vMcxload_install_path,'/usr/bin/which: no mcxload',size('/usr/bin/which: no mcxload',2))
            error('No mcxload installed. Abort.\n');
        else
            vMcxload_install_path = regexprep(vMcxload_install_path,'\n','');
        end

        [~, vMcl_install_path] = system('which mcl');
        if strncmp(vMcl_install_path,'/usr/bin/which: no mcl',size('/usr/bin/which: no mcl',2))
            error('No mcl installed. Abort.\n');
        else
            vMcl_install_path = regexprep(vMcl_install_path,'\n','');
        end

        [~, vMcxdump_install_path] = system('which mcxdump');
        if strncmp(vMcxdump_install_path,'/usr/bin/which: no mcxdump',size('/usr/bin/which: no mcxdump',2))
            error('No mcxdump installed. Abort.\n');
        else
            vMcxdump_install_path = regexprep(vMcxdump_install_path,'\n','');
        end

        [~, ~] = system([vMakeblastdb_install_path ' ' '-in' ' ' vProtein_sequence_FastaFile ' ' '-dbtype prot -out' ' ' vProtein_sequence_FastaFile '.blastDB -hash_index']);
        if vPara > 1
            [~, ~] = system([vBlastp_install_path ' ' '-query' ' ' vProtein_sequence_FastaFile ' ' '-db' ' ' vProtein_sequence_FastaFile '.blastDB -evalue 1 -out' ' ' vProtein_sequence_FastaFile '.blastOUT -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send" -num_threads ' num2str(vPara - 1)]);
        else
            [~, ~] = system([vBlastp_install_path ' ' '-query' ' ' vProtein_sequence_FastaFile ' ' '-db' ' ' vProtein_sequence_FastaFile '.blastDB -evalue 1 -out' ' ' vProtein_sequence_FastaFile '.blastOUT -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send"']);
        end
        [~, ~] = system(['cut -f 1,2,3' ' ' vProtein_sequence_FastaFile '.blastOUT >' ' ' vProtein_sequence_FastaFile '.MCLIN']);
        [~, ~] = system([vMcxload_install_path ' ' '-abc' ' ' vProtein_sequence_FastaFile '.MCLIN --stream-mirror --stream-neg-log10 -stream-tf ''ceil(200)'' -o' ' ' vProtein_sequence_FastaFile '.MCLMATRIX -write-tab' ' ' vProtein_sequence_FastaFile '.MCLTAB']);
        [~, ~] = system([vMcl_install_path ' ' vProtein_sequence_FastaFile '.MCLMATRIX -I 2 -o' ' ' vProtein_sequence_FastaFile '.MCLMATRIX.I20']);
        [~, ~] = system([vMcxdump_install_path ' ' '-icl' ' ' vProtein_sequence_FastaFile '.MCLMATRIX.I20 -tabr' ' ' vProtein_sequence_FastaFile '.MCLTAB -o' ' ' vProtein_sequence_FastaFile '.MCLTAB.dump.I20']);
    end
    
    [vG_MCL, vProtein_IDs] = f_read_MCL_clustering_File([vProtein_sequence_FastaFile '.MCLTAB.dump.I20'], vMapG, vMapP, vMapG_MapP, vRemoveNonProtLocations, vVerbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function reading out information of pgdb flatfiles genes.dat,
% proteins.dat, enzrxns.dat and reactions.dat, maping genes
% to reactions directly.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vPGDB_FlatFileFolder: Full path to the folder where the flatfiles of
%                         the species pgdb is stored.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% Optional Inputs:
%   'GenesAttributes': Attributes of genes.dat that should be read out.
%                      Default is {'UNIQUE-ID','PRODUCT','ACCESSION-1', ...
%                      'ACCESSION-2'}
%   'vProteinsAttributes': Attributes of proteins.dat that should be read
%                          out. Default is {'UNIQUE-ID','CATALYZES', ...
%                          'COMPONENTS'}
%   'EnzrxnsAttributes': Attributes of enzrxns.dat that should be read
%                        out. Default is {'UNIQUE-ID','REACTION','ENZYME'}
%   'ReactionsAttributes': Attributes of reactions.dat that should be read
%                          out. Default is {'UNIQUE-ID',
%                          'ENZYMATIC-REACTION', 'EC-NUMBER'}
%   'G2P': Defines which attributes in genes.dat should be matched to which
%          attributes in proteins.dat. Default is: [2,1];
%   'P2P': Defines which attributes in proteins.dat should be matched to
%          which attributes in proteins.dat (Protein complexes). Default
%          is: [1,3];
%   'P2E': Defines which attributes in proteins.dat should be matched to
%          which attributes in enzrxns.dat. Default is: [2,1];
%   'E2R': Defines which attributes in enzrxns.dat should be matched to
%          which attributes in reactions.dat. Default is: [2,1];
% 
% Outputs: 
%   vGenes_all: Matlab structure containing all the information
%               (Attributes)that was searched for (described in
%               vGenesAttributes, default is {'UNIQUE-ID','PRODUCT', ...
%               'ACCESSION-1','ACCESSION-2'}) in genes.dat of the pgdb.
%   vGenes_all_Att: Attributes of the Matlab structure vGenes_all, can be
%                   displayed by vGenes_all.(char(vGenes_all_Att(1,n)))
%   vReactions_all: Matlab structure containing all the information
%                   (Attributes) that was searched for (described in
%                   vReactionsAttributes, default is {'UNIQUE-ID', ...
%                   'ENZYMATIC-REACTION', 'EC-NUMBER'}) in reactions.dat
%                   of the pgdb.
%   vReactions_all_Att: Attributes of vReactions_all, can be displayed by
%                       vReactions_all.(char(vReactions_all_Att(1,n)))
%   vG2R: Maping of all genes in the Genes_all to the Reactions in
%         Reactions_all
function [vGenes_all, vGenes_all_Att, vReactions_all, vReactions_all_Att, vG2R] = f_map_pgdb_reactions_to_genes(vPGDB_FlatFileFolder, vVerbose, varargin)
    vCurrNargin = 2;
    
    % Attributes of Flatfiles that should be read out
    vGenesAttributes = {'UNIQUE-ID','PRODUCT','ACCESSION-1','ACCESSION-2'};
    vProteinsAttributes = {'UNIQUE-ID','CATALYZES','COMPONENTS'};
    vEnzrxnsAttributes = {'UNIQUE-ID','REACTION','ENZYME'};
    vReactionsAttributes = {'UNIQUE-ID', 'ENZYMATIC-REACTION', 'EC-NUMBER'};
    
    % Which Attributes should be used to match different Flatfiles to each
    % other.
    vG2P = [2,1];
    vP2P = [1,3];
    vP2E = [2,1];
    vE2R = [2,1];
    vR2Pa = [1,2]; %#ok<NASGU>
    
    %Use inputs to write over defaults
    if nargin > vCurrNargin
        vi = 1;
        while vi <= nargin - vCurrNargin
            if ischar(varargin{1,vi})
                if strcmp(varargin{1,vi},'GenesAttributes')
                    vGenesAttributes = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'vProteinsAttributes')
                    vProteinsAttributes = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'EnzrxnsAttributes')
                    vEnzrxnsAttributes = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'ReactionsAttributes')
                    vReactionsAttributes = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'G2P')
                    vG2P = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'P2P')
                    vP2P = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'P2E')
                    vP2E = varargin{1,vi+1};
                    vi = vi + 1;
                elseif strcmp(varargin{1,vi},'E2R')
                    vE2R = varargin{1,vi+1};
                    vi = vi + 1;
                end
            end
            vi = vi + 1;
        end
    end
    
    %make sure that the unique IDs were not forgot.
    if ~strcmp(vGenesAttributes{1,1},'UNIQUE-ID')
        vGenesAttributes = {'UNIQUE-ID', vGenesAttributes};
    end
    if ~strcmp(vProteinsAttributes{1,1},'UNIQUE-ID')
        vProteinsAttributes = {'UNIQUE-ID', vProteinsAttributes};
    end
    if ~strcmp(vEnzrxnsAttributes{1,1},'UNIQUE-ID')
        vEnzrxnsAttributes = {'UNIQUE-ID', vEnzrxnsAttributes};
    end
    if ~strcmp(vReactionsAttributes{1,1},'UNIQUE-ID')
        vReactionsAttributes = {'UNIQUE-ID', vReactionsAttributes};
    end
    
    %Read all attributes and store them
    if vVerbose >= 1
        fprintf('Read in genes.dat\n');
    end
    [vGenes_all, vGenes_all_Att] = f_read_in_PGDB_Flat([vPGDB_FlatFileFolder 'genes.dat'],vGenesAttributes);
    if vVerbose >= 1
        fprintf('Read in proteins.dat\n');
    end
    [vProteins_all, vProteins_all_Att] = f_read_in_PGDB_Flat([vPGDB_FlatFileFolder 'proteins.dat'],vProteinsAttributes);
    if vVerbose >= 1
        fprintf('Read in enzrxns.dat\n');
    end
    [vEnzrxns_all, vEnzrxns_all_Att] = f_read_in_PGDB_Flat([vPGDB_FlatFileFolder 'enzrxns.dat'],vEnzrxnsAttributes);
    if vVerbose >= 1
        fprintf('Read in reactions.dat\n');
    end
    [vReactions_all, vReactions_all_Att] = f_read_in_PGDB_Flat([vPGDB_FlatFileFolder 'reactions.dat'], vReactionsAttributes);
    
    %Match flat files to each other
    if vVerbose >= 1
        fprintf('Match genes to proteins\n');
    end
    [vG2_P] = f_match_pgdb_flat(vGenes_all,vGenes_all_Att,vProteins_all,vProteins_all_Att,vG2P,vVerbose);
    if vVerbose >= 1
        fprintf('Match proteins to protein-complexes\n');
    end
    [vP_P] = f_match_pgdb_flat(vProteins_all,vProteins_all_Att,vProteins_all,vProteins_all_Att,vP2P,vVerbose);
    for vi = 1:size(vP_P,1)
        vP_P(vi,vi) = 1;
    end
    vSUM_old = -1;
    vSUM = sum(sum(vP_P));
    while vSUM ~= vSUM_old
        vSUM_old = vSUM;
        for vi = 1:size(vP_P,1)
            vIDs1 = (vP_P(vi,:) == 1 | vP_P(:,vi)' == 1);
            vP_P(vi,vIDs1) = 1;
            vP_P(vIDs1,vi) = 1;
        end
        vSUM = sum(sum(vP_P));
    end
    if vVerbose >= 1
        fprintf('Match proteins and protein complexes to enzymereactions\n');
    end
    [vP_E] = f_match_pgdb_flat(vProteins_all,vProteins_all_Att,vEnzrxns_all,vEnzrxns_all_Att,vP2E,vVerbose);
    if vVerbose >= 1
        fprintf('Match enzymereactions to reactions\n');
    end
    [vE_R2] = f_match_pgdb_flat(vEnzrxns_all,vEnzrxns_all_Att,vReactions_all,vReactions_all_Att,vE2R,vVerbose);
    
    %Match Genes to reactions directly.
    if vVerbose >= 1
        fprintf('Match genes to to reactions\n');
    end
    vG2R = false(size(vG2_P,1),size(vE_R2,2));
    for Gene_all_ID = 1:size(vG2R,1)
        vP_ID = vG2_P(Gene_all_ID,:)==1;
        vP2_ID = [sum(vP_P(vP_ID,:),1)>0]'; %#ok<NBRAK>
        vE_ID = [sum(vP_E(vP2_ID,:),1)>0]'; %#ok<NBRAK>
        vR_ID = [sum(vE_R2(vE_ID,:),1)>0]'; %#ok<NBRAK>
        vG2R(Gene_all_ID,vR_ID) = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that maps two flat files of a pgdb to each other.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs: 
%   vFlat1: Matlab structure of a read out pgdb Flat file.
%   vFlatAttri1: Attributes of the first structure
%   vFlat2: Matlab structure of a read out pgdb Flat file.
%   vFlatAttri2: Attributes of the second structure
%   vMatchings: Definition of which attributes in structure 1 should be
%               matched to which attributes of structure 2
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vF1_F2: A mapping between items of the first pgdb file to the second
%           one.
function [vF1_F2] = f_match_pgdb_flat(vFlat1,vFlatAttri1,vFlat2,vFlatAttri2,vMatchings,vVerbose)
    vF1_F2 = false(size(vFlat1.(char(vFlatAttri1(1,1))),1),...
                  size(vFlat2.(char(vFlatAttri2(1,1))),1));
    for vMatchings_ID = 1:size(vMatchings,1)
        for vvFlat1_ID = 1:size(vFlat1.(char(vFlatAttri1(1,vMatchings(vMatchings_ID,1)))),1)
            if vVerbose >= 2
                fprintf('Match PGDB FlatFiles ID %i of %i\n',vvFlat1_ID,size(vFlat1.(char(vFlatAttri1(1,1))),1));
            end
            
            vFlat2_temp = false(size(vFlat2.(char(vFlatAttri2(1,vMatchings(vMatchings_ID,2)))),1),1);
            for vFlat1_Dim = 1:size(vFlat1.(char(vFlatAttri1(1,vMatchings(vMatchings_ID,1)))),2)
                for vFlat2_Dim = 1:size(vFlat2.(char(vFlatAttri2(1,vMatchings(vMatchings_ID,2)))),2)
                    vFlat2_temp = vFlat2_temp | strcmp(vFlat2.(char(vFlatAttri2(1,vMatchings(vMatchings_ID,2))))(:,vFlat2_Dim),vFlat1.(char(vFlatAttri1(1,vMatchings(vMatchings_ID,1))))(vvFlat1_ID,vFlat1_Dim));
                end
            end
            vF1_F2(vvFlat1_ID,vFlat2_temp) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that reads in the gene, transcript, and protein ID conversion
% file and stores the IDs as unique lists and mapings between genes and
% transcripts aswell as genes and proteins.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneTranscriptProtein_mapping_File: TSV file with a header describing
%                                        gene, transcript, and protein
%                                        identifier and then for all the
%                                        identfiers the listings (gene-ID,
%                                        transcript-ID, protein-ID).
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vMapG: Unique list of gene IDs stored in the conversion file.
%   vMapT: Unique list of transcript IDs stored in the conversion file.
%   vMapP: Unique list of protein IDs stored in the conversion file.
%   vMapG_MapT: Maping between vMapG and vMapT
%   vMapG_MapP: Maping between vMapG and vMapP
function [vMapG, vMapT, vMapP, vMapG_MapT, vMapG_MapP] = f_read_in_ConversionFile(vGeneTranscriptProtein_mapping_File, vVerbose)%, vPara, varargin)
%     if vPara >= 1
%         vPool = varargin{1,1};
%     end

    vFIO = fopen(vGeneTranscriptProtein_mapping_File);
    vInitialize = 250000;
    vLine = fgetl(vFIO);
    vGTP_Header = regexp(vLine,'\t','split'); %#ok<NASGU>
    vLine = fgetl(vFIO);
    vSplit = regexp(vLine,'\t','split');
    if vVerbose >= 2
        fprintf('Read conversion ids 1\n');
    end
    for vi = 1:size(vSplit,2)
        vGTP.(['Col' num2str(vi)]) = cell(vInitialize,1);
        vGTP.(['Col' num2str(vi)])(1,1) = vSplit(1,vi);
    end
    vLine = fgetl(vFIO);
    vL = 2;
    while ischar(vLine)
        if vVerbose >= 2
            fprintf('Read conversion ids %i\n', vL);
        end
        vSplit = regexp(vLine,'\t','split');
        for vi = 1:size(vSplit,2)
            if vL <= vInitialize
                vGTP.(['Col' num2str(vi)])(vL,1) = vSplit(1,vi);
            else
                vGTP.(['Col' num2str(vi)]) = [vGTP.(['Col' num2str(vi)]);vSplit(1,vi)];
            end
        end
        vLine = fgetl(vFIO);
        vL = vL + 1;
    end
    fclose(vFIO);
    for vi = 1:size(vSplit,2)
        vGTP.(['Col' num2str(vi)]) = vGTP.(['Col' num2str(vi)])(1:vL-1,:);
    end
    
    vMapG = unique(vGTP.Col1(:,1));
    vMapT = unique(vGTP.Col2(:,1));
    vMapP = unique(vGTP.Col3(:,1));
    vMapG_MapT = false(size(vMapG,1),size(vMapT,1));
    vMapG_MapP = false(size(vMapG,1),size(vMapP,1));
    vnvi = size(vGTP.Col1,1);
    for vi = 1:vnvi
        if vVerbose >= 2
            fprintf('Map conversion ids %i of %i to each other\n', vi, vnvi);
        end
%         if vPara == 1
            vMapG_ID = strcmp(vMapG,vGTP.Col1(vi,1));
            vMapT_ID = strcmp(vMapT,vGTP.Col2(vi,1));
            vMapP_ID = strcmp(vMapP,vGTP.Col3(vi,1));
            vMapG_MapT(vMapG_ID,vMapT_ID) = 1;
            vMapG_MapP(vMapG_ID,vMapP_ID) = 1;
%         else
%             vMapG_ID_find = find(strcmp(vMapG,vGTP.Col1(vi,1)));
%             vMapT_ID_find = find(strcmp(vMapT,vGTP.Col2(vi,1)));
%             vMapP_ID_find = find(strcmp(vMapP,vGTP.Col3(vi,1)));
%             parfor vj = vMapG_ID_find
%                 for vk = vMapT_ID_find
%                     vMapG_MapT(vj,vk) = 1;
%                 end
%                 for vk = vMapP_ID_find
%                     vMapG_MapP(vj,vk) = 1;
%                 end
%             end
%         end
    end
    vSum_T = sum(vMapG_MapT,1);
    vSum_P = sum(vMapG_MapT,1);

    if sum(vSum_T>1) > 0 || sum(vSum_P>1) > 0
        if sum(vSum_T>1) > 0
            warning('Some transcripts match to multiple genes in the ID conversion file:\n');
            vSum_T_ID = find(vSum_T>1);
            for vi = 1:size(vSum_T_ID,1)
                fprintf('Please check transcript %s\n', char(vMapT(vSum_T_ID(vi,1))));
            end
        end
        if sum(vSum_P>1) > 0
            warning('Some proteins match to multiple genes in the ID conversion file:\n');
            vSum_P_ID = find(vSum_P>1);
            for vi = 1:size(vSum_P_ID,1)
                fprintf('Please check protein %s\n', char(vMapT(vSum_P_ID(vi,1))));
            end
        end
        error('Please clean Protein ID file.\n');
    end
    if vVerbose >= 1
        fprintf('Number of gene-ID to transcript-ID mappings: %i\n', sum(vSum_T));
        fprintf('Number of Gene to Protein mappings: %i\n', sum(vSum_P));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that maps the reactions stored in the pgdbs to the reactions
% stored in the metabolic domain file.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vG2R: Maping of all genes in the Genes_all to the Reactions in
%         Reactions_all
%   vR_MD: Maping between all Reactions stored in the Metabolic domain file
%          and the metabolic domains that were identified to be of interest
%          in vMD_to_annotate.
%   vReactions_all: Matlab structure containing all the information
%                   (Attributes)that was searched for (described in
%                   vReactionsAttributes, default is {'UNIQUE-ID', ...
%                   'ENZYMATIC-REACTION', 'EC-NUMBER'};) in reactions.dat
%                   of the pgdb.
%   vReactions_all_Att: Attributes of vReactions_all, can be displayed by
%                       vReactions_all.(char(vReactions_all_Att(1,n)))
%   vRxnMetaDom: Matlab structure of all the items stored in the Metabolic
%                domain file
%   vRxnMetaDom_Att: Attributes of vRxnMetaDom (lets you access vRxnMetaDom
%                    with vRxnMetaDom.(char(vRxnMetaDom_Att(1,n))))
% 
% Outputs: 
%   vRpgdb_Rmd: A maping between reactions stored in the pgdb and the
%               reactions given in the metabolic domain file

function [vRpgdb_Rmd] = f_matchPGDBs_to_MDs(vG2R, vR_MD, vReactions_all, vReactions_all_Att, vRxnMetaDom, vRxnMetaDom_Att)
    vRpgdb_Rmd = false(size(vG2R,2),size(vR_MD,1));
    for vi = 1:size(vReactions_all.(char(vReactions_all_Att(1,1))),1)
        vRpgdb_Rmd(vi,strcmp(vRxnMetaDom.(char(vRxnMetaDom_Att(1,1))),vReactions_all.(char(vReactions_all_Att(1,1)))(vi,1))) = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that maps the conversion file gene IDs to the IDs in the
% genes.dat of the pgdb
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vMapG: Unique list of gene IDs stored in the conversion file.
%   vMapT: Unique list of transcript IDs stored in the conversion file.
%   vMapP: Unique list of protein IDs stored in the conversion file.
%   vMapG_MapT: Maping between vMapG and vMapT
%   vMapG_MapP: Maping between vMapG and vMapP
%   vGenes_all: Matlab structure containing all the information
%               (Attributes)that was searched for (described in
%               vGenesAttributes, default is {'UNIQUE-ID','PRODUCT', ...
%               'ACCESSION-1','ACCESSION-2'};) in genes.dat of the pgdb.
%   vGenesDat_GeneInfo: List of Attributes of genes.dat. Defines where gene
%                       information should be searched for. Default is
%                       {'Unique-ID','Accession-1', 'Accession-2'}.
%   vPGDBIdsToMap: A string that contains G and or T and or P and defines
%                  if the entries represent Genes and or transcript ids and
%                  or proteins ids of the conversion file.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vMapG_Gpgdb: A maping between Genes in the conversion file and genes
%                stored in the pgdb genes.dat file.
function vMapG_Gpgdb = f_map_PGDBs_to_conversion_file(vMapG, vMapT, vMapP, vMapG_MapT, vMapG_MapP, vGenes_all, vGenesDat_GeneInfo, vPGDBIdsToMap, vVerbose)
    vPGDBIdsToMap_bin = false(3,1);
    if ~isempty(regexp(vPGDBIdsToMap,'G','once'))
        vPGDBIdsToMap_bin(1,1) = 1;
    end
    if ~isempty(regexp(vPGDBIdsToMap,'T','once'))
        vPGDBIdsToMap_bin(2,1) = 1;
    end
    if ~isempty(regexp(vPGDBIdsToMap,'P','once'))
        vPGDBIdsToMap_bin(3,1) = 1;
    end
    
    vMapG_Gpgdb = false(size(vMapG,1),size(vGenes_all.(char(vGenesDat_GeneInfo(1,1))),1));
    if vPGDBIdsToMap_bin(1,1) == 1
        for vi = 1:size(vMapG,1)
            for vj = 1:size(vGenesDat_GeneInfo,2)
                if vVerbose >= 2
                    fprintf('Conversion gene-ID %i, pgdb-ID %i\n', vi, vj);
                end
                for vk = 1:size(vGenes_all.(char(vGenesDat_GeneInfo(1,vj))),2)
                    vMapG_Gpgdb(vi,strcmpi(vGenes_all.(char(vGenesDat_GeneInfo(1,vj)))(:,vk),vMapG(vi,1))) = 1;
                end
            end
        end
    end
    
    if vPGDBIdsToMap_bin(2,1) == 1
        for vi = 1:size(vMapT,1)
            for vj = 1:size(vGenesDat_GeneInfo,2)
                if vVerbose >= 2
                    fprintf('Conversion transcript-ID %i, pgdb-ID %i\n', vi, vj);
                end
                for vk = 1:size(vGenes_all.(char(vGenesDat_GeneInfo(1,vj))),2)
                    vMapG_Gpgdb(vMapG_MapT(:,vi),strcmpi(vGenes_all.(char(vGenesDat_GeneInfo(1,vj)))(:,vk),vMapT(vi,1))) = 1;
                end
            end
        end
    end
    
    if vPGDBIdsToMap_bin(3,1) == 1
        for vi = 1:size(vMapP,1)
            for vj = 1:size(vGenesDat_GeneInfo,2)
                if vVerbose >= 2
                    fprintf('Conversion protein-ID %i, pgdb-ID %i\n', vi, vj);
                end
                for vk = 1:size(vGenes_all.(char(vGenesDat_GeneInfo(1,vj))),2)
                    vMapG_Gpgdb(vMapG_MapP(:,vi),strcmpi(vGenes_all.(char(vGenesDat_GeneInfo(1,vj)))(:,vk),vMapP(vi,1))) = 1;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that removes protein and transcript information of the gene
% location file.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneLocation: Matlab structure about the gene location information,
%                  generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vMapG: Unique list of gene IDs stored in the conversion file.
%   vMapT: Unique list of transcript IDs stored in the conversion file.
%   vMapP: Unique list of protein IDs stored in the conversion file.
%   vMapG_MapT: Maping between vMapG and vMapT
%   vMapG_MapP: Maping between vMapG and vMapP
% 
% Outputs: 
%   vGeneLocation: Matlab structure about the gene location information,
%                  generated in PlantClusterFinder function.
function [vGeneLocation] = f_remove_protein_transcript_info_from_genepositionfile(vGeneLocation, vGeneLocation_Att, vMapG, vMapT, vMapP, vMapG_MapT, vMapG_MapP, vVerbose)
    vNi = size(vGeneLocation.(char(vGeneLocation_Att(1,1))),1);
    for vi = 1:vNi
        if vVerbose >= 2
            fprintf('Processing transcript %i of %i\n', vi, vNi);
        end
        vCurrID = vGeneLocation.(char(vGeneLocation_Att(1,1)))(vi,1);
        vCurrID_T = strcmpi(vMapT,vCurrID);
        vCurrID_P = strcmpi(vMapP,vCurrID);
        vCurrID_G = sum(vMapG_MapT(:,vCurrID_T,1),2)>0 & sum(vMapG_MapP(:,vCurrID_P,1),2)>0;
        
        if sum(vCurrID_G) > 1
            error('Gene / Transcript / Protein location file (Gene position file) matches several entries to more than one gene.\nPlease clean up entries for %s in your input files.\n', char(vCurrID));
        elseif sum(vCurrID_G) == 1
            vGeneLocation.(char(vGeneLocation_Att(1,1)))(vi,1) = vMapG(vCurrID_G,1);
        else 
            vCurrID_G = strcmpi(vMapG,vCurrID);
            if sum(vCurrID_G) == 0
                error('Gene / Transcript / Protein location file (Gene position file) has an entry %s that is not covered by your gene ID conversion file.\nPlease add the entry.\n', char(vCurrID));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that collapses multiple entries (e.g. because it originally
% contained transcript or protein information) in the gene location data to
% one entry per gene. It takes the outermost end and start bp of all the
% protein or transcripts of a gene.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneLocation: Matlab structure about the gene location information,
%                  generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
% 
% Outputs: 
%   vUniqueChrom: A list with unique entries of all chromosome and scaffold
%                 names.
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
function [vUniqueChrom, vGeneLocation2] = f_consolidate_gene_info(vGeneLocation, vGeneLocation_Att, vVerbose)
    vUniqueChrom = unique(vGeneLocation.(char(vGeneLocation_Att(1,4)))); 
    vNi = size(vUniqueChrom,1);
    for vi = 1:vNi
        if vVerbose >= 2
            fprintf('Consolidate gene information from %i of %i\n', vi, vNi);
        end
        vTemp_IDs1 = strcmp(vGeneLocation.(char(vGeneLocation_Att(1,4))),vUniqueChrom(vi,1));
        vTemp_IDs1_find = find(vTemp_IDs1);
        [~,vTemp_IDs2] = sort(vGeneLocation.(char(vGeneLocation_Att(1,1)))(vTemp_IDs1,1));
        for vj = 1:size(vGeneLocation_Att,2)
            vTempGeneLocation.(char(vGeneLocation_Att(1,vj))) = vGeneLocation.(char(vGeneLocation_Att(1,vj)))(vTemp_IDs1_find(vTemp_IDs2),:);
        end
        vj = 2;
        while vj <= size(vTempGeneLocation.(char(vGeneLocation_Att(1,1))),1)
            if strcmp(vTempGeneLocation.(char(vGeneLocation_Att(1,1)))(vj,1),vTempGeneLocation.(char(vGeneLocation_Att(1,1)))(vj-1,1))
                vTempGeneLocation.(char(vGeneLocation_Att(1,2)))(vj-1,1) = cellstr(num2str(min([str2num(char(vTempGeneLocation.(char(vGeneLocation_Att(1,2)))(vj-1,1))),str2num(char(vTempGeneLocation.(char(vGeneLocation_Att(1,2)))(vj,1)))])));
                vTempGeneLocation.(char(vGeneLocation_Att(1,3)))(vj-1,1) = cellstr(num2str(max([str2num(char(vTempGeneLocation.(char(vGeneLocation_Att(1,3)))(vj-1,1))),str2num(char(vTempGeneLocation.(char(vGeneLocation_Att(1,3)))(vj,1)))])));
                for vk = 1:size(vGeneLocation_Att,2)
                    vTempGeneLocation.(char(vGeneLocation_Att(1,vk))) = vTempGeneLocation.(char(vGeneLocation_Att(1,vk)))([1:(vj-1),(vj+1):end],:);
                end
                vj = vj - 1;
            end
            vj = vj + 1;
        end
        if ~exist('vGeneLocation2', 'var')
            vGeneLocation2 = vTempGeneLocation;
        else
            for vk = 1:size(vGeneLocation_Att,2)
                vGeneLocation2.(char(vGeneLocation_Att(1,vk))) = [vGeneLocation2.(char(vGeneLocation_Att(1,vk)));vTempGeneLocation.(char(vGeneLocation_Att(1,vk)))];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that sorts all teh genes in the gene location structure fist
% accoring to the chromosome scaffold, then according to their start (start
% bp for plus strand, end pb for minus strand)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vUniqueChrom: A list with unique entries of all chromosome and scaffold
%                 names.
% 
% Outputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
function [vGeneLocation2] = f_sort_genes(vGeneLocation2, vGeneLocation_Att, vUniqueChrom, vVerbose)
    vGene_start_bp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(:,1)));
    vGene_end_bp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(:,1)));
    vGene_strand = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,5))),'1');
    vNi = size(vUniqueChrom,1);
    for vi = 1:vNi
        if vVerbose >= 2
            fprintf('Sorting chromosome or scaffold %i of %i\n', vi, vNi);
        end
        vTemp_IDs1 = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4))),vUniqueChrom(vi,1));
        vTemp_IDs1_find = find(vTemp_IDs1);
        vStart = vGene_start_bp;
        vStart(vGene_strand==0,1) = vGene_end_bp(vGene_strand==0,1);
        [~,vTemp_IDs2] = sort(vStart(vTemp_IDs1));
        for vk = 1:size(vGeneLocation_Att,2)
            vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vTemp_IDs1,:) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vTemp_IDs1_find(vTemp_IDs2),:);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function computes all stretches of DNA that are not encoding at least one
% gene model.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vUniqueChrom: A list with unique entries of all chromosome and scaffold
%                 names.
% 
% Outputs: 
%   vInterChrom: A matrix defining chromosome, start bp, end bp, gene
%                before and gene after a intergenic region.
function [vInterChrom] = f_find_intergenic_regions(vGeneLocation2, vGeneLocation_Att, vUniqueChrom, vVerbose)
    vInterChrom = zeros(size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),5);
    vIC = 1;
    vGene_start_bp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(:,1)));
    vGene_end_bp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(:,1)));
    vGene_Chrm = vGeneLocation2.(char(vGeneLocation_Att(1,4)));
    vNi = size(vUniqueChrom,1);
    for vi = 1:vNi
        vTemp_IDs1 = strcmp(vGene_Chrm,vUniqueChrom(vi,1));
        vChrStrand = true(max([vGene_start_bp(vTemp_IDs1,1);vGene_end_bp(vTemp_IDs1,1)]),1);
        vChrStrand(1:min([vGene_start_bp(vTemp_IDs1,1);vGene_end_bp(vTemp_IDs1,1)])-1,1) = 0;
        vTemp_IDs1_find = find(vTemp_IDs1);
        for vj = 1:size(vTemp_IDs1_find,1)
            vChrStrand(vGene_start_bp(vTemp_IDs1_find(vj,1),1):vGene_end_bp(vTemp_IDs1_find(vj,1),1)) = 0;
        end
        vInterGeneEnd_old = 0;
        while sum(vChrStrand) > 0
            if vVerbose >= 2
                fprintf('Chr/Scaff %i of %i, %i left, intergenic region %i\n', vi, vNi, sum(vChrStrand), vIC);
            end
            vInterGeneStart = find(vChrStrand,1,'first');
            vChrStrandNOT = ~vChrStrand;
            vChrStrandNOT(1:vInterGeneStart-1) = 0;
            if sum(vChrStrandNOT) == 0
                vInterGeneEnd = size(vChrStrandNOT,1);
            else
                vInterGeneEnd = find(vChrStrandNOT,1,'first') - 1;
            end
            
            %% I think here the correct gene is not identified (either the front or the latter gene....
            vInterChrom(vIC,:) = [vi, vInterGeneStart + vInterGeneEnd_old, vInterGeneEnd + vInterGeneEnd_old, ...
                find(vTemp_IDs1 & vGene_end_bp + 1 == vInterGeneStart + vInterGeneEnd_old,1,'last'), ...
                find(vTemp_IDs1 & vGene_start_bp - 1 == vInterGeneEnd + vInterGeneEnd_old,1,'first')];
            vIC = vIC + 1;
            vInterGeneEnd_old = vInterGeneEnd_old + vInterGeneEnd;
            vChrStrand = vChrStrand(vInterGeneEnd+1:end,1);
        end
    end
    vInterChrom = vInterChrom(1:(vIC-1),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that maps genes in the gene location file to the genes in the id
% conversion file.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vMapG: Unique list of gene IDs stored in the conversion file.
% 
% Outputs: 
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
function [vBiomG_G] = f_map_GeneLocation_to_conversion_file(vGeneLocation2, vGeneLocation_Att, vMapG, vVerbose)
%     vBiomG_G = false(size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),size(vMapG,1));
    vBiomG_G = logical(spalloc(size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),size(vMapG,1),size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1)*2));
    vNi = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
    for vi = 1:size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1)
        if vVerbose >= 2
            fprintf('Map gene location to id conversion file %i of %i\n', vi, vNi);
        end
        vBiomG_G(vi,:) = strcmp(vMapG,vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Generates a vector with length of all genes in the gene location
% information and looks up if this gene has a reaction in the pgdb.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%   vMapG_Gpgdb: A maping between Genes in the conversion file and genes
%                stored in the pgdb genes.dat file.
%   vG2R: Maping of all genes in the Genes_all to the Reactions in
%         Reactions_all
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   vBiomG_E: A vector describing which gene has at least one reaction
%             annotated in the pgdb (1, if no reaction 0)
function [vBiomG_E] = f_annotate_enyme_info_to_gene_location(vGeneLocation2, vGeneLocation_Att, vBiomG_G, vMapG_Gpgdb, vG2R, vVerbose)
    vBiomG_E = zeros(size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),1);
    vNi = size(vBiomG_E,1);
    for vi = 1:vNi
        if vVerbose >= 2
            fprintf('Annotate enzyme for gene %i of %i\n', vi, vNi);
        end
        vIDs = vBiomG_G(vi,:);
        vIDs = sum(vMapG_Gpgdb(vIDs,:),1)>0;
        if sum(sum(vG2R(vIDs,:),1))>0
            vBiomG_E(vi,1) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that clusters the genome. It lookes at enzymes that are at max
% spaced by vStepsize genes (checking for 0, 1, ..., vMinStepSize intervening
% genes. It continuous to increase the vStepsize until clusters contain
% more non-enzymes than enzymes. Clusters containing only one group of
% locally duplicated genes are turned down. Clusters with only one reaction
% aswell. Clusters need at least vEnzMinForClusters enzymes.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vMinStepSize: Minimal intervening gene size that should be computed
%                 (default is 5).
%   vMaxStepSize: Maximal intervening gene size that should be computed
%                 (default is 20).
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vBiomG_E: A vector describing which gene has at least one reaction
%             annotated in the pgdb (1, if no reaction 0)
%   vEnzMinForClusters: Defines the minimum number of enyzmes that need to
%                       be present in a cluster.
%   vKbRangeLD: Genes that are encoded more than this bp apart, are not
%               called local duplicates.
%   vGeneRangeLD: Number of genes that space two genes from each other.
%                 Genes that have more than this many intervening genes
%                 between themare not called local duplicates.
%   vG_MCL: Maping between the genes in the conversion-ID file and the
%           MCL-clustering file.
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%   vG2R: Maping of all genes in the Genes_all to the Reactions in
%         Reactions_all
%   vMapG_Gpgdb: A maping between Genes in the conversion file and genes
%                stored in the pgdb genes.dat file.
%   vCriterion: Criterion after which stepsize is chosen (default is 1).
%               1: as published in Schlapfer et al 2017 
%                  sum(enzymes in clusters) - sum(non-enzymes in clusters)
%                  > 0.
%               2: stepsize 3 (because most organisms take this one
%               3: sum(enzymes in clusters) - sum(non-enzymes in clusters,
%                  without hypo genes) -2*sum(hypogenes)
%               4: sum(clusters in stepsize n) > sum(clusters in stepsize
%                  n+1).
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs:
%   vClusters: Cluster Information (stepsize, chromosome id, gene start,
%              gene end, number of enzymes)
%   vClusterBackground: Density information about the genome (window size, 
%                       number of enzymes, Occurrence of windowsize with
%                       this many enzymes in the genome, occurence of
%                       this windowsize in the genome, Frequency of seeing 
%                       a windowsize of this size with at least as many
%                       enzymes in the genome)
%   vFinalStepsize: Stepsize that was chosen by the criterion (clusters
%                   contain in average more enzymes than non-enzymes).
%   vStepSize_Max: Maximal stepsize that was computed.
function [vClusters, vClusterBackground, vFinalStepsize, vStepSize_Max] = f_perform_clustering(vMinStepSize, vMaxStepSize, vGeneLocation2, vGeneLocation_Att, vBiomG_E, vEnzMinForClusters, vKbRangeLD, vGeneRangeLD, vG_MCL, vBiomG_G, vG2R, vMapG_Gpgdb, vCriterion, vVerbose)
    vContinueSearch = 1;
    vClusters = zeros(vMinStepSize*size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),5); % Stepsize, start, end, number of enzymes, number of hypos;
    vClusterBackground = zeros(vMinStepSize*size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),5); % Size, number of enzymes, count, number of total;
    vClusterWrite = 1;
    vClusterBackgroundWrite = 1;
    vStepSize = 0;
    vBiom_Hypo = strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,1))),'Hypo',4);
    while vStepSize <= vMinStepSize || vContinueSearch == 1
        % Get all clusters
        vi = 1;
        while vi <= size(vBiomG_E,1)
            if vVerbose >= 2
                fprintf('Potential Gene start %i of %i\n', vi, size(vBiomG_E,1));
            end
            vStart = vi - 1 + find(vBiomG_E(vi:end,1),1,'first');
            vEnd = vStart;
            vAbolishCluster = 0;
            if size(vStart,1) == 0
                vAbolishCluster = 1;
                vStart = size(vBiomG_E,1) + 1;
                vEnd = vStart;
            else
                vChrom_start = vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vStart,1);
                vBiomG_Chrom = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vChrom_start);
                
                vStopSearchEnd = 0;
                while vStopSearchEnd == 0
                    vEnds = vEnd+1;
                    vEnde = vEnd+1+vStepSize;
                    if vEnde > size(vBiomG_E,1)
                        vEnde = size(vBiomG_E,1);
                    end
                    if size(find(vBiomG_E(vEnds:vEnde,1) & vBiomG_Chrom(vEnds:vEnde,1),1,'last'),1)>0
                        vEnd = vEnds - 1 + find(vBiomG_E(vEnds:vEnde,1) & vBiomG_Chrom(vEnds:vEnde,1),1,'last');
                    else
                        vStopSearchEnd = 1;
                    end
                end
            end
            
            %Cluster needs at least 3 enzymes
            if vAbolishCluster == 0
                if sum(vBiomG_E(vStart:vEnd)) < vEnzMinForClusters
                    vAbolishCluster = 1;
                end
            end
            
            if vAbolishCluster == 0
                vLDisNOTaProblem = 0;
                % Get rid of clusters that only contain LD
                vEIDs = vStart+find(vBiomG_E(vStart:vEnd))-1;
%                 vEIDs_find = find(vEIDs);
                % Test if they are encoded within kb range
                vStart_kb = min([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vEIDs,1))); ...
                    str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vEIDs,1)))]);
                vEnd_kb = max([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vEIDs,1))); ...
                    str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vEIDs,1)))]);
                if vEnd_kb-vStart_kb+1 > vKbRangeLD
                    vLDisNOTaProblem = 1;
                end
                % Test if they are encoded witin gene range
                if vLDisNOTaProblem == 0
                    if sum(~strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vStart:vEnd,1),'Hypothetical_Gene_',18)) >= vGeneRangeLD+2
                        vLDisNOTaProblem = 1;
                    end
                end
                %Check if there are at least two local duplication groups
                if vLDisNOTaProblem == 0
                    vMCL = sum(vG_MCL(sum(vBiomG_G(vEIDs(1,1),:),1)>0,:),1)>0;
                    vj = 1;
                    vTwoGroups = 0;
                    while vj <= size(vEIDs,1) && vTwoGroups == 0
                        if sum(vMCL & sum(vG_MCL(sum(vBiomG_G(vEIDs(vj,1),:),1)>0,:),1)==0) > 0
                            vTwoGroups = 1;
                        end
                        vj = vj + 1;
                    end
                    if vTwoGroups == 1
                        vLDisNOTaProblem = 1;
                    end
                end
                if vLDisNOTaProblem == 0
                    vAbolishCluster = 1;
                end
            end
            % Get rid of clusters that only have one reaction set
            if vAbolishCluster == 0
                vRxns = sum(vG2R(sum(vMapG_Gpgdb(sum(vBiomG_G(vEIDs(1,1),:),1)>0,:),1)>0,:),1)>0;
                vj = 1;
                vTwoGroups = 0;
                while vj <= size(vEIDs,1) && vTwoGroups == 0
                    if sum(vRxns & sum(vG2R(sum(vMapG_Gpgdb(sum(vBiomG_G(vEIDs(vj,1),:),1)>0,:),1)>0,:),1)==0) > 0
                        vTwoGroups = 1;
                    end
                    vj = vj + 1;
                end
                if vTwoGroups == 0
                    vAbolishCluster = 1;
                end
            end
            if vAbolishCluster == 0
                vClusters(vClusterWrite,:) = [vStepSize,vEIDs(1,1),vEIDs(end,1),size(vEIDs,1),sum(vBiom_Hypo(vEIDs(1,1):vEIDs(end,1),1))];
                vClusterWrite = vClusterWrite + 1;
            end
            vi = vEnd + 1;
        end
        
        vC_IDs = vClusters(:,1) == vStepSize;
        vC_IDs(vClusterWrite:end,1) = 0;
        
        % Compute background of clusters
        % Get all sizes of clusters of this round
        [vClusterBackground, vClusterBackgroundWrite] = f_compute_background(vClusterBackground, vClusters, vClusterBackgroundWrite, vGeneLocation2, vGeneLocation_Att, vBiomG_E, vC_IDs);

        % Evaluate clusters for criterion (stepsize found?)
        if vVerbose >= 2
            if vCriterion == 1
                fprintf('Number of enzymes - number of non-enzymes in clustered genes: %i\n', sum(vClusters(vC_IDs,4))-(sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)-sum(vClusters(vC_IDs,4))));
            elseif vCriterion == 3
                fprintf('Number of enzymes - number of non-enzymes (no hypos) - 2* number of hypos in clustered genes: %i\n', sum(vClusters(vC_IDs,4))-(sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)-sum(vClusters(vC_IDs,4))-sum(vClusters(vC_IDs,5)))-2*sum(vClusters(vC_IDs,5)))
            elseif vCriterion == 4
                if vStepSize ~= 0
                    fprintf('Number of clustes more than last stepsize: %i\n', size(vClusters(vC_IDs,1)) - size(vClusters(:,1)==vStepSize-1));
                else
                    fprintf('Number of clustes Stepsize 0: %i\n', size(vClusters(vC_IDs,1)));
                end
            end
        end
        if vCriterion == 1
            if vStepSize == vMaxStepSize
                vContinueSearch = 0;
                vFinalStepsize = vStepSize;
                fprintf('Warning: maximal step size reached: %i. Consider increasing this variable.\n', vMaxStepSize)
            elseif sum(vC_IDs)==0
                fprintf('Warning: no clusters found for stepsize %i.\n', vStepSize);
            elseif sum(vClusters(vC_IDs,4))/sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)<0.5
                if vContinueSearch == 1
                    vFinalStepsize = vStepSize-1;
                end
                vContinueSearch = 0;
            end
        elseif vCriterion == 2
            if vStepSize == vMinStepSize
                vFinalStepsize = 3;
                vContinueSearch = 0;
            end
        elseif vCriterion == 3
            if vStepSize == vMaxStepSize
                vContinueSearch = 0;
                vFinalStepsize = vStepSize;
                fprintf('Warning: maximal step size reached: %i. Consider increasing this variable.\n', vMaxStepSize)
            elseif sum(vC_IDs)==0
                fprintf('Warning: no clusters found for stepsize %i.\n', vStepSize);
            elseif sum(vClusters(vC_IDs,4))-(sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)-sum(vClusters(vC_IDs,4))-sum(vClusters(vC_IDs,5)))-2*sum(vClusters(vC_IDs,5)) < 0
                if vContinueSearch == 1
                    vFinalStepsize = vStepSize-1;
                end
                vContinueSearch = 0;
            end
        elseif vCriterion == 4
            if vStepSize == vMaxStepSize
                vContinueSearch = 0;
                vFinalStepsize = vStepSize;
                fprintf('Warning: maximal step size reached: %i. Consider increasing this variable.\n', vMaxStepSize)
            elseif sum(vC_IDs)==0
                fprintf('Warning: no clusters found for stepsize %i.\n', vStepSize);
            elseif size(vClusters(vC_IDs,1)) - size(vClusters(:,1)==vStepSize-1) < 0
                if vContinueSearch == 1
                    vFinalStepsize = vStepSize-1;
                end
                vContinueSearch = 0;
            end
        end
        vStepSize = vStepSize + 1;
    end
    vClusters = vClusters(1:vClusterWrite-1,:);
    vClusterBackground = vClusterBackground(vClusterBackground(:,4) ~= 0,:);
    vStepSize_Max = vStepSize-1;
    fprintf('Stepsize chosen: %i\n', vFinalStepsize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that lookes at all clusters that were found by the function
% f_perform_clustering(), and computes the background information for the
% windowsizes of these clusters (e.g. a cluster has 4 enzymes and two non
% enzymes, then its size is 6. The function goes with a window size 6
% through the genome and checks how many windows exist with 0, 1, 2, ...,
% 5, 6 enzymes. It then calculates the probability to see a cluster with
% more or equal than 0 enzymes, more or equal than 1, ... more or equal
% than 6 enzymes. This is then used to define top n% clusters, e.g.
% clusters that are at least as dense as the top 5% of windowsizes are
% marked to be of higher interest.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vClusterBackground: Density information about the genome (window size, 
%                       number of enzymes, Occurrence of windowsize with
%                       this many enzymes in the genome, occurence of
%                       this windowsize in the genome, Frequency of seeing 
%                       a windowsize of this size with at least as many
%                       enzymes in the genome)
%   vClusters: Cluster Information (stepsize, chromosome id, gene start,
%              gene end, number of enzymes)
%   vClusterBackgroundWrite: Defines where the next entry should be written
%                            in vClusterBackground.
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vBiomG_E: A vector describing which gene has at least one reaction
%             annotated in the pgdb (1, if no reaction 0)
%   vC_IDs: Defines which Entries in vClusters are Clusters of the current
%           Stepsize that was calculated and should be checked if the
%           background of that size of cluster (number of genes) was
%           already calculated.
% 
% Outputs:
%   vClusterBackground: Density information about the genome (window size, 
%                       number of enzymes, Occurrence of windowsize with
%                       this many enzymes in the genome, occurence of
%                       this windowsize in the genome, Frequency of seeing 
%                       a windowsize of this size with at least as many
%                       enzymes in the genome)
%   vClusterBackgroundWrite: Defines where the next entry should be written
%                            in vClusterBackground.
function [vClusterBackground, vClusterBackgroundWrite] = f_compute_background(vClusterBackground, vClusters, vClusterBackgroundWrite, vGeneLocation2, vGeneLocation_Att, vBiomG_E, vC_IDs)
     vBackgroundToCheck = unique(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1);
     % Check which backgrounds have already been produced and compute new backgrounds
     vi = 1;
     while vi <= size(vBackgroundToCheck,1)
         if sum(vClusterBackground(1:vClusterBackgroundWrite-1,1) == vBackgroundToCheck(vi,1)) == 0
             vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),1) = vBackgroundToCheck(vi,1);
             vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),2) = 0:1:vBackgroundToCheck(vi,1);
             for vk = 1:size(vBiomG_E,1)-(vBackgroundToCheck(vi,1)-1)
                 if size(unique(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vk:(vk+vBackgroundToCheck(vi,1)-1),1)),1) == 1
                     vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),4) = vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),4) + 1;
                     vClusterBackground(vClusterBackgroundWrite+sum(vBiomG_E(vk:vk+vBackgroundToCheck(vi,1)-1,1)),3) = vClusterBackground(vClusterBackgroundWrite+sum(vBiomG_E(vk:vk+vBackgroundToCheck(vi,1)-1,1)),3) + 1;
                 end
             end
             for vk = vClusterBackgroundWrite+vBackgroundToCheck(vi,1):-1:vClusterBackgroundWrite
                 vClusterBackground(vk,5) = sum(vClusterBackground(vk:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),3))/vClusterBackground(vk,4);
             end
             vClusterBackgroundWrite = vClusterBackgroundWrite + vBackgroundToCheck(vi,1) + 1;
         end
         vi = vi + 1;
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that (given clusters and a background density information of a
% genome) calculates whether the clusters are within the top n% of enzymes 
% dense clusters.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   vClusters: Cluster Information (stepsize, chromosome id, gene start,
%              gene end, number of enzymes)
%   vClusterBackground: Density information about the genome (window size, 
%                       number of enzymes, Occurrence of windowsize with
%                       this many enzymes in the genome, occurence of
%                       this windowsize in the genome, Frequency of seeing 
%                       a windowsize of this size with at least as many
%                       enzymes in the genome)
%   vTopPercentClusters: Clusters are labeled as being in the top n% if
%                        they are within this top n% of enzyme dense
%                        clusters compared to the background (default is
%                        5).
%
% Outputs:
%   vClusters_Top: Vector with size Number of Clusters x 1. It is zero if
%                  the density of the cluster is above or the same as given
%                  by the input argument vTopPercentClusters*100.
function [vClusters_Top] = f_perform_cutoff(vClusters, vClusterBackground, vTopPercentClusters)
    vClusters_Top = false(size(vClusters,1),1);
    for vj = 1:size(vClusters,1)
        vSize = vClusters(vj,3)-vClusters(vj,2) + 1;
        vEnzymeSize = vClusters(vj,4);
        vTempBackground_ID = vClusterBackground(:,1) == vSize & vClusterBackground(:,2) == vEnzymeSize;
        if vClusterBackground(vTempBackground_ID,5) <= vTopPercentClusters/100
            vClusters_Top(vj,1) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function that writes the genefile of PlantClusterFinder: TSV format,
% first it includes the information that was stored in the gene location
% file that was provided to Plantclusterfinder (Gene IDs, Gene Start, Gene
% end, Chromosome information, strand). Then it attaches Enzyme (1 or 0) EC
% information and Reaction information. Then whether the Gene is considered
% locally duplicated and if so, which LD-clusters it's proteins were
% associated with. Then per stepsize that was investigated (intervening
% gene size) the ClusterID where the gene is incorporated and the
% informatino whether the cluster was found to be in the top n% of
% clsuters. This is followed by the information about the metabolic
% domains, and the signature / tailoring information.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%   vGeneOutputFile: Path and file name to the gene-Outputfile, make sure
%                    you have access to the file and folder. This is the
%                    first Results output file that will be generated.
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vFinalStepsize: Stepsize that was chosen by the criterion (clusters
%                   contain in average more enzymes than non-enzymes).
%   vStepSize_Max: Maximal stepsize that was computed.
%   vTopPercentClusters: Clusters are labeled as being in the top n% if
%                        they are within this top n% of enzyme dense
%                        clusters compared to the background (default is
%                        5).
%   vMD_to_annotate: A list of Metabolic domains that should be analyzed.
%                    For example {'Amines and Polyamines Metabolism'; ...
%                    'Amino Acids Metabolism'
%   vG2R: Maping of all genes in the Genes_all to the Reactions in
%         Reactions_all
%   vMapG_Gpgdb: A maping between Genes in the conversion file and genes
%                stored in the pgdb genes.dat file.
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%   vReactions_all: Matlab structure containing all the information
%                   (Attributes)that was searched for (described in
%                   vReactionsAttributes, default is {'UNIQUE-ID', ...
%                   'ENZYMATIC-REACTION', 'EC-NUMBER'};) in reactions.dat
%                   of the pgdb.
%   vReactions_all_Att: Attributes of vReactions_all, can be displayed by
%                       vReactions_all.(char(vReactions_all_Att(1,n)))
%   vBiomG_E: A vector describing which gene has at least one reaction
%             annotated in the pgdb (1, if no reaction 0)
%   vG_LD_ClustIDs: The MCL clusters of the gene (comma separated if more
%                   than 1 peptide is associated to a gene
%   vG_LD: A vector of length of the number of Genes in vGeneLocation2,
%          being 1 if the gene was considered to be locally duplicated (a
%          gene was within vKbRangeLD bp distance and not spaced more than
%          with vGeneRangeLD genes from another gene with the same MCL
%          clustering given by vG_MCL.
%   vClusters: Cluster Information (stepsize, chromosome id, gene start,
%              gene end, number of enzymes)
%   vClusters_Top: Vector with size Number of Clusters x 1. It is zero if
%                  the density of the cluster is above or the same as given
%                  by the input argument vTopPercentClusters*100.
%   vR_MD: Maping between all Reactions stored in the Metabolic domain file
%          and the metabolic domains that were identified to be of interest
%          in vMD_to_annotate.
%   vRpgdb_Rmd: A maping between reactions stored in the pgdb and the
%               reactions given in the metabolic domain file
%   vSTClasses: Is a list of signature and tailoring enzyme classes.
%   vR_STClasses: Mapping of Reactions in vReactions_all to the enzyme
%                 classes of vSTClasses.
%   vSigOrTail: Is a list of types of enzymes, e.g. signature, tailoring
%               and so on.
%   vR_SigOrTail: Mapping of Reactions in vReactions_all to the enzyme
%                 types of vSigOrTail.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
%   vOutputType: Defines which type of output is given: simple (default
%                setting) gives only top n% for the final stepsize. verbose
%                gives all computed stepsize until the final stepsize. old
%                gives at least 5 until the final stepsize.
%
% Outputs:
%   None (File is generated with the path of the first input argument)
function f_write_Gene_Output(vGeneOutputFile, vGeneLocation2, vGeneLocation_Att, vFinalStepsize, vStepSize_Max, vTopPercentClusters, vMD_to_annotate, ...
        vG2R, vMapG_Gpgdb, vBiomG_G, vReactions_all, vReactions_all_Att, vBiomG_E, ...
        vG_LD_ClustIDs, vG_LD, vClusters, vClusters_Top, vR_MD, vRpgdb_Rmd, ...
        vSTClasses, vR_STClasses, vSigOrTail, vR_SigOrTail, vVerbose, vOutputType)
    vFIO = fopen(vGeneOutputFile, 'w');
    
    %% Print Header
    %Gene Location
    fprintf(vFIO, '%s\t%s\t%s\t%s\t%s', vGeneLocation_Att{1,1:5});
    %Gene information
    fprintf(vFIO, '\tEnzyme\tEC\tRXN\tLocally Duplicated\tLD_ClusterIds');
    %Clustering information
    if strcmp(vOutputType,'old')
        for vi = 0:vStepSize_Max
            if vi == vFinalStepsize
                fprintf(vFIO, '\tCHOSEN Stepsize %i\tCHOSEN Stepsize %i, Top %i%%', vi, vi, vTopPercentClusters);
            else
                fprintf(vFIO, '\tStepsize %i\tStepsize %i, Top %i%%', vi, vi, vTopPercentClusters);
            end
        end
    elseif strcmp(vOutputType,'simple')
        fprintf(vFIO, '\tStepsize %i, Top %i%%', vFinalStepsize, vTopPercentClusters);
    elseif strcmp(vOutputType,'verbose')
        for vi = 0:vFinalStepsize
            fprintf(vFIO, '\tStepsize %i\tStepsize %i, Top %i%%', vi, vi, vTopPercentClusters);
        end
    else
        error('Outputtype %s not supported.',vOutputType);
    end
    %Metabolic Domain info
    for vi = 1:size(vMD_to_annotate,1)
        fprintf(vFIO, '\t%s', char(vMD_to_annotate(vi,1)));
    end
    %Signature Tailoring info
    fprintf(vFIO, '\tSignature Tailoring Enzyme Class\tSignature or Tailoring');
    fprintf(vFIO,'\n');
    
    %% Print Data
    vNvi = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
    for vi = 1:vNvi %Print all genes
        if vVerbose >= 2
            fprintf('Write gene %i of %i\n', vi, vNvi);
        end
        %Print Gene Location
        fprintf(vFIO, '%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1)));
        fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vi,1)));
        fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vi,1)));
        fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,6)))(vi,1)));
        fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,5)))(vi,1)));
        
        %Compose Gene information (reactions, enzymes, local duplications...)
        vRXN_IDs = sum(vG2R(sum(vMapG_Gpgdb(vBiomG_G(vi,:),:),1)>0,:),1)>0;
        vRXN_IDs_find = find(vRXN_IDs);
        vECstring = '';
        vRXNstring = '';
        if size(vRXN_IDs_find,2)>0
            vRXNstring = char(vReactions_all.UNIQUE_ID(vRXN_IDs_find(1,1),1));
            if ~isempty(vReactions_all.(char(vReactions_all_Att(1,3))){vRXN_IDs_find(1,1),1})
                vECstring = char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,1),1));
            end
            for vk = 2:size(vReactions_all.(char(vReactions_all_Att(1,3))),2)
                if ~isempty(vReactions_all.(char(vReactions_all_Att(1,3))){vRXN_IDs_find(1,1),vk})
                    vECstring = [vECstring ', ' char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,1),vk))];
                end
            end
            for vj = 2:size(vRXN_IDs_find,2)
                vRXNstring = [vRXNstring ', ' char(vReactions_all.UNIQUE_ID(vRXN_IDs_find(1,vj)))];
                vECstring = [vECstring ', ' char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,vj),1))];
                for vk = 2:size(vReactions_all.(char(vReactions_all_Att(1,3))),2)
                    if ~isempty(vReactions_all.(char(vReactions_all_Att(1,3))){vRXN_IDs_find(1,vj),vk})
                        vECstring = [vECstring ', ' char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,vj),vk))];
                    end
                end
            end
        end
        vOld = size(vRXNstring,2)+1;
        vNew = size(vRXNstring,2);
        while vOld > vNew
            vOld = size(vRXNstring,2);
            vRXNstring = regexprep(vRXNstring, ', , ', ', ');
            vNew = size(vRXNstring,2);
        end
        vOld = size(vECstring,2)+1;
        vNew = size(vECstring,2);
        while vOld > vNew
            vOld = size(vECstring,2);
            vECstring = regexprep(vECstring, ', , ', ', ');
            vNew = size(vECstring,2);
        end
        vRXNstring = regexprep(vRXNstring, ', $', '');
        vRXNstring = regexprep(vRXNstring, '^, ', '');
        vECstring = regexprep(vECstring, ', $', '');
        vECstring = regexprep(vECstring, '^, ', '');
        % Print gene information.
        fprintf(vFIO, '\t%i\t%s\t%s', vBiomG_E(vi,1), vECstring, vRXNstring);
        if isempty(vG_LD_ClustIDs{vi,1})
            vG_LD_ClustIDs(vi,1) = cellstr('');
        end
        fprintf(vFIO, '\t%i\t%s', vG_LD(vi,1), char(vG_LD_ClustIDs(vi,1)));
        % Print cluster information
        if strcmp(vOutputType,'old')
            for vj = 0:vStepSize_Max
                vClustIDs = (vClusters(:,1) == vj) & (vClusters(:,2) <= vi) & (vi <= vClusters(:,3));
                if sum(vClustIDs)>0
                    fprintf(vFIO, '\t%s\t%i', ['C' num2str(find(vClustIDs))], vClusters_Top(vClustIDs));
                else
                    fprintf(vFIO, '\t\t');
                end
            end
        elseif strcmp(vOutputType,'simple')
            vClustIDs = (vClusters(:,1) == vFinalStepsize) & (vClusters(:,2) <= vi) & (vi <= vClusters(:,3));
            if sum(vClustIDs)>0
                if vClusters_Top(vClustIDs) == 1
                    fprintf(vFIO, '\t%s', ['C' num2str(find(vClustIDs))]);
                else
                    fprintf(vFIO, '\t');
                end
            else
                fprintf(vFIO, '\t');
            end 
        elseif strcmp(vOutputType,'verbose')
            for vj = 0:vFinalStepsize
                vClustIDs = (vClusters(:,1) == vj) & (vClusters(:,2) <= vi) & (vi <= vClusters(:,3));
                if sum(vClustIDs)>0
                    fprintf(vFIO, '\t%s\t%i', ['C' num2str(find(vClustIDs))], vClusters_Top(vClustIDs));
                else
                    fprintf(vFIO, '\t\t');
                end
            end
        end
        vMD_IDs = sum(vR_MD(sum(vRpgdb_Rmd(vRXN_IDs,:),1)>0,:),1)>0;
        for vj = 1:size(vMD_IDs,2)
            fprintf(vFIO, '\t%i',vMD_IDs(1,vj));
        end
        vSTClassesString = '';
        vSTClasses_find = find(sum(vR_STClasses(vRXN_IDs,:),1)>0);
        if size(vSTClasses_find,2)>0
            vSTClassesString = char(vSTClasses(vSTClasses_find(1,1)));
            for vj = 2:size(vSTClasses_find,2)
                vSTClassesString = [vSTClassesString ', ' char(vSTClasses(vSTClasses_find(1,vj)))];
            end
        end
        vSigOrTailString = '';
        
        vSIGCol = strcmpi(vSigOrTail,'signature');
        vSIGTailCOls = sum(vR_SigOrTail(vRXN_IDs,:),1)>0;
        if sum(vSIGTailCOls(1,vSIGCol))>0
            vSigOrTailString = 'signature';
        elseif sum(vSIGTailCOls(1,:))>0
            vSigOrTailString = 'tailoring';
        end
        fprintf(vFIO, '\t%s\t%s', vSTClassesString, vSigOrTailString);
        fprintf(vFIO, '\n');
    end
    fclose(vFIO);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that writes the clusterfile of PlantClusterFinder: TSV format,
% first it includes the stepsize of the cluster (intervening genes), then
% cluster ID, chromosome/scaffold the cluster is located on, start bp, end
% bp, physical size, start gene ID, end gene ID, and all the gene IDs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%   vClusterOutputFile: Path and file name to the cluster-Outputfile, make
%                       sure you have access to the file and folder. This
%                       is the second Results output file that will be
%                       generated.
%   vClusters: Cluster Information (stepsize, chromosome id, gene start,
%              gene end, number of enzymes)
%   vClusters_Top: Vector with size Number of Clusters x 1. It is zero if
%                  the density of the cluster is above or the same as given
%                  by the input argument vTopPercentClusters*100.
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields of the Matlab structure
%                      vGeneLocation and vGeneLocation2
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
%   vFinalStepsize: Stepsize that was chosen by the criterion (clusters
%                   contain in average more enzymes than non-enzymes).
%   vOutputType: Defines which type of output is given: simple (default
%                setting) gives only top n% for the final stepsize. verbose
%                gives all computed stepsize until the final stepsize. old
%                gives at least 5 until the final stepsize.
% 
% Outputs:
%   None (File is generated with the path of the first input argument)
function f_write_Cluster_Output(vClusterOutputFile, vClusters, vClusters_Top, vGeneLocation2, vGeneLocation_Att, vVerbose, vFinalStepsize, vOutputType)
    %Check for output type:
    if ~strcmp(vOutputType,'old') && ~strcmp(vOutputType,'verbose') && ~strcmp(vOutputType,'simple')
        error('OutputType has to either be simple or verbose')
    end
    vFIO = fopen(vClusterOutputFile, 'w');
    
    %Print Header
    fprintf(vFIO, 'Stepsize\tClusterID\tTop n%% cluster\tChromosome\tStart bp\tEnd bp\tPhysical size bp\tNumber of genes\tNumber of Hypothetical genes\tnumber of Enzymes\tStart Gene\tEnd Gene\tGene IDs\n');
    
    %Print Data
    vNvi = size(vClusters,1);
    for vi = 1:vNvi
        if strcmp(vOutputType,'old') || ...
            (strcmp(vOutputType,'verbose') && vClusters(vi,1) <= vFinalStepsize) || ...
            (strcmp(vOutputType,'simple') && vClusters(vi,1) == vFinalStepsize && vClusters_Top(vi) == 1)
            if vVerbose >= 2
                fprintf('Write cluster %i of %i\n', vi, vNvi);
            end
            fprintf(vFIO, '%i\t%s\t%i\t%s', vClusters(vi,1), ['C' num2str(vi)], vClusters_Top(vi,1), char(vGeneLocation2.(char(vGeneLocation_Att(1,6)))(vClusters(vi,2),1)));
            vStartbp = min([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vClusters(vi,2):vClusters(vi,3),1)));str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vClusters(vi,2):vClusters(vi,3),1)))]);
            vEndbp = max([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vClusters(vi,2):vClusters(vi,3),1)));str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vClusters(vi,2):vClusters(vi,3),1)))]);
            fprintf(vFIO, '\t%i\t%i\t%i',vStartbp, vEndbp, vEndbp-vStartbp+1);
            fprintf(vFIO, '\t%i\t%i\t%i',vClusters(vi,3)-vClusters(vi,2)+1, sum(strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,2):vClusters(vi,3)),'Hypothetical_Gene_',18)), vClusters(vi,4));
            vClustersString = char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,2),1));
            for vj = vClusters(vi,2)+1:vClusters(vi,3)
                vClustersString = [vClustersString ', ' char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vj,1))];
            end
            fprintf(vFIO, '\t%s\t%s\t%s\n',char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,2),1)),char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,3),1)),vClustersString);
        end
    end
    fclose(vFIO);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that checks if the user has access to read (and write) the data
% that serve as in- and ouputfiles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: 
%   vPGDB_FlatFileFolder: Full path to the folder where the flatfiles of
%                         the species pgdb is stored.
%   vMD_reactions_File: TSV of the metabolic domains of reactions. Check if
%                       file is outdated (needs to contain all reactions
%                       that are present in the PGDB that have genes
%                       annotated to).
%   vProtein_sequence_FastaFile: Protein sequence file of the genome.
%                                  Should only contain Protein IDs as
%                                  header e.g. >protein1-ID
%   vGeneTranscriptProtein_mapping_File: TSV file with a header describing
%                                        gene, transcript, and protein
%                                        identifier and then for all the
%                                        identfiers the listings (gene-ID,
%                                        transcript-ID, protein-ID).
%   vGeneLocation_File: TSV file with gene ID, start bp, end bp, chromosome
%                       / scaffold, strand (encoded as 1 or -1). This file
%                       can be generated with a biomart or out of a gff3
%                       file.
%   vDNA_FastaFile: Fasta file of the masked genome nucleotide sequences.
%   vSignatureTailorFile: TSV file describing the classes of enzymes that
%                         should be classified as signature or tailoring
%                         enzymes to identify clusters containing such.
%   vGeneOutputFile: Path and file name to the gene-Outputfile, make sure
%                    you have access to the file and folder. This is the
%                    first Results output file that will be generated.
%   vClusterOutputFile: Path and file name to the cluster-Outputfile, make
%                       sure you have access to the file and folder. This
%                       is the second Results output file that will be
%                       generated.
% 
% Outputs: 
%   None (File checks only)
function f_check_files_for_read_write(vPGDB_FlatFileFolder, vMD_reactions_File, vProtein_sequence_FastaFile, vGeneTranscriptProtein_mapping_File, vGeneLocation_File, vDNA_FastaFile, vSignatureTailorFile, vGeneOutputFile, vClusterOutputFile)
    vFIO = fopen([vPGDB_FlatFileFolder 'genes.dat']);
    if vFIO < 0
        error('Can not read %sgenes.dat\n', vPGDB_FlatFileFolder);
    end
    fclose(vFIO);
    
    vFIO = fopen([vPGDB_FlatFileFolder 'proteins.dat']);
    if vFIO < 0
        error('Can not read %sproteins.dat\n', vPGDB_FlatFileFolder);
    end
    fclose(vFIO);
    
    vFIO = fopen([vPGDB_FlatFileFolder 'enzrxns.dat']);
    if vFIO < 0
        error('Can not read %senzrxns.dat\n', vPGDB_FlatFileFolder);
    end
    fclose(vFIO);
    
    vFIO = fopen([vPGDB_FlatFileFolder 'reactions.dat']);
    if vFIO < 0
        error('Can not read %sreactions.dat\n', vPGDB_FlatFileFolder);
    end
    fclose(vFIO);
    
    vFIO = fopen(vMD_reactions_File);
    if vFIO < 0
        error('Can not read %s\n', vMD_reactions_File);
    end
    fclose(vFIO);
    
    vFIO = fopen(vProtein_sequence_FastaFile);
    if vFIO < 0
        error('Can not read %s\n', vProtein_sequence_FastaFile);
    end
    fclose(vFIO);
    
    vFIO = fopen(vGeneTranscriptProtein_mapping_File);
    if vFIO < 0
        error('Can not read %s\n', vGeneTranscriptProtein_mapping_File);
    end
    fclose(vFIO);
    
    vFIO = fopen(vGeneLocation_File);
    if vFIO < 0
        error('Can not read %s\n', vGeneLocation_File);
    end
    fclose(vFIO);
    
    vFIO = fopen(vDNA_FastaFile);
    if vFIO < 0
        error('Can not read %s\n', vDNA_FastaFile);
    end
    fclose(vFIO);
    
    vFIO = fopen(vSignatureTailorFile);
    if vFIO < 0
        error('Can not read %s\n', vSignatureTailorFile);
    end
    fclose(vFIO);
    
    vFIO = fopen(vGeneOutputFile, 'w');
    if vFIO < 0
        error('Can not read / write %s\n', vGeneOutputFile);
    end
    fclose(vFIO);
    
    vFIO = fopen(vClusterOutputFile, 'w');
    if vFIO < 0
        error('Can not read / write %s\n', vClusterOutputFile);
    end
    fclose(vFIO);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that prints the help output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%   None
% 
% Outputs:
%   None (Screen output)
function f_display_help()
    vCellHelp = f_help_content();
    for vi = 1:size(vCellHelp,1)
        fprintf('%s\n', char(vCellHelp(vi,1)));
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that defines the help output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%   None
% 
% Outputs:
%   vCellHelp: A cell array matrix containing the text that should be
%              printed out by the function f_display_help()
function vCellHelp = f_help_content()
vCellHelp = {'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'; ...
             '% '; ...
             '% PlantClusterFinder 1.3'; ...
             '% PlantClusterFinder (PCF) detects metabolic gene clusters in a sequenced'; ...
             '% genome. It uses a gene location file provided by the user (see below) and'; ...
             '% a PGDB created with Pathway Tools as well as further information (see'; ...
             '% below) to identify enzyme-coding genes (metabolic genes) located together'; ...
             '% on a chromosome. Initially only continous stretches of metabolic genes '; ...
             '% lying directly next to each other are allowed. This condition is relaxed'; ...
             '% by iteratively increasing the intervening (non-metabolic) gene size by'; ...
             '% one. Several criteria to select for clusters are provided. In addition to'; ...
             '% this, clusters can be prevented from forming by a section of criteria.'; ...
             '% Details of PCF (version 1.0) can be found in PMID: 28228535.  '; ...
             '% '; ...
             '% The major differences between this version (1.3) and previous versions '; ...
             '% (1.0 and 1.2) are:'; ...
             '% '; ...
             '%   1) Physical breaks of the genome or sequencing gaps of unknown size are'; ...
             '%      typically encoded by stretches of Ns in the genome assembly fasta'; ...
             '%      file. Previously we inserted N hypothetical genes. This however '; ...
             '%      diluted the background of low quality genomes with non-enzymes, and '; ...
             '%      hence the likelyhood of a cluster to be classified as top x% of '; ...
             '%      enzyme dense regions was better than in a genome that had good '; ...
             '%      quality. In version 1.3 we identify these breaks and prevent '; ...
             '%      formation of a cluster over these gaps. '; ...
             '%   2) Any sequencing information that is missing is typically hard masked '; ...
             '%      with Ns. Previously, any intergenic region affected by at least one '; ...
             '%      N was evaluated for its length, and hypothetical genes were inserted '; ...
             '%      accordingly (See Schlapfer et al, PMID: 28228535). This led '; ...
             '%      sometimes to unrealistic prevention of detecting gene clusters. It'; ... 
             '%      is unlikely that missing information about a single nucleotide would '; ...
             '%      (if it would be known) lead to the finding of multiple gene models. '; ...
             '%      Thus here we changed the code to insert 2 hypothetical genes only if '; ...
             '%      a strech of unknown sequence is larger than nth percentile of gene '; ...
             '%      sequences (set to 5). We also provide the option to NOT insert any '; ...
             '%      hypothetical genes all together. Instead we by default we use '; ...
             '%      MaxSeqGapSize set to 100000 and MaxInterGeneDistByMedian set to 50 '; ...
             '%      resulting in similar clusters as in PCF version 1.0.'; ...
             '%   3) Large gene poor intergenic regions are present in genomes. In this '; ...
             '%      version we provide the option of several parameters to prevent '; ...
             '%      clusters from spanning such large gene poor regions.  '; ...
             '% '; ...
             '% Authors: Pascal Schlapfer, December 2017'; ...
             '%          Bo Xue, December 2017'; ...
             '% '; ...
             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'; ...
             '% '; ...
             '% Mandatory inputs: all files need full path, order does not matter.'; ...
             '%  '; ...
             '%     -pgdb fullpath: Give full path of PGDB flat file folder, where the'; ...
             '%                     flat files of the species pgdb is stored.'; ...
             '%     -rmdf fullpath: TSV file of the metabolic domains of reactions. Check'; ...
             '%                     if file is outdated (needs to contain all reactions'; ...
             '%                     that are present in the PGDB that have genes'; ...
             '%                     annotated to).'; ...
             '%     -md list: A list of Metabolic domains that should be analyzed. For'; ...
             '%               example {''Amines and Polyamines Metabolism''; ''Amino Acids Metabolism''}'; ...
             '%     -psf fullpath: Protein sequence file of the genome. Can only contain'; ...
             '%                    Protein IDs as header e.g. >protein1-ID'; ...
             '%     -gtpf fullpath: TSV file with a header describing gene, transcript,'; ...
             '%                     and protein identifier and then for all the'; ...
             '%                     identfiers the listings (gene-ID, transcript-ID,'; ...
             '%                     protein-ID). This is used to map information'; ...
             '%                     regarding transcripts (currently not used) and'; ...
             '%                     proteins (used for MCL clustering, potentially also'; ...
             '%                     used in PGDBs) to the genes (e.g. the gene location'; ...
             '%                     information).'; ...
             '%     -glof fullpath: TSV file with gene ID, start bp, end bp, chromosome /'; ...
             '%                     scaffold, strand (encoded as 1 or -1). This file can'; ...
             '%                     be generated with a biomart or out of a gff3 file.'; ...
             '%     -dnaf fullpath: Fasta file of the HARDMASKED genome nucleotide'; ...
             '%                     sequences.'; ...
             '%     -sitf fullpath: TSV file describing the classes of enzymes that'; ...
             '%                     should be classified as signature or tailoring'; ...
             '%                     enzymes to identify clusters containing such.'; ...
             '%     -gout fullpath: Path and file name to the gene-Outputfile, make sure'; ...
             '%                     you have access to the file and folder. This is the'; ...
             '%                     first results output file that will be generated.'; ...
             '%     -cout fullpath: Path and file name to the cluster-Outputfile, make'; ...
             '%                     sure you have access to the file and folder. This is'; ...
             '%                     the second Results output file that will be'; ...
             '%                     generated.'; ...
             '%  '; ...
             '% Optional Inputs:'; ...
             '%     ''GenesDat_GeneInfo'': needs to be followed by a list of Attributes in'; ...
             '%                          genes.dat. Defines where gene ID information'; ...
             '%                          should be searched for. Default is ''Unique-ID'','; ...
             '%                          ''Accession-1'', ''Accession-2'''; ...
             '%     ''MCLClusterFile'': needs to be followed by either 1 or 0. With a default'; ...
             '%                       of 0. 1 indicates that the input file is not a'; ...
             '%                       protein fasta file, but a precomputed MCL'; ...
             '%                       clustering file (precomputing can be usefull for'; ...
             '%                       speed). If such a file is used, then nothing but'; ...
             '%                       gene IDs and tabs and new lines are allowed in this'; ...
             '%                       file.'; ...
             '%     ''HypoGenePercentile'': needs to be followed by a number. Defines the'; ...
             '%                           length of sequencing gap that should be'; ...
             '%                           populated by two hypothetical genes. If this is'; ...
             '%                           10, then two hypothetical genes are introduced'; ...
             '%                           as soon as the sequencing gap exceeds the'; ...
             '%                           (100-10) lower percentile of the background'; ...
             '%                           size of all genes of the genome (default is 5).'; ...
             '%     ''MaxSeqGapSize'': needs to be followed by a number. Maximal masked'; ...
             '%                      nucleotide region (in bp) that a cluster is allowed'; ...
             '%                      to bridge (default is 100000). '; ...
             '%     ''MaxInterGeneDist'': needs to be followed by a number. Maximal'; ...
             '%                         intergenic distance ( in bp) that still can be'; ...
             '%                         crossed by a gene cluster. (default is -1, thus'; ...
             '%                         inactive).'; ...
             '%     ''MaxInterGeneDistByMedian'': needs to be followed by a number. Maximal'; ...
             '%                                 intergenic distance that still can be'; ...
             '%                                 crossed by a gene cluster, here defined'; ...
             '%                                 by this number times the median of the'; ...
             '%                                 gene sizes. (default is 50. Set it to -1'; ...
             '%                                 if it should be inactive).'; ...
             '%     ''PercentileForMaxInterGeneDist'': needs to be followed by a number.'; ...
             '%                                      Calculates the maximal intergenic'; ...
             '%                                      distance that still can be crossed'; ...
             '%                                      by a gene cluster. This number'; ...
             '%                                      determines the percentile to choose'; ...
             '%                                      from all the distances. If this is'; ...
             '%                                      99.9, then a cluster is allowed to'; ...
             '%                                      bridge if the largest intergenic'; ...
             '%                                      distance in the cluster is below'; ...
             '%                                      99.9% of the background size of all'; ...
             '%                                      intergenic distances of the genome'; ...
             '%                                      (default is -1, thus inactive).'; ...
             '%     ''SeqGapSizesChromBreak'': needs to be followed by a list of numbers.'; ...
             '%                              The masked nucleotides sequences of these'; ...
             '%                              specific lengths will prevent cluster'; ...
             '%                              formation (Default is an empty list: []).'; ...
             '%     ''PreScreenGaps'': needs to be followed by either 1 or zero. Defines if'; ...
             '%                      sequencing gaps should be prescreened to delete them'; ...
             '%                      if they are anyway too small to matter. (default is'; ...
             '%                      0)'; ...
             '%     ''OverwriteSeqGapInfo'': needs to be followed by either 1 or 0. Forces'; ...
             '%                            masked nucleotide analysis for finding'; ...
             '%                            sequencing gaps (default is 0, no'; ...
             '%                            enforcement).'; ...
             '%     ''OverwriteMCLClustering'': needs to be followed by either 1 or 0.'; ...
             '%                               Forces (if set to 1) that the MCL'; ...
             '%                               clustering is redone even if the result'; ...
             '%                               files are already present (default is 0).'; ...
             '%     ''Verbose'': needs to be followed by a number (0 is default, 1 gives'; ...
             '%                more screen-output, 2 gives exhaustive screen-output).'; ...
             '%     ''EnzMinForClusters'': needs to be followed by a number: Defines the'; ...
             '%                          minimum number of enyzmes that need to be'; ...
             '%                          present in a cluster (default is 3).'; ...
             '%     ''KbRangeLD'': needs to be followed by a number. Genes that are'; ...
             '%                  separated by more than this size (in bp) are not'; ...
             '%                  considered as local duplicates. (default is 100000)'; ...
             '%     ''GeneRangeLD'': needs to be followed by a number. Number of'; ...
             '%                    intervening genes that are allowed to separate two'; ...
             '%                    local duplicated genes. Two genes that are separated'; ...
             '%                    by more than this many intervening genes are not'; ...
             '%                    called local duplicates. (default is 10)'; ...
             '%     ''TopPercentClusters'': needs to be followed by a number. Clusters are'; ...
             '%                           labeled if they are within this top X% of'; ...
             '%                           enzyme dense clusters compared to the'; ...
             '%                           background. (default is 5)'; ...
             '%     ''MinStepSize'': needs to be followed by a number. Minimal intervening'; ...
             '%                    gene size that should be computed (default is 5).'; ...
             '%     ''MaxStepSize'': needs to be followed by a number. Maximal intervening'; ...
             '%                    gene size that should be computed (default is 20).'; ...
             '%     ''Criterion'': needs to be followed by a number. Criterion after which'; ...
             '%                  stepsize (intervening gene size) is chosen (default is'; ...
             '%                  2).'; ...
             '%         1: as published in Schlapfer et al 2017 (PMID:28228535)'; ...
             '%             sum(enzymes in clusters) - sum(non-enzymes in clusters) > 0.'; ...
             '%         2: stepsize 3 (because most organisms take this one)'; ...
             '%         3: sum(enzymes in clusters) - sum(non-enzymes in clusters, without hypo genes) -2*sum(hypogenes)'; ...
             '%         4: sum(clusters in stepsize n) > sum(clusters in stepsize n+1).'; ...
             '%     -help: Display help'; ...
             '%     -para: needs to be followed by a number: how many cpus can be used by'; ...
             '%            the algorithm (used for MCL clustering and Sequencing gap'; ...
             '%            search). The default is 1.'; ...
             '%     ''UnmaskedDNA'': needs to be followed by a number. Defines if plant'; ...
             '%                    Cluster finder should continue despite a non-masked'; ...
             '%                    genome sequence file was used (default is 0, not'; ...
             '%                    continue).'; ...
             '%     ''PGDBIdsToMap'': needs to be followed by a string. Can contain the'; ...
             '%                     Letters G and/or T and/or P. In that case (the pgdb'; ...
             '%                     is mapped to Gene (G) and/or Transcript (T) and/or'; ...
             '%                     Protein (P) Ids of the Gene conversion file. The'; ...
             '%                     default value is ''G''.'; ...
             '%     ''RunAfterPart'': Needs to be followed by a number. Sets the parts of'; ...
             '%                     the algorithm that should be recomputed. (default is'; ...
             '%                     0)'; ...
             '%     ''Tempsaves'': Defines if Temporary saves should be made. (default is'; ...
             '%                  0)'; ...
             '%     ''TempsavesOverwrite'': Defines if Temporary saves should be'; ...
             '%                           overwritten. (Default is 0).'; ...
             '%     ''RemoveNonProtLocations'': Defines if the protein-fasta file should be'; ...
             '%                               used to clean out the gene position file'; ...
             '%                               from sequences that do not show up on the'; ...
             '%                               protein fasta file, and thus do not have'; ...
             '%                               had a chance to be predicted to be an'; ...
             '%                               enzyme. (Default 0)'; ...
             '%     ''InsertHypos'': needs to be followed by 0 or 1 (can be extended with'; ...
             '%                    code). Defines if the genome should be populated by'; ...
             '%                    hypothetical genes in sequencing gaps that are not'; ...
             '%                    physical breaks. (Default 0)'; ...
             '%     ''HypoAmount'': needs to be followed by a number. Defines by how many'; ...
             '%                   hypothetical genes a sequencing gap that is not a'; ...
             '%                   physical break should be populated by artificial genes'; ...
             '%                   to punish such sequencing gaps from becoming clusters.'; ...
             '%                   (Default 0)'; ...
             '%     ''OutputType'': Defines the format of the Outputfiles: ''old'' is using'; ...
             '%                   the old format to report all clusters of all stepsizes,'; ...
             '%                   ''verbose'' defines that top percent clusters of all'; ...
             '%                   stepsizes are reported. ''simple'' just reports the top'; ...
             '%                   percent clusters of the chosen stepsize (default'; ...
             '%                   simple).'; ...
             '% '; ...
             '%  Outputs:'; ...
             '%     None (Two files are generated: vGeneOutputFile, containing all the'; ...
             '%          information about genes and vClusterOutputFile, containing all'; ...
             '%          the information about clusters)'; ...
             '% '; ...
             '%  USAGE:'; ...
             '%  '; ...
             '%     MATLAB (paths with windows, please be aware that certain parts need'; ...
             '%            to be run in linux!):'; ...
             '%         PlantClusterFinder(''-pgdb'', vPGDB_FlatFileFolder, ''-rmdf'', ...'; ...
             '%                            vMD_reactions_File, ''-md'', ...'; ...
             '%                            vMD_to_annotate, ''-psf'', ...'; ...
             '%                            vProtein_sequence_FastaFile, ''-gtpf'', ...'; ...
             '%                            vGeneTranscriptProtein_mapping_File, ...'; ...
             '%                            ''-glof'', vGeneLocation_File, ''-dnaf'', ...'; ...
             '%                            vDNA_FastaFile, ''-sitf'', ...'; ...
             '%                            vSignatureTailorFile, ''-gout'', ...'; ...
             '%                            vGeneOutputFile, ''-cout'', ...'; ...
             '%                            vClusterOutputFile, varargin)'; ...
             '%     '; ...
             '%     Example: (use absolute paths!)'; ...
             '%         PlantClusterFinder(''-pgdb'', ''[PlantClusterFinder]\csubellipsoidea\pgdb\csubellipsoideacyc\1.0\data\'', ...'; ...
             '%         ''-rmdf'', ''[PlantClusterFinder]\Inputs\ReactionMetabolicDomainClassification.txt'', ...'; ...
             '%         ''-md'', {''Amines and Polyamines Metabolism''; ''Amino Acids Metabolism''; ''Carbohydrates Metabolism''; ''Cofactors Metabolism''; ''Detoxification Metabolism''; ''Energy Metabolism''; ''Fatty Acids and Lipids Metabolism''; ''Hormones Metabolism''; ''Inorganic Nutrients Metabolism''; ''Nitrogen-Containing Compounds''; ''Nucleotides Metabolism''; ''Phenylpropanoid Derivatives''; ''Polyketides''; ''Primary-Specialized Interface Metabolism''; ''Redox Metabolism''; ''Specialized Metabolism''; ''Sugar Derivatives''; ''Terpenoids''}, ...'; ...
             '%         ''-psf'', ''[PlantClusterFinder]\csubellipsoidea\CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa'', ...'; ...
             '%         ''-gtpf'', ''[PlantClusterFinder]\csubellipsoidea\gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt'', ...'; ...
             '%         ''-glof'', ''[PlantClusterFinder]\csubellipsoidea\glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt'', ...'; ...
             '%         ''-dnaf'', ''[PlantClusterFinder]\csubellipsoidea\CsubellipsoideaC_169_227_v2.0.hardmasked.fa'', ...'; ...
             '%         ''-sitf'', ''[PlantClusterFinder]\Inputs\scaffold-tailoring-reactions-05082016.tab'', ...'; ...
             '%         ''-gout'', ''[PlantClusterFinder]\csubellipsoidea\csubellipsoidea_Gene_v1_3.txt'', ...'; ...
             '%         ''-cout'', ''[PlantClusterFinder]\csubellipsoidea\csubellipsoidea_Clust_v1_3.txt'', ...'; ...
             '%         ''SeqGapSizesChromBreak'', [10000], ''PGDBIdsToMap'', ''GTP'');'; ...
             '%         % Note, PGDBIdsToMap is only needed here because the pgdb'; ...
             '%         % contains protein Ids in place of gene identifiers.'; ...
             '%         '; ...
             '%     Shell, standalone (linux, but you can compile a windows version, see'; ...
             '%                       how to compile new version file):'; ...
             '%         To use standalone, download matlab runtime from matlab website.'; ...
             '%         Make sure, you download the v91. Then replace'; ...
             '%         "/share/apps/MATLAB/MATLAB_Runtime/v91" int the example with your'; ...
             '%         path to your runtime.'; ...
             '%         '; ...
             '%         USE ABSOLUTE PATHS!'; ...
             '% '; ...
             '%         ./run_PlantClusterFinder.sh /share/apps/MATLAB/MATLAB_Runtime/v91 -pgdb "./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/" -rmdf "./Inputs/ReactionMetabolicDomainClassification.txt" -md "{''Amines and Polyamines Metabolism''; ''Amino Acids Metabolism''; ''Carbohydrates Metabolism''; ''Cofactors Metabolism''; ''Detoxification Metabolism''; ''Energy Metabolism''; ''Fatty Acids and Lipids Metabolism''; ''Hormones Metabolism''; ''Inorganic Nutrients Metabolism''; ''Nitrogen-Containing Compounds''; ''Nucleotides Metabolism''; ''Phenylpropanoid Derivatives''; ''Polyketides''; ''Primary-Specialized Interface Metabolism''; ''Redox Metabolism''; ''Specialized Metabolism''; ''Sugar Derivatives''; ''Terpenoids''}" -psf "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa" -gtpf "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt" -glof "./csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt" -dnaf "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa" -sitf "./Inputs/scaffold-tailoring-reactions-05082016.tab" -gout "./csubellipsoidea/ARAGene1_3_memex.txt" -cout "./csubellipsoidea/ARAClust1_3_memex.txt" SeqGapSizesChromBreak ''[10000]'' PGDBIdsToMap GTP'; ...
             '%         % Note, PGDBIdsToMap is only needed here because the pgdb'; ...
             '%         % contains protein Ids in place of gene identifiers. '; ...
             '% '};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that reads in TSV files with headers. The headers are used as
% fields in matlab structures. These field names are not allowed to have
% certain special characters. Follwoing Characters are currently handeled:
% ' ' (space), ')', '(', '#', '/', '-', '.'. It is possible that others
% exist that need to either be avoided as part of header or the code needs
% to be adapted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: 
%   vFile: Full path to TSV file with header.
%   vUnmaskedDNA: Option to ignore the error if you use unmasked DNA files.
% Optional Inputs:
%   'vSplitter': Use an alternative sign than tab ('\t') to split text into
%                columns.
% 
% Outputs: 
%   vOutput: Matlabstructure with all the content of the TSV file.
%            Different Columns are stored in different fields of the matlab
%            structure.
%   vHeader: Cell array of the fields of vOutput:
function [vOutput, vHeader] = f_extract_results_with_header(vFile,varargin)
    vSplitter = '\t';
    if nargin > 1
        for vi = 1:2:size(varargin,2)
            if strcmp(varargin{vi},'vSplitter')
                vSplitter = varargin{vi+1};
            end
        end
    end
    vFIO = fopen(vFile);
    vLine = fgetl(vFIO);
    vSplit = regexp(vLine,vSplitter,'split');
    vHeader = vSplit;
    for vi = 1:size(vHeader,2)
        vHeader(1,vi) = cellstr(strtrim(char(vHeader(1,vi))));
        vHeader(1,vi) = regexprep(vHeader(1,vi),' ','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'\)','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'\(','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'\(','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'#','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'\/','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'\-','_');
        vHeader(1,vi) = regexprep(vHeader(1,vi),'\.','_');
        
        if sum(strcmp(vHeader(1,1:vi-1),vHeader(1,vi)))>0
            vj = 1;
            while sum(strcmp(vHeader(1,1:vi-1),[char(vHeader(1,vi)) '_' num2str(vj)]))>0
                vj = vj + 1;
            end
            vHeader(1,vi) = cellstr([char(vHeader(1,vi)) '_' num2str(vj)]);
        end
        if strcmp(vHeader(1,vi),'')
            vHeader(1,vi) = cellstr(['Attribute_' num2str(vi)]);
        end
        vOutput.(char(vHeader(1,vi))) = cell(1,1);
    end
        
    vLine = fgetl(vFIO);
    vW = 0;
    while ischar(vLine)
        vSplit = regexp(vLine,vSplitter,'split');
        for vi = 1:size(vSplit,2)
            if exist('vSpacer','var')
                vOutput.(char(vHeader(1,vi)))(vW + 1,1) = regexprep(regexprep(vSplit(1,vi),['^' vSpacer],''),[vSpacer '$' ],'');
            else
                vOutput.(char(vHeader(1,vi)))(vW + 1,1) = vSplit(1,vi);
            end
        end
        vLine = fgetl(vFIO);
        vW = vW + 1;
    end
    fclose(vFIO);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that analyzes the masked DNA sequences of a _GAPOutput file, and
% counts the occurence of every size of gap.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: 
%   vFileIn: Full path to the _GAPOutput file that whould be analyzed
%            (produced within PlantClusterFinder.
%   vUnmaskedDNA: Option to ignore the error if you use unmasked DNA files.
%   vVerbose: Tells the code on how much of output should be printed out to
%             the screen. (0 = nothing, 1 = little, 2 = more, ...)
% 
% Outputs: 
%   None (File with list of sequencing gaps numbers is generated with
%   ending _count)
function f_analyze_PlantClusterGapFile(vFileIn, vUnmaskedDNA, vVerbose)
    vInit = 7000000;
    [vFIN,vFINERR] = fopen(vFileIn);
    if vFIN == -1
        error('Tried to open %s file, error: %s\n', vFileIn, vFINERR);
    end
    vLine = fgetl(vFIN);
    if ~ischar(vLine)
        error('%s is empty.\n',vFileIn);
    end
    vCounter = zeros(vInit,2);
    vCounter_write = 0;
    vL = 1;
    while ischar(vLine)
        if vVerbose >= 2
            fprintf('Analyze Sequencing gap file line %i\n', vL)
        end
        vSplit = regexp(vLine,' ','split');
        if strcmp(vSplit(1,2),'N')
            if vCounter_write+1 > vInit
                vCounter = [vCounter;zeros(10000,2)];
                vInit = vInit + 10000;
            end
            if sum(vCounter(:,1) == str2num(char(vSplit(1,4)))) == 0
                vCounter(vCounter_write+1,:) = [str2num(char(vSplit(1,4))),1];
                vCounter_write = vCounter_write + 1;
            else
                vPlace = vCounter(:,1) == str2num(char(vSplit(1,4)));
                vCounter(vPlace,2) = vCounter(vPlace,2) + 1;
            end
        end
        vLine = fgetl(vFIN);
        vL = vL + 1;
    end
    fclose(vFIN);
    vCounter = vCounter(1:vCounter_write,:);
    
    if ~exist('vCounter','var') && vUnmaskedDNA == 0
        error('DNA fasta File does not have masked DNA.\nIf you want to ignore this error, use Option ''UnmaskedDNA 1''\n');
    elseif ~exist('vCounter','var')
        warning('DNA fasta File does not have masked DNA.\n You used Option ''UnmaskedDNA 1'' with the intention to use Unmasked DNA file.\n');
    end
    
    [vFOUT,vFOUTERR] = fopen([vFileIn '_count.txt'],'w');
    if vFOUT == -1
        error('Tried to generate %s_count.txt file, error: %s\n', vFileIn, vFOUTERR);
    end
    fprintf(vFOUT, 'Gapsize in bp\tNumber of Gaps\n');
    if exist('vCounter','var')
        [~, vSortIds] = sort(vCounter(:,2),'descend');
        vCounter = vCounter(vSortIds,:);
        vi = 1;
        while vi <= size(vCounter,1)
            vIDs = vCounter(:,2) == vCounter(vi,2);
            vCounter2 = vCounter(vIDs,:);
            [~, vSortIds2] = sort(vCounter2(:,1));
            vCounter2 = vCounter2(vSortIds2,:);
            vCounter(vIDs,:) = vCounter2;
            vi = vi + sum(vIDs);
        end
        for vi = 1:size(vCounter,1)
            fprintf(vFOUT, '%i\t%i\n',vCounter(vi,1),vCounter(vi,2));
        end
    end
    fclose(vFOUT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that reads in a pgdb flat file (e.g. genes.dat, proteins.dat,
% enzrxns.dat, reactions.dat) saving Attributes given by vAttrib.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: 
%   vFile: PGDB Flat file that should be read in.
%   vAttrib: Attributes of the PGDB Flat file that should be saved (all
%            others are discarded.
% 
% Outputs: 
%   vOutput: Matlab structure with the contents of each attribute as a
%            field of the matlab structure.
%   vAttrib_Head: Fields of vOutput as cell arrays.
function [vOutput,vAttrib_Head] = f_read_in_PGDB_Flat(vFile,vAttrib)
    vFIO = fopen(vFile);
    vAttrib_Head = vAttrib;
    for vi = 1:size(vAttrib,2)
        vAttrib_Head(1,vi) = regexprep(vAttrib_Head(1,vi),'\-','_');
        vOutput.(char(vAttrib_Head(1,vi))) = cell(1,1);
    end
    vLine = fgetl(vFIO);
    vW = 0;
    while ischar(vLine)
        for vi = 1:size(vAttrib,2)
            if strncmp(vLine,vAttrib(1,vi),size(vAttrib{1,vi},2))
                if vi == 1
                    for vk = 2:size(vAttrib,2)
                        if size(vOutput.(char(vAttrib_Head(1,1))),1) ~= size(vOutput.(char(vAttrib_Head(1,vk))),1)
                            vOutput.(char(vAttrib_Head(1,vk)))(vW,vAttrib_nr(vk)+1) = cellstr('');
                        end
                    end
                    vAttrib_nr = zeros(1,size(vAttrib,2));
                    vW = vW + 1;
                end
                vOutput.(char(vAttrib_Head(1,vi)))(vW,vAttrib_nr(vi)+1) = cellstr(vLine(1,size(vAttrib_Head{1,vi},2)+4:end));
                vAttrib_nr(vi) = vAttrib_nr(vi)+1;
            end
        end
        vLine = fgetl(vFIO);
    end
    fclose(vFIO);
    for vk = 2:size(vAttrib,2)
        if size(vOutput.(char(vAttrib_Head(1,1))),1) ~= size(vOutput.(char(vAttrib_Head(1,vk))),1)
            vOutput.(char(vAttrib_Head(1,vk)))(vW,vAttrib_nr(vk)+1) = cellstr('');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that uses the protein fasta file to search for genes that encode
% proteins and then deletes all other entries in the gene location
% variable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vGeneLocation_Att: Describing the fields attributes of the Matlab
%                      structure vGeneLocation and vGeneLocation2.
%   vMapP: Unique list of protein IDs stored in the conversion file.
%   vMapG_MapP: Maping between vMapG and vMapP
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
%   vProtein_IDs: A cell array of all protein-IDs that are stored in the
%                 protein fasta file. only active if
%                 vRemoveNonProtLocations is 1. Else the array has no rows.
% 
% Outputs: 
%   vGeneLocation2: Matlab structure about the gene location information,
%                   generated in PlantClusterFinder function.
%   vBiomG_G: Maping between the Gene location file and the gene-IDs in the
%             conversion file.
function [vGeneLocation2, vBiomG_G] = f_Remove_non_protein_fasta_gene_locations(vGeneLocation2, vGeneLocation_Att, vMapP, vMapG_MapP, vBiomG_G, vProtein_IDs)
    vGeneLoc_in_Prot_F = false(size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),1);
    for vi = 1:size(vProtein_IDs,1)
        vGeneLoc_in_Prot_F(sum(vBiomG_G(:,sum(vMapG_MapP(:,strcmp(vMapP, vProtein_IDs(vi,1))),2)>0),2)>0) = 1;
    end
    for vi = 1:size(vGeneLocation_Att,2)
        vGeneLocation2.(char(vGeneLocation_Att(1,vi))) = vGeneLocation2.(char(vGeneLocation_Att(1,vi)))(vGeneLoc_in_Prot_F,:);
    end
    vBiomG_G = vBiomG_G(vGeneLoc_in_Prot_F,:);
end

