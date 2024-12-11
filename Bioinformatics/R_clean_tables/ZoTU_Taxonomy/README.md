## Scripts

- [Martineau_AlexARMS_bioinformatics_cleantable_bothseqruns.R](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/scripts/Martineau_AlexARMS_bioinformatics_cleantable_bothseqruns.R) : script performing merging of two sequencing runs outputs, sample rarefaction, filtering of BLASTn resulst, merging of BLASTn outputs from two different data bases and taxonomic assignment of ZOTUs.


## Outputs

- [both_seqruns_merged_abundance.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/both_seqruns_merged_abundance.csv) : raw abundance table of both sequencing runs merged
- [both_seqruns_merged_abundance_rarefied_nosg.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/both_seqruns_merged_abundance_rarefied_nosg.csv) : abundance table of both sequencing runs merged with singletons removed and rarefaction of samples performed
- [both_seqruns_merged_abundance_rarefied_nosg_noseq](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/both_seqruns_merged_abundance_rarefied_nosg_noseq.csv) : same as both_seqruns_merged_abundance_rarefied_nosg.csv without sequences information
- [both_seqruns_rarefied_nosg.fasta](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/both_seqruns_rarefied_nosg.fasta) : ZOTu fasta file of rarefied abundance matrix of both sequencing runs
- [bothseqruns_histogram_raw_read_abundance.png](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/bothseqruns_histogram_raw_read_abundance.png) : histogram of raw read abundance for each sample for both sequencing runs outputs merged
- [sponge_both_seqruns_raw_blastn_output.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/sponge_both_seqruns_raw_blastn_output.csv) : raw BLASTn outputs of both_seqruns_rarefied_nosg.fasta against sponge data base
- [sponge_both_seqruns_sequences_not_retained.fasta](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/sponge_both_seqruns_sequences_not_retained.fasta) : fasta file containing sequences that did not match sponges in the custom data base, based on parameters
- [sponge_nt_abundance_table.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/sponge_both_seqruns_taxonomy_table.csv) : taxonomic assignments of filtered BLASTn results from sponge data base (sponge_both_seqruns_raw_blastn_output.csv)
- [nt_both_seqruns_raw_blastn_output.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/nt_both_seqruns_raw_blastn_output.csv) : raw BLASTn outputs of sponge_both_seqruns_sequences_not_retained.fasta against GenBank nt data base
- [nt_both_seqruns_taxonomy_table.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/nt_both_seqruns_taxonomy_table.csv) : taxonomic assignments of filtered BLASTn results from nt data base (nt_both_seqruns_raw_blastn_output.csv)
- [sponge_nt_abundance_table.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/sponge_nt_abundance_table.csv) : final abundance table for both sequencing runs merging taxonomic results from both data bases
- [sponge_nt_taxonomy_table.csv](https://github.com/martineauG/Alex_ARMS/blob/main/Bioinformatics/R_clean_tables/ZoTU_Taxonomy/outputs/sponge_nt_taxonomy_table.csv) : final taxonomy table for both sequencing runs merging taxonomic results from both data bases
