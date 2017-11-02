#!/bin/bash

vsearch -fastx_filter total.fastq -fastaout total.1.filtered.fasta -fastq_maxee 1 -threads 8

#usearch -derep_fulllength total.fasta -fastaout total.2.derep.fasta -sizeout -uc total.2a.derep.uc -threads 8
vsearch -derep_fulllength total.1.filtered.fasta -output total.2.derep.fasta
# # dereplicate seqs (headers become: >pool426.677C.RLP(28).bmt_21;size=38122)
# # use vsearch instead of usearch (32-bit) for large files.
#
usearch -sortbysize total.2.derep.fasta -minsize 2 -fastaout total.3.sizesorted.fasta
# # sort by size (by looking at size=xxxx) and remove singletons
#
# #usearch -cluster_otus total.3.sizesorted.fasta -otus total.4.repset-raw.fasta -uparseout total.4a.otu-summary.txt -relabel OTU_ -sizein -sizeout
usearch -cluster_otus total.3.sizesorted.fasta -otus total.4.repset-raw.fasta -uparseout total.4a.otu-summary.txt -relabel OTU_ -sizein -sizeout
# # automatically removes de novo chimeras here (usearch7 and above)
# # provides sequence rep set and otu classification summary
# # needs to have ';size=xxxx;' from derep step
#
usearch -uchime2_ref total.4.repset-raw.fasta -db ~/Databases/gold.fa -strand plus -nonchimeras total.5.repset.fasta -threads 8
# #fasta_number.py total.5.repset-nochimeras.fasta OTU_ > total.6.repset.fasta
# # label OTUs using UPARSE python script: new headers: OTU_###
#
#usearch -usearch_global total.fasta -db total.5.repset.fasta -strand plus -id 0.97 -uc total.5a.otu.uc -threads 8
vsearch -usearch_global total.1.filtered.fasta -db total.5.repset.fasta -strand plus -id 0.97 -uc total.5a.otu.uc -threads 8
# #map the _original_ quality filtered reads back to OTUs
# #creates a table showing otu membership (otu.map.uc)
# #non-matches are reported as 'N' at beginning of line (instead of H)
# #Some reads may not match any OTU for these reasons:
# #(1) the read is chimeric,
# #(2) the read has more than 3% errors,
# #(3) the read has a singleton sequence so was discarded.
#
python ~/scripts/uc2otutab_mod.py total.5a.otu.uc > total.6.otu-table.txt
# # make OTU table. I modified the function 'GetSampleID' in the script 'uc2otutab.py' and renamed the script 'uc2otutab_mod.py':
# # The modified function is: function is:
# # def GetSampleId(Label):
# #    SampleID = Label.split()[0].split('_')[0]
# #    return SampleID
#
biom convert --table-type="OTU table" -i total.6.otu-table.txt -o total.7.otu-table.biom --to-json
# # convert to biom format
#
# # classification, ref=gg99
# mothur "#classify.seqs(fasta=total.5.repset.fasta, taxonomy=/Users/matthewwipperman/Databases/Gg_13_8_99.taxonomy/gg_13_8_99.gg.tax, template=/Users/matthewwipperman/Databases/Gg_13_8_99.taxonomy/gg_13_8_99.fasta, iters=1000, processors=8)"
assign_taxonomy.py -i total.5.repset.fasta -o tax --id_to_taxonomy_fp /Users/matthewwipperman/Databases/Gg_13_8_99.taxonomy/gg_13_8_99.gg.tax --reference_seqs_fp /Users/matthewwipperman/Databases/Gg_13_8_99.taxonomy/gg_13_8_99.fasta --assignment_method mothur --confidence 0
blastn.py total.5.repset.fasta -db refseq_rna
# # assign taxonomy
#
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp tax/total.5.repset_tax_assignments.txt -i total.7.otu-table.biom -o total.8.otu-tax.biom --output-as-json
# #add taxonomy metadata to biom
#
align_seqs.py -i total.5.repset.fasta -o aligned --template_fp ~/Databases/core_set_aligned.fasta.imputed
# ##align_seqs.py -i total.5.repset.fasta -o aligned -m muscle
filter_alignment.py -i aligned/total.5.repset_aligned.fasta -o aligned
make_phylogeny.py -i aligned/total.5.repset_aligned_pfiltered.fasta -o total.10.tree --tree_method fasttree --log_fp aligned/make_phylo.log
# # create phylogenetic tree
#
# #alpha_diversity.py -i total.8.otu-tax.biom -o alpha_diversity.txt -t total.10.tree -m simpson_reciprocal,goods_coverage,observed_species,shannon,simpson,simpson_e,PD_whole_tree,chao1
# #beta_diversity.py -i total.8.otu-tax.biom -o beta_diversity -t total.10.tree -m euclidean,binary_euclidean,bray_curtis,bray_curtis_faith,bray_curtis_magurran,manhattan,pearson,spearman_approx,unifrac,unifrac_g,unifrac_g_full_tree,unweighted_unifrac,unweighted_unifrac_full_tree,weighted_normalized_unifrac,weighted_unifrac
# #calculate alpha and beta diversity
#
#send text when done
echo 'Now you can make pretty pictures in R :)' | mail -s 'TBRU run complete!' 'matthew.wipperman@gmail.com'
#
