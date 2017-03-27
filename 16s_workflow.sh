#fastq to fasta
awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' 1.fastq > 1.fasta
#test map.tsv
validate_mapping_file.py -o vmf-map/ -m map.tsv
#OTU分析
pick_open_reference_otus.py -o otus/ -i 1.fasta -p ../uc_fast_params.txt
'##################
uc_fast_params.txt 内容为
pick_otus:enable_rev_strand_match True
###################'

#counts计数
biom summarize-table -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom

#########################################################
'
Num samples: 12
Num observations: 13033
Total count: 392643
Table density (fraction of non-zero values): 0.244

Counts/sample summary:
 Min: 27310.0
 Max: 38910.0
 Median: 32806.000
 Mean: 32720.250
 Std. dev.: 3838.757
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
 waterqb: 27310.0
 waterqa: 27948.0
 waterqc: 28595.0
 waterxc: 29069.0
 INqa: 30944.0
 Inxa: 31656.0
 waterxa: 33956.0
 Inxc: 34397.0
 waterxb: 35669.0
 INqc: 36577.0
 Inxb: 37612.0
 INqb: 38910.0
'
############################################################
#根据各个sample的counts数确定差异分析参数，取各样本最小值
core_diversity_analyses.py -o cdout/ -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom -m map.tsv -t otus/rep_set.tre -e 27310
根据map.tsv各栏参数作组间差异分析
core_diversity_analyses.py -o cdout/ --recover_from_failure -c "SampleType" -i otus/otu_table_mc2_w_tax_no_pynast_failures.biom -m map.tsv -t otus/rep_set.tre -e 27310