perl ../module_yellow/adddotone.pl yellow.fa out.fa
sed -i 's/\r//g' out.fa 
~/software_database/TransDecoder-3.0.1/TransDecoder.LongOrfs -t out.fa
cd out.fa.transdecoder_dir/
blastp -query longest_orfs.pep  -db ../../uniprot_sprot.pep  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6
hmmscan --cpu 24 --domtblout pfam.domtblout ../../Pfam-A.hmm longest_orfs.pep > pfam.log
cd ..
~/software_database/TransDecoder-3.0.1/TransDecoder.Predict -t out.fa --retain_pfam_hits ./out.fa.transdecoder_dir/pfam.domtblout  --retain_blastp_hits out.fa.transdecoder_dir/blastp.outfmt6 
blastx -query out.fa -db ../uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query out.fa.transdecoder.pep -db ../uniprot_sprot.pep -num_threads 24 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 24 --domtblout TrinotatePFAM.out ../Pfam-A.hmm out.fa.transdecoder.pep > pfam.log
~/software_database/signalp-4.1/signalp -f short -n signalp.out out.fa.transdecoder.pep
~/software_database/tmhmm-2.0c/bin/tmhmm --short < out.fa.transdecoder.pep > tmhmm.out
~/software_database/Trinotate-3.0.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome out.fa --path_to_rnammer ~/software_database/rnammer-1.2.src/rnammer
~/software_database/Trinotate-3.0.1/Trinotate ../Trinotate_V3.sqlite init --gene_trans_map yellow_gene_trans.txt --transcript_fasta out.fa --transdecoder_pep longest_orfs.pep
~/software_database/Trinotate-3.0.1/Trinotate ../Trinotate_V3.sqlite LOAD_swissprot_blastp blastp.outfmt6 
~/software_database/Trinotate-3.0.1/Trinotate ../Trinotate_V3.sqlite LOAD_swissprot_blastx blastx.outfmt6 
~/software_database/Trinotate-3.0.1/Trinotate ../Trinotate_V3.sqlite LOAD_pfam TrinotatePFAM.out
~/software_database/Trinotate-3.0.1/Trinotate ../Trinotate_V3.sqlite LOAD_tmhmm tmhmm.out
~/software_database/Trinotate-3.0.1/Trinotate ../Trinotate_V3.sqlite LOAD_signalp signalp.out 
