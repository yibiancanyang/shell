#!/bin/bash
#TransDecoder make transdecoder.pep
TransDecoder.LongOrfs -t target_transcripts.fasta
blastp -query transdecoder_dir/longest_orfs.pep  -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
hmmscan --cpu 8 --domtblout pfam.domtblout /path/to/Pfam-A.hmm transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t yellow_mRNA.fa --retain_pfam_hits pfam.domblout --retain_blastp_hits blastp.outfmt6
#trinotate
blastx -query yellow_mRNA.fa -db uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query transdecoder.pep -db uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm transdecoder.pep > pfam.log
signalp -f short -n signalp.out transdecoder.pep

