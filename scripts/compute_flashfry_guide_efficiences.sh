java -Xmx4g -jar /iblm/netapp/home/karthik/FlashFry/FlashFry-assembly-1.15.jar \
 score \
 --input /iblm/netapp/home/karthik/FlashFry/Gasperini2019/guide_sequences.fasta.off_targets \
 --output /iblm/netapp/home/karthik/FlashFry/Gasperini2019/guide_sequences.fasta.off_targets.scores \
 --scoringMetrics doench2016cfd \
 --database /iblm/netapp/home/karthik/FlashFry/databases/hg19
