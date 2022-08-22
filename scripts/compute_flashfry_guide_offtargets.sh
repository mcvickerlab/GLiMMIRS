java -Xmx4g -jar /iblm/netapp/home/karthik/FlashFry/FlashFry-assembly-1.15.jar \
 discover \
 --database /iblm/netapp/home/karthik/gasperini_project/flashfry/hg19 \
 --fasta /iblm/netapp/home/karthik/FlashFry/Gasperini2019/guide_sequences.fasta \
 --output /iblm/netapp/home/karthik/FlashFry/Gasperini2019/guide_sequences.fasta.off_targets \
 --flankingSequence 0
