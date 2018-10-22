#!bash

usage()
{
cat << EOF
eg. bash $0 -i -r utr -g *.fa -v *.vcf -t t1 -s s1

Options:
    -i       STR           input bam eg. *.bam
    -d       STR           known SNP directory eg. /native/database_NGS/hg19/gatk_bundle
    -g       STR           reference genome eg. /native/database_NGS/hg19/reference/ucsc.hg19.fasta
    -t       INT           thread
    -r       STR           target bed region eg. /nextcode/nfs_test/users/pengfeizh/rarecode_local_sentieon/BED/SureSelect_Human_Exon_V5_merge150.bed
    -o       STR           output directory
    -m       STR           directory of genomic map for shapeit eg.~/bin/shapeit_v2/genetic_map_b37
    -c       STR           chrlist file eg. /titan2/shen_qing/sentieon_test/chr.list
    -h       help
EOF
    exit 0
}
[ $1 ] || usage

set -x

while getopts "hi:d:g:t:r:o:m:c:" OPTION
do
    case $OPTION in
        h) usage;;
        i) input=$OPTARG;;
        d) knownSNP_dir=$OPTARG;;
        g) ref=$OPTARG;;
        t) thread=$OPTARG;;
        r) targetbed=$OPTARG;;
        o) output_dir=$OPTARG;;
        m) gmap=$OPTARG;;
        c) chrlist=$OPTARG;;
        ?) usage;;
    esac
	done
shift $((OPTIND - 1))

input=`readlink -f $input`
knownSNP_dir=`readlink -f $knownSNP_dir`
ref=`readlink -f $ref`
targetbed=`readlink -f $targetbed`
chrlist=`readlink -f $chrlist`
gmap=`readlink -f $gmap`
output_dir=`readlink -f $output_dir`


## Ali Cloud related
export PYTHONPATH=$PYTHONPATH:/nextcode/nfs_test/users/dai_dandan/bin/pythonmodel/lib/python2.7/site-packages
#source ~/dx-toolkit/environment
#bcs login cn-shanghai LTAIP4LvBjKfolva f9wxOqlcGwAlIvNTYfH5qA4Gv8v15ubcs

## Sentieon related
export SENTIEON_LICENSE=nn-02.nextcode.pharmatechs.com:8990
export Sentieon_PATH=/nextcode/sge_software/sentieon/currentVersion/bin/sentieon

#knownSNP_dir=/native/database_NGS/hg19/gatk_bundle
#ref=/native/database_NGS/hg19/reference/ucsc.hg19.fasta
#input=/titan2/shen_qing/sentieon_test/WGC053164U_dedup_realign_recal.bam
#thread=10
#targetbed=/nextcode/nfs_test/users/pengfeizh/rarecode_local_sentieon/BED/SureSelect_Human_Exon_V5_merge150.bed
#output_dir=/titan2/shen_qing/sentieon_test
#gmap=~/bin/shapeit_v2/genetic_map_b37
#chrlist=/titan2/shen_qing/sentieon_test/chr.list
#qualcal=


db_1000G_phase1_indel=${knownSNP_dir}/1000G_phase1.indels.hg19.vcf
db_mills=${knownSNP_dir}/Mills_and_1000G_gold_standard.indels.hg19.vcf
dbsnp=${knownSNP_dir}/dbsnp_138.hg19.vcf
db_1000G_omni=${knownSNP_dir}/1000G_omni2.5.hg19.vcf
db_1000G_phase1_snp=${knownSNP_dir}/1000G_phase1.snps.high_confidence.hg19.vcf
hapmap=${knownSNP_dir}/hapmap_3.3.hg19.vcf

	
mkdir tmp
#bam index
for i in `cat $input`
do
name=`basename $i|cut -d '.' -f 1`
java -Djava.io.tmpdir=tmp/ -jar /nextcode/sge_software/picard/currentVersion/picard.jar \
					  BuildBamIndex I="$i"

#qualCal table
$Sentieon_PATH driver -i $i -r $ref -t $thread --interval $targetbed --interval_padding 150 \
--algo QualCal -k $dbsnp -k $db_1000G_phase1_indel -k $db_mills "$output_dir"/"$name".qualCal 2>>"$output_dir"/"$name".qualCal.log
qualcal="$output_dir"/"$name".qualCal

#Haplotyper
$Sentieon_PATH driver -r $ref -i $i -q $qualcal -t $thread --interval $targetbed --interval_padding 150 \
--algo Haplotyper --emit_mode gvcf -d $dbsnp --call_conf 30 --emit_conf 10 "$output_dir"/"$name".gvcf.gz 2>>"$output_dir"/"$name".gvcf.log
done

cd $output_dir
gvcflist=`ls *.gvcf.gz|paste -s -d ' '|perl -ne s'/ / -v /g;print;' `
echo "this is joint gvcf list"
echo $gvcflist

#joint gvcf
$Sentieon_PATH driver -r $ref -t $thread --interval $targetbed --interval_padding 150 \
--algo GVCFtyper -d $dbsnp -v $gvcflist --call_conf 30 --emit_conf 10 output.vcf.gz

mkdir joint_vcf && cd joint_vcf
ln -s ../output.vcf.gz
ln -s ../output.vcf.gz.tbi

#VQSR
$Sentieon_PATH driver -r $ref -t $thread --algo VarCal --vcf output.vcf.gz --var_type SNP \
--resource $hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0 \
--resource $db_1000G_omni --resource_param omni,known=false,training=true,truth=true,prior=12.0 \
--resource $db_1000G_phase1_snp --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 \
--resource $dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 \
--max_gaussians 8 --tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
--annotation QD --annotation DP --annotation MQRankSum --annotation ReadPosRankSum --annotation FS \
--tranches_file output_snp.tranches output_snp.recal

$Sentieon_PATH driver -r $ref -t $thread --algo ApplyVarCal --vcf output.vcf.gz \
--var_type SNP --recal output_snp.recal --tranches_file output_snp.tranches snp_recal_indel_raw.vcf.gz

$Sentieon_PATH driver -r $ref -t $thread --algo VarCal --vcf snp_recal_indel_raw.vcf.gz --var_type INDEL \
--resource $db_mills --resource_param mills,known=true,training=true,truth=true,prior=12.0 \
--max_gaussians 4 --tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
--annotation DP --annotation FS --annotation MQRankSum --annotation ReadPosRankSum \
--tranches_file output_indel.tranches output_indel.recal 

$Sentieon_PATH driver -r $ref -t $thread --algo ApplyVarCal --vcf snp_recal_indel_raw.vcf.gz --var_type INDEL \
--recal output_indel.recal --tranches_file output_indel.tranches recalibrated.vcf.gz

rm -f snp_recal_indel_raw.vcf.gz
rm -f output_snp* output_indel*

#remove bi-allelic?
zcat recalibrated.vcf.gz | perl -lane 'print if $F[4]!~/,/;'|gzip > recalibrated.no_bi_allelic.vcf.gz
ln -s recalibrated.no_bi_allelic.vcf.gz final.vcf.gz


#shapeit phased 
mkdir shapeit && cd shapeit
ln -s ../final.vcf.gz

for i in `cat $chrlist`
do
	vcftools --gzvcf final.vcf.gz --chr "$i" --recode --out output."$i"
done

mkdir shapeit_log 
for i in `cat $chrlist`
do
shapeit --thread $thread --input-vcf output."$i".recode.vcf -O output.phased."$i".vcf \
--input-map $gmap/genetic_map_"$i"_combined_b37.txt -L shapeit_log/shapeit_"$i".log >shapeit_log/"$i".log 2>shapeit_log/"$i".err
done

#bash ~/scripts/FR_script/vcf2tmp.sh . chr10 4000
#bash ~/scripts/FR_script/tmp2gor_totalvcf_3_sq.sh chr10.tmp1.gor prefix