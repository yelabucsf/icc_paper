##script to run cellranger
## set to use m5a4xlarge

while getopts s:l:o: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
        l) library=${OPTARG};;
        o) out=${OPTARG};;
    esac
done


echo $sample
echo $library
echo $out

ref="/data/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"

cd $out

/home/ubuntu/tools/cellranger-arc-2.0.2/cellranger-arc count --id=${sample} \
                       --reference=${ref} \
                       --libraries=${library} \
                       --localcores=32 \
                       --localmem=128
