if [ $# -ge 1 ]; then 
  INPUT=$1
else 
  INPUT="Runlist/input_test.txt"
fi

if [ ! -f ${INPUT} ]; then eosmount ~/eos; fi

outdir=`basename ${INPUT}`
OUTPUT="output/${outdir}"

if [ -d ${OUTPUT} ]; then
  rm -rf ${OUTPUT}
fi

if [ ! -d `dirname ${OUTPUT}` ]; then
  mkdir `dirname ${OUTPUT}` -p
fi

echo "runFourLep -outdir $OUTPUT -in $INPUT"
#runFourLep -outdir $OUTPUT -in $INPUT > log.txt
runFourLep -outdir $OUTPUT -in $INPUT
