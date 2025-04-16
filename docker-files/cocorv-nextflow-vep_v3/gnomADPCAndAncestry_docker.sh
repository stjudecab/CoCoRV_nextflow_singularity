# export MAMBA_ROOT_PREFIX='/research/rgs01/home/clusterHome/wchen1/micromamba'
# __mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
# if [ $? -eq 0 ]; then
#     eval "$__mamba_setup"
# else
#     alias micromamba="$MAMBA_EXE"  # Fallback on help from mamba activate
# fi
# unset __mamba_setup
# eval "$(micromamba shell hook --shell bash)"
# micromamba activate cocorv-env-gnomad

main() {
  CoCoRVFolder=$1
  loadingPath=$2
  rfModelPath=$3
  VCFForPrediction=$4
  build=$5
  PCPopulationFile=$6
  threshold=$7
  populationFile=$8

  python ${CoCoRVFolder}/utilities/gnomADPCAndAncestry.py ${loadingPath} ${rfModelPath} ${VCFForPrediction} ${build} ${PCPopulationFile} ${threshold}
  zcat ${PCPopulationFile} | cut -f2 --complement  > ${populationFile}
}

main "$@"