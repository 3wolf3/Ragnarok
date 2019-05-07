#!/bin/sh
###source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc62-opt/setup.sh
export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_9/biglib/slc6_amd64_gcc630:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_9/lib/slc6_amd64_gcc630:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_9/external/slc6_amd64_gcc630/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/llvm/4.0.1/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/lib:$LD_LIBRARY_PATH
cd ${_CONDOR_SCRATCH_DIR}
xrdcp -s root://cmseos.fnal.gov//store/user/3wolf3/condor_tarballs/$4 .
tar -xf $4
rm $4
./main ${1} ${2} ${3}
#root -l -q -b 'EventReadOut_v22.C("'${1}'","'${2}'","'${3}'")'
#root -l -q -b 'trigEffStudy_2017data.C("'${1}'","'${2}'","'${3}'")'
### Now that the run is over, there is one or more root files created
echo "List all root files = "
ls *.root
### Change this to just *.root maybe, or $1/*.root....
echo "List all files"
ls 
echo "*******************************************"
OUTDIR=root://cmseos.fnal.gov//store/user/3wolf3/May6/Data/$1/
echo "xrdcp output for condor"
#for FILE in $1/*/*.root
for FILE in *.root
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  echo "${FILE}" 
  echo "${OUTDIR}"
 xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done





