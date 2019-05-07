# Ragnarok

#RECIPE:

source /cvmfs/cms.cern.ch/cmsset_default.csh

setenv SCRAM_ARCH slc6_amd64_gcc630

cmsrel CMSSW_9_4_9

cd CMSSW_9_4_9/src

cmsenv

git cms-init

git clone https://github.com/btannenw/MiniAOD.git -b CMSSW_940

git clone https://github.com/btannenw/ttH-LeptonPlusJets.git -b 94x

git cms-merge-topic cms-egamma:EgammaPostRecoTools_940

git cms-merge-topic guitargeek:EgammaID_9_4_X

git cms-merge-topic cms-met:METFixEE2017_949_v2

git cms-merge-topic yrath:deterministicSeeds

scram b -j 8


#Common Classifier Code


git clone https://gitlab.cern.ch/ttH/CommonClassifier.git TTH/CommonClassifier

source TTH/CommonClassifier/setup/install_mem.sh

source TTH/CommonClassifier/setup/install_recoLikelihood.sh

scram b -j 8


#Ragnarok


git clone https://github.com/3wolf3/Ragnarok.git








Future Sanity
g++ -o main -Wl,-rpath . EventReadOut.C `root-config --cflags --glibs`  ./libTTHCommonClassifier.so
