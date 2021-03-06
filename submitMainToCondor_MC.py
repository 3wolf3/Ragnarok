# /usr/bin/python

#Author: Ben Tannenwald
#Date: March 28, 2018
#Purpose: Script to submit condor jobs for all files in data/mc sample

import os,sys, argparse

# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--outputDir", help="output directory for processed histograms and roofiles")
#TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
parser.add_argument("--inputTXTfile", help=".txt file containing list of input files for a given sample")
#inputTXTFile.txt
parser.add_argument("--isMC", help="flag whether sample is Monte Carlo (true) or data (false)")
parser.add_argument("--trigSF", help="flag whether run is for trigger SF study (true) or analysis plots (false)", action='store_false')
args = parser.parse_args()

if (len(vars(args)) != 4): # 4 --> four: one for each options
    os.system('python submitMainToCondor.py -h')
    quit()

# ** A. Test output directory existence and create if DNE
if(args.outputDir is None):
    print "#### Need to specify output directory using --outputDir <desired output directory> ####\nEXITING"
    quit()
else:
    if ( not os.path.exists(args.outputDir) ):
        print "Specified input file ({0}) DNE.\nCREATING NOW".format(args.inputTXTfile)
        os.system("mkdir {0}".format(args.outputDir))

    if ( not os.path.exists( (args.outputDir + '/condor_logs/') ) ):
        os.system("mkdir {0}".format( (args.outputDir + '/condor_logs/')) )
    if ( not os.path.exists( (args.outputDir + '/condor_err/') ) ):
        os.system("mkdir {0}".format( (args.outputDir + '/condor_err/')) )
    if ( not os.path.exists( (args.outputDir + '/condor_out/') ) ):
        os.system("mkdir {0}".format( (args.outputDir + '/condor_out/')) )

    print '-- Setting outputDir = {0}'.format(args.outputDir)

# ** B. Test input .txt file and exit if DNE
if(args.inputTXTfile is None):
    print "#### Need to specify input .txt file using --inputTXTfile <address to .txt file> ####\nEXITING"
    quit()
else:
    if ( not os.path.exists(args.inputTXTfile) ):
        print "#### Specified input file ({0}) DNE ####.\nEXITING".format(args.inputTXTfile)
        quit()
    else:
        print '-- Setting inputTXTfile = {0}'.format(args.inputTXTfile)

# ** C. Test isMC flag and exit if not sensible
if(args.isMC is None):
    print "#### Need to set isMC flag using --isMC <true/false> ####\nEXITING"
    quit()
else:
    if( not(args.isMC == "true" or args.isMC == "True" or args.isMC == "false" or args.isMC == "False") ):
        print "#### Please use true/True/false/False when setting --isMC <option>. Supplied value ({0}) does not match ####\nEXITING".format(args.isMC)
    else:
        print '-- Setting isMC = {0}'.format(args.isMC)

# ** D. Test trigSF flag and exit if not sensible
if(args.trigSF is False):
    print "#### Running in analysis plotting mode ####\n"
else:
    print "#### Running in triggerSF mode ####\n"


# ** E. Exit if no grid proxy
#if ( not os.path.exists(os.path.expandvars("$X509_USER_PROXY")) ):
#    print "#### No GRID PROXY detected. Please do voms-proxy-init -voms cms before submitting Condor jobs ####.\nEXITING"
#    quit()


# *** 1. Create .tar of directory and store in personal EOS
print "##########     Tarring workdir     ##########"
#tarball_name = "{0}.tar.gz".format(args.outputDir)
tarball_name = "MasterBall.tar.gz".format(args.outputDir)
#os.system("tar -cvzf {0} ./ --exclude 'plots*' --exclude 'runOver*'  --exclude '.git' --exclude 'test*' --exclude 'submitOneFile_' --exclude '*.tar.gz' --exclude 'ttbar*' --exclude '*-18' --exclude '*2018' --exclude 'MET*' --exclude 'single*' --exclude 'pass*' --exclude 'quick*' --exclude 'Output*.root'  --exclude 'oldFilelists' --exclude 'fastPlots*' --exclude 'syncExercise*.csv'   --exclude 'jetHT*'".format(tarball_name))
#if ( not os.path.exists("/eos/uscms/store/user/3wolf3/{0}/".format(args.outputDir)) ):
#    os.system("mkdir /eos/uscms/store/user/3wolf3/{0}/".format(args.outputDir))
#os.system("xrdcp {0} root://cmseos.fnal.gov//store/user/3wolf3/{1}/".format(tarball_name, args.outputDir))
#os.system("xrdcp {0} root://cmseos.fnal.gov//store/user/benjtann/{0}/{1}".format(tarball_name, args.outputDir))
#os.system("rm {0}".format(tarball_name))

# *** 2. Create temporary .pdl file for condor submission
print "\n##########     Submitting Condor jobs     ##########\n"
txtfile = open(args.inputTXTfile, 'r')

for line in txtfile:
    number = (line.rsplit('_',1))[1].split('.root')[0] # get number of file
    infile = line.split('\n')[0]
    jdl_filename = "submitOneFile_{0}_{1}.jdl".format(args.outputDir, number)

    os.system("touch {0}".format(jdl_filename))
    os.system("echo universe = vanilla > {0}".format(jdl_filename))
    os.system("echo Executable = runOneMain_MC.csh >> {0}".format(jdl_filename))
    os.system("echo Should_Transfer_Files = YES >> {0}".format(jdl_filename))
    os.system("echo WhenToTransferOutput = ON_EXIT >> {0}".format(jdl_filename))
    os.system("echo Transfer_Input_Files = runOneMain_MC.csh, {0} >> {1}".format(tarball_name, jdl_filename))
    #os.system("echo notify_user = benjamin.tannenwald@CERN.CH >> {0}".format(jdl_filename))
    #os.system("notify_user = benjtann@FNAL.GOV >> {0}".format(jdl_filename))
    os.system("echo Output = {0}/condor_out/outfile_{1}.out  >> {2}".format(args.outputDir, number, jdl_filename))
    os.system("echo Error = {0}/condor_err/outfile_{1}.err >> {2}".format(args.outputDir, number, jdl_filename))
    os.system("echo Log = {0}/condor_logs/outfile_{1}.log >> {2}".format(args.outputDir, number, jdl_filename))
    #os.system("echo x509userproxy = ${{X509_USER_PROXY}} >> {0}".format(jdl_filename))
    os.system("echo Arguments = {0} {1} {2} {3} >> {4}".format(args.outputDir, args.isMC, infile, tarball_name, jdl_filename))
    os.system("echo Queue 1 >> {0}".format(jdl_filename))   

    
    os.system("condor_submit {0}".format(jdl_filename))


# *** 3. Cleanup submission directory
print "\n##########     Cleanup submission directory     ##########\n"
#os.system("rm submitOneFile_{0}_*.jdl".format(args.outputDir))
#os.system("mv {0} {1}/".format(tarball_name, args.outputDir) )





