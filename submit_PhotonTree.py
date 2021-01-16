import os,glob

from SLURMWorker.SLURMWorker import SLURMWorker

# ______________________________________________________________________ ||
input_file_pattern  = "/cmsuf/data/store/user/t2/users/klo/Delphes/ALP_HToZaTo2l2g_M1/2020-06-02/2020-06-02_ALP_HToZaTo2l2g_M1_13TeV_*.root"
job_name            = "2020-09-22_PhotonTreeProducer"
out_dir             = "/cmsuf/data/store/user/t2/users/klo/Delphes/ALP_HToZaTo2l2g_M1/2020-09-22/PhotonTreeProducer/"

cmssw_dir           = "/ufrc/avery/kinho.lo/Delphes/CMSSW_10_0_5/src/"
delphes_dir         = os.getcwd()

# ______________________________________________________________________ ||
input_file_list = [f for f in glob.glob(input_file_pattern)]
input_file_list.sort()
n_file = len(input_file_list)
for ijob in range(0,n_file):
    each_job_name = job_name+"_"+str(ijob)
    out_file_name = each_job_name+".root"
    commands = "\n".join([
        "cd "+cmssw_dir,
        "scramv1 runtime -sh",
        "cd "+delphes_dir,
        "root -b -q \'PhotonTreeProducer.C(\"{inputFile}\",\"{outputFile}\")\'".format(inputFile=input_file_list[ijob],outputFile=os.path.join(out_dir,out_file_name)),
        ],)
    script_file_name = os.path.join(out_dir,each_job_name+".cfg")
    worker = SLURMWorker()
    worker.make_sbatch_script(
        script_file_name,
        job_name,
        "kin.ho.lo@cern.ch",
        "1",
        "1gb",
        "03:00:00",
        os.path.join(out_dir,each_job_name),
        commands,
        )
    worker.sbatch_submit(script_file_name)
