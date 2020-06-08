import os,glob

from SLURMWorker.SLURMWorker import SLURMWorker

# ______________________________________________________________________ ||
input_file_pattern  = "/cmsuf/data/store/user/klo/ALP_HToZaTo2l2g/GEN-SIM/ALP_HToZaTo2l2g_M1/94X_mc2017_2020May/200510_211100/000*/*.root"
job_name            = "2020-06-02_ALP_HToZaTo2l2g_M1_13TeV"
n_file_per_job      = 10
cmssw_dir           = "/ufrc/avery/kinho.lo/Delphes/CMSSW_10_0_5/src/"
out_dir             = "/cmsuf/data/store/user/t2/users/klo/Delphes/ALP_HToZaTo2l2g_M1/2020-06-02/"
delphes_dir         = os.getcwd()

# ______________________________________________________________________ ||
input_file_list = [f for f in glob.glob(input_file_pattern) if "inLHE" not in f]
input_file_list.sort(key=lambda x: int(os.path.basename(x).split("_")[-1].replace(".root","")))
n_file = len(input_file_list)
for ijob in range(0,n_file,n_file_per_job):
    start_index = ijob
    if start_index + n_file_per_job > n_file - 1:
        rnd_index = n_file - 1
    else:
        end_index = start_index + n_file_per_job 
    each_job_name = job_name+"_"+str(start_index)+"-"+str(end_index)
    out_file_name = each_job_name+".root"
    commands = "\n".join([
        "cd "+cmssw_dir,
        "scramv1 runtime -sh",
        "cd "+delphes_dir,
        "./DelphesCMSFWLite cards/delphes_card_CMS.tcl "+os.path.join(out_dir,out_file_name)+" "+" ".join(input_file_list[start_index:end_index]),
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
