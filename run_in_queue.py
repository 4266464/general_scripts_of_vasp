import scc_lib
import time
import random


def queue_avail(queue_list):
    wait_time = []
    for queue in queue_list:
        if True:
            scc_lib.linux_command(
                f"bsub -q {self.detail['queue']} -n {self.detail['core']} -o %J.log {version_dir[version]}")[0]
            return
        else:
            wait_time.append([queue, jobs_before / node_dict])
    return sorted(wait_time, key=lambda queue: queue[1])[0][0]


queue and node


def run_in_queue(command):
    # 默认处于工作目录下
    queue_list = ["ckduan", "smallopa", "smallib"]

    job = {"soc": {"tc4600v3": "mpijob /opt/vasp/5.4.4/bin/vasp_ncl",
                   "login01": "mpirun /share/softwares/packages/vasp.5.4.4_oneAPI/bin/vasp_ncl"},
           "gamma": {"tc4600v3": "mpijob /opt/vasp/5.4.4/bin/vasp_gam",
                     "login01": "mpirun /share/softwares/packages/vasp.5.4.4_oneAPI/bin/vasp_gam"},
           "std": {"tc4600v3": "mpijob /opt/vasp/5.4.4/bin/vasp_std",
                   "login01": "mpirun /share/softwares/packages/vasp.5.4.4_oneAPI/bin/vasp_std"}}

    login01_head = '#!/bin/bash\n' \
                   '#SBATCH -J test\n' \
                   '#SBATCH -p compute\n' \
                   '#SBATCH -n 48\n' \
                   '#SBATCH -N 1\n' \
                   '#SBATCH --time=24:00:00\n' \
                   '#SBATCH -o %j.log\n' \
                   '#SBATCH --qos=cpu\n'

    hostname = scc_lib.linux_command("hostname")[0].replace("\n", "")
    if hostname == "login01":
        login01_head += f"{job[command][hostname]}\n"
        f"cat {login01_head} > vasp.sh"


    "vasp.sh"
    "sbatch vasp.sh"
    command = f"bsub -q {self.detail['queue']} -n {self.detail['core']} -o %J.log {version_dir[version]}"

    job_info = scc_lib.linux_command()[0]

    print(job_info, end="")

    with open("jobid", "a+") as q:
        q.write(job_info)

    home_dir = {"tc"}
    with open("/home/phys/" + scc_lib.linux_command("whoami")[0].replace("\n", "") + "/.jobs", "a+") as q:
        q.write(job_info)
    # q.write(time.asctime(time.localtime(time.time()))+"\n")
    q.write(os.getcwd() + "\n")

    with open("/share/home/ckduan/.jobs", "a+") as q:
        q.write(job_info)
    q.write(os.getcwd() + "\n")

    os.chdir(cwd)


return None
