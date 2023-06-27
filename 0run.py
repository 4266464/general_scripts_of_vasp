#!/home/phys/qif/anaconda3/bin/python
#!/share/home/ckduan/anaconda3/bin/python

import os
import argparse


def linux_command(command):
    try:
        with os.popen(command, "r") as p:
            command_return = p.readlines()
    except:
        print(command)
    return command_return


def runvasp_tc4600v3():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", nargs="+",default=["gam"], help="the command will run on queue")
    parser.add_argument("-q", "--queue", default="ckduan", choices=["ckduan", "test", "smallib"],
                        help="submit to which queue, eg ckduan test smallib")
    parser.add_argument("-n", "--ncore", default=0, type=int,
                        help="number of core")
    # parser.add_argument("-v", "--version", default="gam", choices=["gam", "soc", "std"],
    #                     help="the version of code")

    args = parser.parse_args()
    queue_dic = {"ckduan": [24, ""],
                 "smallib": [24, ""],
                 "test": [20, " -m \"k802 k804\""]}
    version_dic = {"soc": "mpijob /opt/vasp/5.4.4/bin/vasp_ncl",
                   "gam": "mpijob /opt/vasp/5.4.4/bin/vasp_gam",
                   "std": "mpijob /opt/vasp/5.4.4/bin/vasp_std"}

    queue = args.queue + queue_dic[args.queue][1]
    ncore = queue_dic[args.queue][0] if args.ncore == 0 else args.ncore

    if args.command[0]=="python":
        command = "~/anaconda3/bin/python "+" ".join([i.replace("=","-") for i in args.command[1:]])
    else:
        command = version_dic[args.command[0]]

    #job_info = linux_command("bsub -q " + queue + " -n " + str(ncore) + " -o %J.log " + version)[0]
    print("bsub -q " + queue + " -n " + str(ncore) + " -o %J.log " + command)
    #print(job_info, end="")
    # with open("jobid", "a+") as q:
    #     q.write(job_info)
    # with open("/home/phys/" + linux_command("whoami")[0].replace("\n", "") + "/.jobs", "a+") as q:
    #     q.write(job_info)
    #     q.write(os.getcwd() + "\n")
    return None


def runvasp_login01():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--ncore", default=48, type=int,
                        help="number of core")
    parser.add_argument("-v", "--version", default="gam", choices=["gam", "soc", "std"],
                        help="the version of VASP")
    args = parser.parse_args()
    sbatch = ["#!/bin/bash",
              "#SBATCH -J vasp",
              "#SBATCH -p compute",
              "#SBATCH -n " + str(args.ncore),
              "#SBATCH -N " + str(args.ncore // 48),
              "#SBATCH --time=36:00:00",
              "#SBATCH -o %j.log",
              "#SBATCH --qos=cpu",
              "mpirun /share/softwares/vasp/from_SCC/vasp_" + args.version.replace("soc", "ncl")]

    if os.path.isfile("vasp.sh"):
        linux_command("rm vasp.sh")
    with open("vasp.sh", "w") as q:
        for i in sbatch:
            q.write(i + "\n")
    job_info = linux_command("sbatch vasp.sh")[0]
    print(job_info, end="")

    with open("jobid", "a+") as q:
        q.write(job_info)
    with open("/share/home/ckduan/.jobs", "a+") as q:
        q.write(job_info)
        q.write(os.getcwd() + "\n")
    return None

if __name__ == '__main__':
    hostname = linux_command("hostname")[0].replace("\n", "")
    if hostname == "tc4600v3":
        linux_command("module purge")
        linux_command("module load intel/2016.3.210")
        runvasp_tc4600v3()
    elif hostname == "login01":
        runvasp_login01()
