#!/home/phys/qif/anaconda3/bin/python

import argparse
import scc_lib
import logging
import os
import re

logging.basicConfig(level=logging.DEBUG, format="%(levelname)s:%(funcName)s:%(message)s")

def job_state(job_list=[]):
    # return {'R': ['7502', '7816'],
    # 'PD': ['7819', '7820']}
    # state = {"tc4600v3": {"PEND", "RUN"},
    #          "login01": {"PD", "R"}, }
    bjobs = {"tc4600v3": "bjobs",
             "login01": "squeue -u ckduan"}
    hostname = scc_lib.hostname()
    logging.debug(f"hostname = <{hostname}>")

    # state = state[hostname]
    state = {}
    log=""
    # JOBID    USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
    # 4340177  qif     RUN   smallopa   tc4600v3    28*node391  */vasp_gam Feb 22 16:40
    for i in scc_lib.linux_command(bjobs[hostname])[1:]:
        job_split = [j for j in i.split(" ") if not j == ""]
        if len(job_split)<3:
            continue
        log += f"<{job_split[0]}> {job_split[2]}\n"
        path = scc_lib.linux_command(f"grep {job_split[0]} ~/.jobs -A1")[1]
        log += path
        if job_split[2][0]=="R":
            path = path.replace("\n", "")
            for j in scc_lib.linux_command(f"ls -ltr {path}"):
                log += j
        if job_split[2] in state:
            state[job_split[2]].append(job_split[0])
        else:
            state[job_split[2]] = [job_split[0]]
    logging.debug(log)
    return state


def incar_para(incar_path, para):
    value = False
    for i in scc_lib.linux_command("grep " + para + " " + incar_path + "/INCAR"):
        key, temp_value = i.replace("\n", "").split("=")
        if key.replace(" ", "") == para:
            value = temp_value
            break
    if value:
        logging.debug(f"<{para}> is found <{value}> in <{incar_path}>.")
    else:
        logging.debug(f"no <{para}> is found in <{incar_path}>.")
    return value


class vasp_job:
    def __init__(self, path):
        self.path = path
        self.state = "undefined"
        jobid = scc_lib.linux_command("grep " + path + " ~/.jobs -B1")
        log = f"<{path}> is initiated. "
        self.jobid = []
        if jobid:
            log += "associated jobid "
            for i in jobid:
                match_result = jobid_template.search(i)
                if match_result:
                    self.jobid.append(match_result.group(1))
                    log += f"<{self.jobid[-1]}> "
            log += "found"
        else:
            log += "no associated jobid found. "
        logging.debug(log)

    def check_log(self):
        log_list = []
        for i in os.listdir(self.path):
            if i[-4:] == ".log":
                log_list.append(i.split(".")[0])
        if not log_list:
            logging.debug(f"no log file found in <{self.path}>. ")
            self.state = "wait"
        else:
            log = "log file " + " ".join(log_list) + f" found in <{self.path}>. "
            log_list = list(set(log_list).intersection(set(self.jobid)))
            log += (" ".join(log_list) if log_list else "none") + " belong to this job. "
            logging.debug(log)
            if not log_list:
                self.state = "wait"
            else:
                log_list.sort()
                log_file = log_list[-1] + ".log"
                logging.debug(f"current log file is <{log_file}>. ")
                if "writing wavefunctions" not in \
                        scc_lib.linux_command("tail " + self.path + "/" + log_file + " -n1")[0]:
                    logging.debug(f"fail since wf not write according to <{log_file}>. ")
                    self.state = "fail"
                elif scc_lib.linux_command(
                        "grep 'reached required accuracy - stopping structural energy minimisation' " + self.path + "/" + log_file):
                    logging.debug(f"good relax according to <{log_file}>. ")
                    self.state = "done"
                else:
                    nsw = incar_para(self.path, "NSW")
                    nsw = 0 if not nsw else nsw
                    nelm = incar_para(self.path, "NELM")
                    nelm = 60 if not nelm else nelm
                    if int(nsw) > 1:
                        bad_step = 0
                        nsw_step = len(scc_lib.linux_command("grep F= " + self.path + "/OSZICAR"))
                        for i in scc_lib.linux_command("grep F= " + self.path + "/OSZICAR -B1 | grep :"):
                            if [j for j in i.split(" ") if j][1] == nelm:
                                bad_step += 1
                        logging.debug(f"unfinished relax with all ion step = <{nsw_step}>, "
                                      f"bad step = <{bad_step}> according to <{log_file}>. ")
                        self.state = "rerun"
                    else:
                        if [i for i in
                            scc_lib.linux_command("grep F= " + self.path + "/OSZICAR -B1 | grep :")[0].split(" ") if
                            i][1] == nelm:
                            logging.debug(f"insufficient electric step according to <{log_file}>. ")
                            self.state = "rerun"
                        else:
                            logging.debug(f"good static according to <{log_file}>. ")
                            self.state = "done"
        return None


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument("path_file", help="t")
    # parser.add_argument("-a", "--action", action="store_true",
    #                     help="")
    # args = parser.parse_args()
    # logging.debug(f"args.path_file = <{args.path_file}>, args.action = <{args.action}>")
    jobid_template = re.compile("(\d{4,})")
    job_state = job_state()

    # job_list = []
    # with open(args.path_file, "r") as f:
    #     for i in f.readlines():
    #         apath = i.replace("\n", "")
    #         if os.path.isdir(apath):
    #             ajob = vasp_job(apath)
    #             job_intersection = list(set(ajob.jobid).intersection(set(job_state.values())))
    #             if len(job_intersection) == 1:
    #                 job_intersection = job_intersection[0]
    #             if job_intersection:
    #                 logging.debug(f"<{job_intersection}> is running in <{ajob.path}>. ")
    #                 ajob.state = "run"
    #             else:
    #                 ajob.check_log()
    #             job_list.append(ajob)
    #         else:
    #             logging.info(f"<{apath}> not exist. ")
    #
    # for i in job_list:
    #     print(i.state, i.path)
