# python调用并回显linux命令
import os
def linux_command(command):
    try:
        with os.popen(command, "r") as p:
            command_return=p.readlines()
    except:
        print(command)
    return command_return


def hostname():
    return linux_command("hostname")[0].replace("\n", "")