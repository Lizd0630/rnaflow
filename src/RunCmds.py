#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-15 02:19:41
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

from multiprocessing import Pool
import subprocess
import time


def run_sys(func, silent=True):
    t1 = time.asctime(time.localtime(time.time()))
    print(f"{t1} perfoming:\n    {func}\n")
    p = subprocess.Popen(args=func,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         shell=True)
    p.wait()
    t2 = time.asctime(time.localtime(time.time()))
    print(f"{t2} finished:\n    {func}\n")
    stderr = p.stderr.read().decode("utf-8")
    if len(stderr) > 0:
        print(f"Errors:\n{stderr}\n")
    else:
        pass
    if silent:
        pass
    else:
        print("Logs:\n")
        stdout = p.stdout.read().decode("utf-8")
        if len(stdout) > 0:
            print(stdout)
        else:
            print("null")


class RunCmds:
    def __init__(self, cmds, silent):
        self.cmds = cmds
        self.silent = silent

    def running(self, n_jobs):
        p = Pool(n_jobs)
        for cmd in self.cmds:
            p.apply_async(run_sys, (cmd, self.silent))
        p.close()
        p.join()


# import threading
# import subprocess
# import time
#
#
# class RunCmds:
#     def __init__(self, cmds):
#         self.cmds = cmds
#
#     def running(self, n_jobs):
#         threads = []
#         for cmd in self.cmds:
#             t = threading.Thread(target=subprocess.Popen, kwargs={"args": cmd, "shell": "True"})
#             threads.append(t)
#         for t in threads:
#             t.start()
#             while True:
#                 time.sleep(10)
#                 # 判断正在运行的线程数量,如果小于5则退出while循环,
#                 # 进入for循环启动新的进程.否则就一直在while循环进入死循环
#                 if(len(threading.enumerate()) < n_jobs):
#                     break
