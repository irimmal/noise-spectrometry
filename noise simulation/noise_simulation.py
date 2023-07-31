#!/usr/bin/env python3
'''
指定磁场 B，特征时间 T1 和涨落的方差 var，模拟噪声的观测
'''
import argparse
import numpy as np
import scipy
import numpy.random as rd

psr = argparse.ArgumentParser()
psr.add_argument("opt", type=str, help="output file")
psr.add_argument("--var", type=float, help="涨落的方差")
psr.add_argument("--T1", type=float, help="特征时间")
psr.add_argument("-B", type=float, help="磁感应强度")
args = psr.parse_args()

# 在以下的程序中，使用
# args.opt 为输出文件名
# args.var 为涨落的方差
# args.T1 为特征时间
# args.B 为磁感应强度

# 以下为模拟噪声的程序
# 单一磁矩测试，由于无相互作用项，可先对单一磁子进行演化