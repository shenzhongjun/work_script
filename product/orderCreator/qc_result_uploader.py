#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
在Windows本地将“实验室质检结果反馈.xlsx”表自动上传至集群
因存在堡垒机，未成功实现
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-8-5 09:47:11"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import sys
import time

import paramiko

passwd = 'zhouyj@123'
blip = "172.16.22.251"
blport = 2222
bluser = "zhouyj"
blpasswd = 'zhouyj@123..'
passinfo = "\'s password: "  # ssh 登陆输入密码时的前缀

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(hostname=blip, port=blport, username=bluser, password=blpasswd)
# 直接执行命令无效
# stdin, stdout, stderr = ssh.exec_command('1')
# print(stdout.read())
# stdin2, stdout2, stderr2 = ssh.exec_command('zhouyj@123')
# print(stdout2.read())
# stdin3, stdout3, stderr3 = ssh.exec_command('ls')
# print(stdout3.read())

# 交互式
channel = ssh.invoke_shell()
channel.settimeout(10)
channel.send('1'.encode() + '\n'.encode())
channel.send('\n'.encode())
time.sleep(10)
buff = channel.recv(9999).decode()
print(buff)

# channel2 = ssh.invoke_shell()
# channel2.settimeout(10)
# channel2.send('1'.encode() + '\n'.encode())
# channel2.send('\n'.encode())
# time.sleep(2)
# buff2 = channel2.recv(9999).decode()
# print(buff2)

# while not buff.endswith(passinfo):  # 是否以字符串 's password 结尾
#     channel.send('1'.encode() + '\n'.encode())
#     resp = channel.recv(1024).decode()
#     print(resp)
#     # try:
#     #     resp = channel.recv(1024)
#     #     print(resp.decode())
#     # except Exception as e:
#     #     print(e)
#     #     print('Error info:%s connection time.' % (str(e)))
#     #     channel.close()
#     #     ssh.close()
#     #     sys.exit()
#     buff += str(resp)
#     # if not buff.find('yes/no') == -1:  # 模拟ssh登陆是输入yes
#     #     channel.send('yes\n')
#     # buff = ''

channel.send(f'{passwd}\n'.encode())  # 发送密码
channel.send('ping www.qq.com -c 4\n'.encode())
buff=''
# try:
#     while buff.find('# ')==-1:
#         resp = channel.recv(9999)
#         buff += resp
# except Exception as e:
#     print("error info:" + str(e))

print(buff)
# channel.close()
# ssh.close()
