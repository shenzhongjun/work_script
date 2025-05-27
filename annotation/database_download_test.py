#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
varsome数据库等注释相关数据库网络自动化爬取、下载。双因子验证麻烦，不弄。
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-5-25 16:33:43"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import requests
import argparse
import os
import re
import time
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import requests
from tqdm import tqdm


if __name__ == "__main__":
    # 设置logging模块
    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.INFO)
    # # 创建控制台logging处理器
    # console_handler = logging.StreamHandler()
    # console_handler.setLevel(logging.INFO)
    # console_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    # logger.addHandler(console_handler)
    #
    # args = get_args()
    # options = Options()
    # # service = Service()
    # # options = webdriver.ChromeOptions()
    # options.add_experimental_option("debuggerAddress", "127.0.0.1:9222")
    # driver = webdriver.Chrome(options=options)
    # print(driver.title)
    # target_url = "https://varsome.com/"
    # driver.get(target_url)
    # # time.sleep(120)
    # print(driver.title)
    #
    # # 等待搜索框显示出来
    # WebDriverWait(driver, 300).until(
    #     EC.presence_of_element_located((By.XPATH, '//*[@id="search"]'))
    # )
    # print('输入变异')
    #
    # time.sleep(2)
    # gene_input = driver.find_element_by_xpath(
    #     '//*[@id="search"]')
    # gene_input.send_keys(f'NM_004333.6:c.1406G>A')
    #
    # time.sleep(2)
    # submit_query = driver.find_element_by_xpath('//*[@id="varsome-search-btn"]')
    # submit_query.click()
    # time.sleep(5)  # 等待oncoprint页面出现加载元素
    #
    # # 等待oncoprint页面完全加载完成：
    # WebDriverWait(driver, 900).until(
    #     EC.presence_of_element_located(  # 点击somatic按钮
    #         (By.XPATH, '//*[@id="start"]/div[16]/div/div[2]/div[1]/div/form/div[1]/div[2]/div[2]'))
    # )
    #
    # somatic_btn = driver.find_element_by_xpath('//*[@id="start"]/div[16]/div/div[2]/div[1]/div/form/div[1]/div[2]/div[2]')
    # somatic_btn.click()
    # time.sleep(5)  # 等待oncoprint页面出现加载元素
    #
    # search_btn = driver.find_element_by_xpath('//*[@id="start"]/div[16]/div/div[2]/div[1]/div/form/div[2]/div[3]/button[2]/div/div')
    # search_btn.click()
    # time.sleep(5)
    #
    # # 等待oncoprint页面完全加载完成：
    # WebDriverWait(driver, 900).until(
    #     EC.presence_of_element_located(  # 点击somatic按钮
    #         (By.XPATH, '//*[@id="variant-info"]/div/div[1]/div'))
    # )
    #
    # apilink = driver.find_element_by_xpath('//*[@id="variant-info"]/div/div[2]/div/div/div[4]/span')
    # apilink.click()
    # time.sleep(5)
    #
    # # logger.info(f'{gene}:输入待查询基因')
    # time.sleep(2)

    # 账号：chinabioinfo@126.com 密码：china@varsome
    # 验证码：k4idfv77
    # fu77wy4r
    # soddzfwd

    sites_path = 'E:/work/01_development/15_新merge表/VarSome数据库/varsome_70sites.txt'
    out_path = 'E:/work/01_development/15_新merge表/VarSome数据库/varsome_70sites_geturl.txt'
    varlist = []

    with open(sites_path) as f:
        for line in f:
            varlist.append(line.strip().split('\t'))
    print(varlist)

    with open(out_path, 'w') as w:
        for var in varlist:
            # url = f"https://api.varsome.com/lookup/NM_000321.3:c.1660G%3EA/1038"
            url = f"https://api.varsome.com/lookup/{var[0]}:{var[1].replace('>', '%3E')}/1038"
            params = {
                "add-all-data": "1",
                "add-ACMG-annotation": "1",
                "add-AMP-annotation": "1",
                "annotation-mode": "somatic"
            }

            # 发送 GET 请求
            response = requests.get(url, params=params)

            # 检查请求是否成功
            if response.status_code == 200:
                data = response.json()  # 解析 JSON 数据
                print(data)  # 打印返回的数据
                w.write(data)
            else:
                print(f"请求失败，状态码: {response.status_code}")
                print(response.text)  # 打印错误信息
            break

