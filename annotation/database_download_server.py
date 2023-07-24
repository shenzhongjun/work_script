#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
cBioPortal数据库等注释相关数据库网络自动化爬取、下载。服务器版
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-5-25 16:33:43"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import re
import time
import logging
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def get_args():
    parser = argparse.ArgumentParser(description='cBioPortal数据库等注释相关数据库网络自动化爬取、下载。')
    parser.add_argument('--cbgenes', help='用于cBioPortal数据下载的基因列表', default='cBioportal_download_genelist-20230621.txt')
    parser.add_argument('--download_dir', help='cBioPortal下载路径', default=os.getcwd())
    return parser.parse_args()


def down_cbioportal_by_gene(gene, studies, download_dir, driver, logger):
    """
    通过模拟浏览器操作的方式，从cbioportal获取基因变异的致病性等注释信息
    """
    logger.info(f'开始下载{gene}基因数据，保存至{download_dir}')

    try:
        # 打开 cBioPortal 网站
        logger.info(f'{gene}:打开cBioPortal网站')
        driver.get('https://www.cbioportal.org/')
        # 等待最后一个study显示出来
        WebDriverWait(driver, 300).until(
            EC.presence_of_element_located(
                (By.XPATH,
                 '//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[1]/div/div/div/div[2]/div[2]/ul/ul[34]/li/label/span[1]'))
            )

        logger.info(f'{gene}:点击studies按钮')
        time.sleep(2)
        if studies == 'tcga':
            studies_bt = driver.find_element_by_xpath('//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[1]/div/div/div/div[2]/div[2]/div[2]/div/div/button[1]')
        else:
            studies_bt = driver.find_element_by_xpath('//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[1]/div/div/div/div[2]/div[2]/div[2]/div/div/button[2]')
        studies_bt.click()

        # 等待选择的总样本数显示出来
        WebDriverWait(driver, 300).until(
            EC.presence_of_element_located((By.XPATH,
                                            '//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[1]/div/div/div/div[1]/div[2]/a[2]'))
        )

        logger.info(f'{gene}:点击查询按钮')
        time.sleep(2)
        query_bt = driver.find_element_by_xpath('//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[2]/a[1]')
        query_bt.click()
        # 等待基因输入框显示出来
        WebDriverWait(driver, 300).until(
            EC.presence_of_element_located((By.XPATH,
                                            '//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[5]/div[2]/div[2]/div[1]/textarea'))
        )

        logger.info(f'{gene}:输入待查询基因')
        time.sleep(2)
        gene_input = driver.find_element_by_xpath(
            '//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[5]/div[2]/div[2]/div[1]/textarea')
        gene_input.send_keys(f'{gene}')     # f'{gene}:SOMATIC'
        # 等待基因格式确认结果
        WebDriverWait(driver, 180).until(
            EC.presence_of_element_located((By.XPATH,
                                            '//*[@id="geneBoxValidationStatus"]/div/div/span[2]'))
        )

        logger.info(f'{gene}:点击提交查询按钮')
        time.sleep(2)
        submit_query = driver.find_element_by_xpath('//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[6]/button')
        submit_query.click()

        # 等待oncoprint页面完全加载完成：
        WebDriverWait(driver, 900).until(
            EC.presence_of_element_located(     # 大的"等待中"突变出现
                (By.XPATH, '//*[@id="mainColumn"]/div/div[1]/div[1]'))
        )

        logger.info(f'{gene}:等待OncoPrint页面加载')
        WebDriverWait(driver, 1800).until(
            EC.presence_of_element_located(     # 小的"等待中"图标出现
                (By.XPATH, '//*[@id="mainColumn"]/div/div/div[2]/div/div/div/div[1]/div[1]')) and
            EC.presence_of_element_located(     # 注释中出现
                (By.XPATH, '//*[@id="mainColumn"]/div/div/div[2]/div/div[1]/div/div[1]/div[2]/div/div/div/span[3]')) and
            EC.presence_of_element_located(     # 渲染中出现
                (By.XPATH, '//*[@id="mainColumn"]/div/div/div[2]/div/div/div/div[1]/div[2]/div/div/div/span[4]'))
        )

        logger.info(f'{gene}:OncoPrint页面注释、渲染中')
        WebDriverWait(driver, 3600).until(
            EC.invisibility_of_element_located(  # 小的"等待中"图标消失
                (By.XPATH, '//*[@id="mainColumn"]/div/div/div[2]/div/div/div/div[1]/div[1]')) and
            EC.invisibility_of_element_located(  # 注释中消失
                (By.XPATH, '//*[@id="mainColumn"]/div/div/div[2]/div/div[1]/div/div[1]/div[2]/div/div/div/span[3]')) and
            EC.invisibility_of_element_located(  # 渲染中消失
                (By.XPATH, '//*[@id="mainColumn"]/div/div/div[2]/div/div/div/div[1]/div[2]/div/div/div/span[4]'))
        )

        logger.info(f'{gene}:OncoPrint页面加载完成，切换到Mutation界面')
        time.sleep(10)
        mutation_bt = driver.find_element_by_xpath('//*[@id="mainColumn"]/div/div/div[2]/ul/li[4]/a')
        mutation_bt.click()
        # 等待页面全部加载完成
        WebDriverWait(driver, 1800).until(
            EC.presence_of_element_located(
                (By.XPATH,
                 '//*[@id="mutationsPageTabs"]/div/div/div/div[2]/div/span/div/div[1]/div/span/span/div/button[2]')) and
            EC.presence_of_element_located((By.ID, 'showMoreButton'))
        )

        logger.info(f'{gene}:Mutation加载完成，开始下载')
        time.sleep(5)
        # 先删除可能存在的本地文件
        downloaded_file = os.path.join(download_dir, 'table.tsv')
        if os.path.exists(downloaded_file):
            os.remove(downloaded_file)

        download_bt = driver.find_element_by_xpath(
            '//*[@id="mutationsPageTabs"]/div/div/div/div[2]/div/span/div/div[1]/div/span/span/div/button[2]')
        download_bt.click()

        # 等待本地文件下载完成
        WebDriverWait(driver, 3600).until(lambda driver: os.path.exists(downloaded_file))

        # 下载文件判断
        try:
            with open(downloaded_file) as f:
                if f.readline().strip() and f.readline().strip():
                    logger.info(f'{gene}:本地文件已确认，下载成功！！！')
                else:
                    raise Exception(f'{gene}:本地文件为空！！！')
            # 下载文件重命名
            os.renames(downloaded_file, os.path.join(download_dir, f'cBioPortal_{studies}_{gene}.txt'))
            logger.info(f'{gene}:文件已重命名，退出浏览器')
        except BaseException:
            logger.info(f'{gene}:本地文件下载错误！！！')
            raise
    except BaseException:
        logger.info(f'{gene}:cBioPortal网站数据加载失败！！！')
        raise


def init_webdriver(download_dir):
    # 创建独立的WebDriver对象便于多进程并行
    options = Options()
    options.add_argument('--headless')           # 设置为后台运行，不显示前台界面，但貌似速度也没快多少，还不方便查错。脚本完善后可使用
    options.add_argument('--ignore-certificate-errors')
    options.add_argument('--disable-dev-shm-usage')     # 禁用/dev/shm目录，使用/tmp目录，否则docker运行时内存不足报错
    options.add_argument('--disable-extensions')
    options.add_argument("--disable-popup-blocking")
    options.add_argument("--safebrowsing-disable-download-protection")
    options.add_argument("--no-sandbox")
    options.add_argument("--verbose")
    options.add_argument("--log-path=chrome.log")
    prefs = {"download.default_directory": download_dir}
    options.add_experimental_option("prefs", prefs)
    return webdriver.Chrome(options=options)


def get_mutation_data(gene, download_dir):
    # 创建下载文件夹
    download_dir = fr'{download_dir}/{gene}'
    os.makedirs(download_dir, exist_ok=True)
    if os.path.exists(fr'{download_dir}/SUCCESSED'):
        return

    # 设置logging模块
    logger = logging.getLogger(gene)
    logger.setLevel(logging.INFO)
    # 创建控制台logging处理器
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    logger.addHandler(console_handler)

    retry_limit = 3
    retry_count = 0
    while retry_count < retry_limit:
        driver = init_webdriver(download_dir)

        file_handler = logging.FileHandler(f'{download_dir}/try_curated_count{retry_count + 1}_log.txt')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
        logger.addHandler(file_handler)
        try:
            down_cbioportal_by_gene(gene, 'curated', download_dir, driver, logger)
            logger.removeHandler(file_handler)
            with open(f'{download_dir}/SUCCESSED', 'w'):
                pass
            return
        except BaseException as ex:
            logger.error(f'{gene}:第{retry_count+1}次下载失败！')
            logger.error(ex)
            logger.removeHandler(file_handler)
            retry_count += 1
        finally:            # 即便try成功，return了，finally里的代码仍会被执行
            time.sleep(2)
            driver.quit()

    logger.info(f'{gene}:使用curated studies下载失败，换成tcga studies继续尝试！')
    retry_count = 0
    while retry_count < retry_limit:
        driver = init_webdriver(download_dir)

        file_handler = logging.FileHandler(f'{download_dir}/try_tcga_count{retry_count + 1}_log.txt')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
        logger.addHandler(file_handler)
        try:
            down_cbioportal_by_gene(gene, 'tcga', download_dir, driver, logger)
            logger.removeHandler(file_handler)
            with open(f'{download_dir}/SUCCESSED', 'w'):
                pass
            return
        except BaseException as ex:
            logger.error(f'{gene}:第{retry_count + 1}次下载失败！')
            logger.info(ex)
            logger.removeHandler(file_handler)
            retry_count += 1
        finally:
            time.sleep(2)
            driver.quit()

    logger.info(f'基因{gene}下载失败次数过多，本次放弃下载！')


if __name__ == "__main__":
    # 设置logging模块
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # 创建控制台logging处理器
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    logger.addHandler(console_handler)

    args = get_args()
    if os.path.isfile(args.cbgenes):
        logging.basicConfig(level=logging.INFO, format='[%(asctime)s %(levelname)s] %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.info('输入cBioPortal基因下载列表，开始下载')
        genelist = []
        with open(args.cbgenes) as f:
            for i in f:
                genelist.append(i.strip())
        download_dir = args.download_dir

        with ProcessPoolExecutor(max_workers=16) as executor:
            # 提交任务给 ProcessPoolExecutor
            for gene in genelist:
                executor.submit(get_mutation_data, gene, download_dir)
                # time.sleep(1)     # 每投一个任务sleep10s，全部投上以后才开始执行！不是投一个执行一个！



