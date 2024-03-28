#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
cBioPortal数据库等注释相关数据库网络自动化爬取、下载
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


def get_args():
    parser = argparse.ArgumentParser(description='cBioPortal数据库等注释相关数据库网络自动化爬取、下载。')
    parser.add_argument('--mane', help='MANE的summary文件，用于获取每个蛋白的长度', default='MANE.GRCh38.v1.0.summary.txt')
    parser.add_argument('--cbgenes', help='用于cBioPortal数据下载的基因列表', default='cBioportal_download_genelist-20230619.xlsx')
    parser.add_argument('--download_dir', help='cBioPortal下载路径', default=os.getcwd())
    return parser.parse_args()


def get_protein_length(protein):
    """通过MANE RefSeq蛋白id获取蛋白长度"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "id": protein,
        "rettype": 'gb',
        # "retmode": "xml"       # 提供的json不标准，没法用。加retmode=xml反馈xml格式，否则是类似网页默认的genebank形式的txt
    }
    retry_limit = 100   # 最大重试次数
    retry_count = 0
    protein_length = ''
    while retry_count < retry_limit:
        time.sleep(1)
        try:
            # 等价于浏览器输入https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=xx&rettype=gb
            response = requests.get(base_url, params=params, proxies=None, verify=False)
            if response.status_code == 200:
                protein_length = re.search(r'source\s+1\.\.([0-9]+)', response.text).group(1)   # 正则获取蛋白长度
                print(f'蛋白质{protein}查询成功，长度为{protein_length}')
                break
            else:
                retry_count += 1
                print(f"蛋白质{protein}查询失败{response.status_code}，进行第{retry_count}次重试")
        except BaseException as ex:
            print(ex)
            retry_count += 1
            print(f"蛋白质{protein}连接失败，进行第{retry_count}次重试")

    if retry_count >= retry_limit:
        print(f"蛋白质{protein}查询失败次数过多，放弃查询")
        protein_length = 'NA'
    return protein_length


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
        time.sleep(5)      # 等待oncoprint页面出现加载元素
        # 等待oncoprint页面完全加载完成：
        WebDriverWait(driver, 900).until(
            EC.presence_of_element_located(     # 大的"等待中"突变出现
                (By.XPATH, '//*[@id="mainColumn"]/div/div[1]/div[1]'))
        )

        logger.info(f'{gene}:等待OncoPrint页面加载')
        WebDriverWait(driver, 900).until(
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
        WebDriverWait(driver, 900).until(
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
        WebDriverWait(driver, 1800).until(lambda driver: os.path.exists(downloaded_file))

        # 下载文件判断
        try:
            # with open(downloaded_file) as f:
            #     if f.readline().strip() and f.readline().strip():
            #         logger.info(f'{gene}:本地文件已确认，下载成功！！！')
            df = pd.read_table(downloaded_file)
            if (df['Annotation'].str.split(';', expand=True)[0].str.split(':', expand=True)[1].str.strip() == 'NA').all():
                logger.error(f'{gene}:文件下载不完整！！！')
                raise Exception(f'{gene}:文件下载不完整！！！')
            else:
                logger.info(f'{gene}:本地文件已确认，下载成功！！！')
                # 下载文件重命名
                os.renames(downloaded_file, os.path.join(download_dir, f'cBioPortal_{studies}_{gene}.txt'))
                logger.info(f'{gene}:文件已重命名，退出浏览器')
        except BaseException as ex:
            logger.info(f'{gene}:本地文件下载错误！！！{ex}')
            raise
    except BaseException as ex:
        logger.info(f'{gene}:cBioPortal网站数据加载失败！！！{ex}')
        raise


def init_webdriver(download_dir):
    # 创建独立的WebDriver对象便于多进程并行
    options = Options()
    options.add_argument('--headless')           # 设置为后台运行，不显示前台界面，但貌似速度也没快多少，还不方便查错。脚本完善后可使用
    options.add_argument('--ignore-certificate-errors')
    options.add_argument('--disable-extensions')
    options.add_argument("--disable-popup-blocking")
    options.add_argument("--safebrowsing-disable-download-protection")
    options.add_argument("--no-sandbox")
    options.add_argument("--verbose")
    options.add_argument("--log-path=chrome.log")
    prefs = {"download.default_directory": download_dir}
    options.add_experimental_option("prefs", prefs)
    return webdriver.Chrome(ChromeDriverManager().install(), options=options)


def get_mutation_data(gene, download_dir):
    # 创建下载文件夹
    # download_dir = fr'E:\work\01_development\15_新merge表\功能数据\cBioPortal_download\{gene}'
    download_dir = fr'{download_dir}\{gene}'
    os.makedirs(download_dir, exist_ok=True)
    if os.path.exists(fr'{download_dir}\SUCCESSED'):
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


def get_oncokb_by_gene(gene, driver, logger):
    """
    通过模拟浏览器操作的方式，从OncoKB直接获取基因变异的致病性等注释信息
    """
    pass


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
    if os.path.isfile(args.mane):
        mane_df = pd.read_table(args.mane)
        protein_lens = []
        protein_list = mane_df['RefSeq_prot'].tolist()
        # 并行处理会被服务器拒绝造成response.status_code 429，直接for循环
        for protein in tqdm(protein_list):
            protein_length = get_protein_length(protein)
            protein_lens.append(protein_length)
        mane_df['Protein_len'] = protein_lens
        mane_df.to_csv('MANE.GRCh38.v1.1.summary.txt', sep='\t', index=False)
    if os.path.isfile(args.cbgenes):
        logging.basicConfig(level=logging.INFO, format='[%(asctime)s %(levelname)s] %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.info('输入cBioPortal基因下载列表，开始下载')
        # cbdf = pd.read_excel(args.cbgenes)
        # cbdf = cbdf[cbdf['Downloaded'] != 'YES']
        # genelist = cbdf['Gene'].tolist()
        genelist = []
        with open(args.cbgenes) as f:
            for i in f:
                genelist.append(i.strip())
        print(genelist)
        download_dir = args.download_dir
        # genelist = ['ACOD1', 'ACSL3', 'ACSL6']
        with ProcessPoolExecutor(max_workers=2) as executor:
            # 提交任务给 ProcessPoolExecutor
            for gene in genelist:
                executor.submit(get_mutation_data, gene, download_dir)
                # time.sleep(1)





# ======================= 弃用代码======================
# ---------------多进程调用API接口，因太过频繁弃用-------------
# 创建进程池执行器，最大进程数为12
# with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
#     # 提交任务给执行器
#     protein_list = mane_df['RefSeq_prot'].tolist()
#     results = executor.map(get_protein_length, protein_list)
#
#     # 处理结果
#     for protein, length in zip(protein_list, results):
#         protein_lens.append(length)
#         testw.write(protein+'\t'+length+'\n')

# 解析XML数据
# xml_data = response.text
# start_index = xml_data.find("<GBSeq_length>") + len("<GBSeq_length>")
# end_index = xml_data.find("</GBSeq_length>")
# protein_length = int(xml_data[start_index:end_index])

# """
# 使用webdriver模块通过基因名和转录本号获得蛋白号和蛋白长度
# 因MANE自带蛋白号且API更快，弃用此方法。
# """
# # 使用Chrome浏览器
# driver = webdriver.Chrome(ChromeDriverManager().install())
# # 打开ncbi网站
# driver.get("https://www.ncbi.nlm.nih.gov/nuccore")
# # 搜索EGFR基因
# # 方法：找到搜索框，F12或右键检查，<div class="jig-ncbiclearbutton-wrap ui-ncbiclearbutton-wrap"里面的
# # <input type="text" name="term" id="term" disableU title="Use up and down arrows to choose an item from the autocomplete."
# transcript_id = 'NM_005228.5'
# search_box = driver.find_element_by_id("term")
# search_box.send_keys(transcript_id)
# driver.implicitly_wait(3)  # 等待3秒钟，否则可能卡住
# search_box.send_keys(Keys.RETURN)
# # search_button = driver.find_element_by_id("search")   # 点击搜索按钮，时灵时不灵的
# # search_button.click()
#
# # 等待搜索结果页面加载完成
# WebDriverWait(driver, 20).until(
#     EC.presence_of_element_located((By.CLASS_NAME, "genbank"))  # genbank是大class，包含页面下大部分信息
# )
#
# cds_element = driver.find_element_by_id(f'feature_{transcript_id}_CDS_0')  # 根据id查询CDS部分的全部信息
# protein_id = re.search('/protein_id="(NP_.*)"', cds_element.text).group(1)  # 正则获取蛋白id
#
# protein_link = driver.find_element_by_link_text(protein_id)  # 根据查询到的蛋白id获取该超链接
# protein_link.click()  # 点击该蛋白id的超链接
#
# # 等待蛋白结果页面加载完成
# WebDriverWait(driver, 20).until(
#     EC.presence_of_element_located((By.CLASS_NAME, "genbank"))  # genbank是大class，包含页面下大部分信息
# )
#
# locus_element = driver.find_element_by_id(f'feature_{protein_id}_source_0')  # 根据id查询FEATURES-source部分的全部信息
# protein_length = re.search(r'source\s+1\.\.([0-9]+)', locus_element.text).group(1)  # 正则获取蛋白长度


# -----------cbioportal获取变异的API实现------------
# 因仅能获取基本变异信息，无任何注释，故弃用。改用模拟浏览器操作的方式
# 设置API请求参数
# gene_id = 1956
# study_id = "pan_origimed_2020"
# url = f"https://www.cbioportal.org/api/molecular-profiles/{study_id}_mutations/mutations?entrezGeneId={gene_id}&sampleListId={study_id}_all"
#
# # 发送API请求并下载数据
# response = requests.get(url)
# if response.status_code == 200:
#     # 将数据保存到文件中
#     print(response.content.decode("utf-8"))
#     df = pd.read_json(response.content.decode("utf-8"))
#     df.to_csv('other.txt', sep='\t', index=False)
#     print(df)
# else:
#     print(f"Error: {response.status_code}")

# -----------cbioportal获取变异的selenium实现------------
# 创建Chrome浏览器实例
# driver = webdriver.Chrome(r'C:\Program Files\Google\Chrome\Application\chrome.exe')   # 会闪退，不如用新初始化的
# 直接开到mutations页面的方法行不通，换成curated_studies就加载不了了，只能从头开始
# base_url = 'https://www.cbioportal.org/results/mutations'
# tcga_sudies = 'laml_tcga_pan_can_atlas_2018%2Cacc_tcga_pan_can_atlas_2018%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cchol_tcga_pan_can_atlas_2018%2Ccoadread_tcga_pan_can_atlas_2018%2Cdlbc_tcga_pan_can_atlas_2018%2Cesca_tcga_pan_can_atlas_2018%2Cgbm_tcga_pan_can_atlas_2018%2Chnsc_tcga_pan_can_atlas_2018%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Clihc_tcga_pan_can_atlas_2018%2Cluad_tcga_pan_can_atlas_2018%2Clusc_tcga_pan_can_atlas_2018%2Cmeso_tcga_pan_can_atlas_2018%2Cov_tcga_pan_can_atlas_2018%2Cpaad_tcga_pan_can_atlas_2018%2Cpcpg_tcga_pan_can_atlas_2018%2Cprad_tcga_pan_can_atlas_2018%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_tcga_pan_can_atlas_2018%2Cstad_tcga_pan_can_atlas_2018%2Ctgct_tcga_pan_can_atlas_2018%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Cucs_tcga_pan_can_atlas_2018%2Cucec_tcga_pan_can_atlas_2018%2Cuvm_tcga_pan_can_atlas_2018'
# curated_studies = 'paac_jhu_2014%2Cmel_tsam_liang_2017%2Call_stjude_2016%2Caml_ohsu_2018%2Claml_tcga_pan_can_atlas_2018%2Cmnm_washu_2016%2Cacyc_fmi_2014%2Cacyc_jhu_2016%2Cacyc_mda_2015%2Cacyc_mskcc_2013%2Cacyc_sanger_2013%2Cacc_2019%2Cacbc_mskcc_2015%2Cacc_tcga_pan_can_atlas_2018%2Campca_bcm_2016%2Cbcc_unige_2016%2Cblca_mskcc_solit_2014%2Cblca_mskcc_solit_2012%2Cblca_bgi%2Cblca_dfarber_mskcc_2014%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_hta9_htan_2022%2Cbrca_metabric%2Cbrca_mskcc_2019%2Cbrca_smc_2018%2Cbfn_duke_nus_2015%2Cbrca_bccrc%2Cbrca_broad%2Cbrca_sanger%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cpan_origimed_2020%2Cchol_icgc_2017%2Cchol_nccs_2013%2Cchol_nus_2012%2Cchol_tcga_pan_can_atlas_2018%2Clcll_broad_2013%2Ccll_broad_2015%2Ccll_iuopa_2015%2Ccllsll_icgc_2011%2Cccrcc_dfci_2019%2Ccoad_caseccc_2015%2Ccoad_cptac_2019%2Ccoadread_dfci_2016%2Ccoadread_genentech%2Ccoadread_tcga_pan_can_atlas_2018%2Ccoadread_mskcc%2Chccihch_pku_2019%2Ccscc_dfarber_2015%2Ccscc_hgsc_bcm_2014%2Ccscc_ucsf_2021%2Cctcl_columbia_2015%2Cpact_jhu_2011%2Cdesm_broad_2015%2Cdifg_glass_2019%2Cdlbcl_dfci_2018%2Cdlbcl_duke_2017%2Cdlbc_tcga_pan_can_atlas_2018%2Cnhl_bcgsc_2013%2Ccrc_nigerian_2020%2Cucec_cptac_2020%2Cucec_ccr_msk_2022%2Cucec_ccr_cfdna_msk_2022%2Cesca_broad%2Cesca_tcga_pan_can_atlas_2018%2Cescc_icgc%2Cescc_ucla_2014%2Ces_iocurie_2014%2Cgbc_shanghai_2014%2Cegc_tmucih_2015%2Cstad_oncosg_2018%2Cgbm_cptac_2021%2Cgbm_columbia_2019%2Cgbm_tcga_pan_can_atlas_2018%2Cglioma_msk_2018%2Chnsc_broad%2Chnsc_jhu%2Chnsc_tcga_pan_can_atlas_2018%2Cliad_inserm_fr_2014%2Chcc_meric_2021%2Chcc_inserm_fr_2015%2Chistiocytosis_cobi_msk_2019%2Cpanet_shanghai_2013%2Cchol_jhu_2013%2Cihch_ismms_2015%2Cihch_smmu_2014%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_bgi%2Cccrcc_irc_2014%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Chcc_msk_venturaa_2018%2Clihc_amc_prv%2Clihc_riken%2Clihc_tcga_pan_can_atlas_2018%2Clgg_ucsf_2014%2Clgsoc_mapk_msk_2022%2Cluad_broad%2Cluad_cptac_2020%2Cluad_oncosg_2020%2Cluad_tcga_pan_can_atlas_2018%2Cluad_tsp%2Clung_smc_2016%2Clung_nci_2022%2Clusc_cptac_2021%2Clusc_tcga_pan_can_atlas_2018%2Cmsk_impact_2017%2Cmixed_allen_2018%2Cmpnst_mskcc%2Cmcl_idibips_2013%2Cmbl_broad_2012%2Cmbl_dkfz_2017%2Cmbl_pcgp%2Cmbl_sickkids_2016%2Cskcm_mskcc_2014%2Cmng_utoronto_2021%2Cmeso_tcga_pan_can_atlas_2018%2Cbiliary_tract_summit_2022%2Cbrca_igr_2015%2Cmel_dfci_2019%2Cskcm_dfci_2015%2Cskcm_vanderbilt_mskcc_2015%2Cmel_ucla_2016%2Cprad_mich%2Cprad_su2c_2019%2Cmetastatic_solid_tumors_mich_2017%2Cmixed_selpercatinib_2020%2Cmm_broad%2Cmds_tokyo_2011%2Cmds_iwg_2022%2Cmpn_cimr_2013%2Cnpc_nusingapore%2Cnbl_amc_2012%2Cnbl_ucologne_2015%2Cnepc_wcm_2016%2Cnhl_bcgsc_2011%2Cnsclc_mskcc_2018%2Cnsclc_tracerx_2017%2Cnsclc_unito_2016%2Chnsc_mdanderson_2013%2Cov_tcga_pan_can_atlas_2018%2Cpog570_bcgsc_2020%2Cpancan_pcawg_2020%2Cpaad_qcmg_uq_2016%2Cpaad_tcga_pan_can_atlas_2018%2Cpaad_utsw_2015%2Cpaad_cptac_2021%2Cpanet_jhu_2011%2Cpanet_arcnet_2017%2Call_phase2_target_2018_pub%2Caml_target_2018_pub%2Cbrain_cptac_2020%2Ces_dfarber_broad_2014%2Cnbl_target_2018_pub%2Cpediatric_dkfz_2017%2Cmixed_pipseq_2017%2Cpptc_2019%2Crt_target_2018_pub%2Cwt_target_2018_pub%2Cpcpg_tcga_pan_can_atlas_2018%2Cplmeso_nyu_2015%2Cpcnsl_mayo_2015%2Cprad_broad%2Cprad_fhcrc%2Cprad_mskcc%2Cprad_eururol_2017%2Cprad_tcga_pan_can_atlas_2018%2Cprad_mskcc_cheny1_organoids_2014%2Cprostate_dkfz_2018%2Cprad_msk_2019%2Cprostate_pcbm_swiss_2019%2Cbrca_cptac_2020%2Cccrcc_utokyo_2013%2Cnccrcc_genentech_2014%2Cmrt_bcgsc_2016%2Crms_nih_2014%2Csummit_2018%2Csarc_mskcc%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_broad%2Cskcm_tcga_pan_can_atlas_2018%2Cskcm_yale%2Cskcm_broad_brafresist_2012%2Cscco_mskcc%2Csclc_jhu%2Csclc_ucologne_2015%2Csclc_cancercell_gardner_2017%2Cvsc_cuk_2018%2Cstad_pfizer_uhongkong%2Cstad_tcga_pan_can_atlas_2018%2Cstad_utokyo%2Ctgct_tcga_pan_can_atlas_2018%2Cangs_painter_2020%2Cangs_project_painter_2018%2Cbrca_mbcproject_wagle_2017%2Cmpcproject_broad_2021%2Ctet_nci_2014%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Curcc_mskcc_2016%2Cutuc_mskcc_2015%2Cutuc_cornell_baylor_mdacc_2019%2Cutuc_igbmc_2021%2Cutuc_msk_2019%2Cblca_bcan_hcrn_2022%2Cblca_cornell_2016%2Cucs_jhu_2014%2Cucs_tcga_pan_can_atlas_2018%2Cuccc_nih_2017%2Cucec_tcga_pan_can_atlas_2018%2Cum_qimr_2016%2Cuvm_tcga_pan_can_atlas_2018'
# final_url = f'{base_url}?cancer_study_list={curated_studies}&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations&case_set_id=w_mut&gene_list={gene}%253ASOMATIC&geneset_list=%20&tab_index=tab_visualize&Action=Submit'

# 因抓取不到网页下载状态，改为只判断本地文件
# 等待下载中元素出现
# WebDriverWait(driver, 10).until(
#     EC.presence_of_element_located((By.XPATH,
#                                     '/html/body/div[4]/div[2]/div/div/div/span[2]')) or
#     EC.presence_of_element_located((By.CSS_SELECTOR,
#                                     'body > div:nth-child(9) > div.undefined.fade.in.modal > div > div > div > span:nth-child(2)')) or
#     EC.presence_of_element_located((By.XPATH,
#                                     '/html/body/div[4]/div[2]/div/div'))    # 大下载框
# )
# logger.info(f'{gene}:网站已开始下载')
# # 等待下载中元素消失
# WebDriverWait(driver, 3600).until(
#     EC.invisibility_of_element_located((By.XPATH, '/html/body/div[4]/div[2]/div/div/div/span[2]')) or
#     EC.invisibility_of_element_located((By.CSS_SELECTOR, 'body > div:nth-child(9) > div.undefined.fade.in.modal > div > div > div > span:nth-child(2)')) or
#     EC.invisibility_of_element_located((By.XPATH, '/html/body/div[4]/div[2]/div/div'))
# )
# logger.info(f'{gene}:网站已下载完成')

# 决定保留sv和cnv
# logger.info(f'{gene}:取消sv和cnv')
# time.sleep(5)  # 可能有卡顿，多等待一下
# sv_ckbox = driver.find_element_by_xpath(
#     '//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[3]/div[2]/label[2]/input')
# sv_ckbox.click()
# time.sleep(1)
# cnv_ckbox = driver.find_element_by_xpath(
#     '//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[3]/div[2]/label[3]/input')
# cnv_ckbox.click()
# select_sample = Select(driver.find_element_by_xpath('//*[@id="mainColumn"]/div/div[2]/div[2]/div/div/div/div[4]/div[2]/div/div'))
# select_sample.select_by_visible_text('')      # 无法选中，暂不实现


