#!/home/xuefei/.conda/envs/RawCopy/bin/python3.9
# -*- coding: UTF-8 -*-

__author__ = "Fei Xue"
__email__ = "xuefei@zhuanhuayixue.org"
__version__ = "1a"
__date__ = "9/3/21 08:58 "
__copyright__ = "Copyright (C) 2021 All rights reserved"

import os
import time
import glob
import json
import subprocess
import pandas as pd


class DingDing:

	def __init__(self):
		self.CONFIG = json.load(open('CONFIG.JSON'))
		self.webhook = self.CONFIG['DingTalk']['数据下机群']['webhook']
		self.secret = self.CONFIG['DingTalk']['数据下机群']['secret']

	def send_dingding_by_file_info(self, email_file, head_info=''):
		shell = f""" \
ssh cancer-master '/mnt/share05/data/product_raw_data/rawdata/script/DingTalk.py \
	--webhook "{self.webhook}" \
	--secret "{self.secret}" \
	--msg_file {email_file} \
	--head_info {head_info} '
"""
		run_code, run_info = subprocess.getstatusoutput(shell)
		print(run_code, run_info)

	def send_dingding_by_string_info(self, msg_info):
		run_code, run_info = subprocess.getstatusoutput(f""" \
ssh cancer-master '/mnt/share05/data/product_raw_data/rawdata/script/DingTalk.py \
	--webhook "{self.webhook}" \
	--secret "{self.secret}" \
	--msg_info "{msg_info}" \
	--atMobiles "17696073046" '
""")
		print(run_code, run_info)


class RawCopy(DingDing):

	def __init__(self):
		super().__init__()
		self.localtime_for_log = time.strftime("%Y.%m.%d_%H:%M:%S")
		self.localtime_for_tag = time.strftime("%Y%m%d_%H%M%S")
		self.complete_json = self.CONFIG['complete_json']
		self.faild_json1 = glob.glob('/mnt/writefq-01/WFlog/fail/*.json')
		self.faild_json2 = glob.glob('/mnt/writefq-02/WFlog/fail/*.json')
		self.raw_path_1 = '/mnt/writefq-01/rawdata'
		self.raw_path_2 = '/mnt/writefq-02/rawdata'
		self.cancer_code = open(self.CONFIG['cancer_code']).read().strip().split('\n')
		self.PATH_RAW = self.CONFIG['PATH_RAW']
		self.COMPLETE = self.CONFIG['COMPLETE']
		self.list_new_down_data = []
		# 程序log
		self.run_log = open('/mnt/rhd/cancer_storage/data/product_raw_data/rawdata/script/run.log', 'at')
		# cp log
		self.cp_log = f'/mnt/rhd/cancer_storage/data/product_raw_data/rawdata/log/cp_log/cp_{self.localtime_for_tag}.log'

	def find_need_deal_finsh_json(self):
		set_pass_json_file = set(open(self.complete_json).read().split('\n'))
		finsh_json1 = set(glob.glob('/mnt/writefq-01/WFlog/finsh/*.json'))
		finsh_json2 = set(glob.glob('/mnt/writefq-02/WFlog/finsh/*.json'))
		need_deal_finsh_json1 = finsh_json1 - set_pass_json_file
		need_deal_finsh_json2 = finsh_json2 - set_pass_json_file
		return need_deal_finsh_json1, need_deal_finsh_json2

	def main(self):
		self.run_log.write(f'\n本轮开始扫描：{self.localtime_for_log}\n')
		need_deal_finsh_json1, need_deal_finsh_json2 = self.find_need_deal_finsh_json()
		info_list = [[need_deal_finsh_json1, self.faild_json1, self.raw_path_1, 'writefq-01'],
					 [need_deal_finsh_json2, self.faild_json2, self.raw_path_2, 'writefq-02']]

		for finsh_json, faild_json, raw_path, flag in info_list:
			print('json', len(finsh_json))
			self.deal_each_finsh_json(raw_path, finsh_json)

		# print(len(self.list_new_down_data))
		# t = pd.DataFrame(self.list_new_down_data, columns=['each_finsh_json', 'slide_id', 'each_library', 'file'])
		# t.to_csv('1.txt', sep='\t', index=False)

		if len(self.list_new_down_data) != 0:
			# columns=['each_finsh_json', 'slide_id', 'each_library', 'file']
			pd_need_copy, pd_cancer_pooling = self.find_need_cp_list(self.list_new_down_data)
			self.run_log.write(f'\t需要拷贝的文件个数为：{str(len(pd_need_copy))}..\n')
			if len(pd_need_copy) != 0:
				#
				# 处理拷贝数据
				#
				self.run_log.write(f'\t需要拷贝的文件信息\n')
				info = [f"\t{x['slide_id']}\t{x['each_library']}\t{os.path.basename(x['file'])}" for index, x in pd_need_copy.iterrows()]
				self.run_log.write('\n'.join(info) + '\n')
				# ******************************************************************************************************
				# 抓取待拷贝数据对应的邮件信息
				# ******************************************************************************************************

				list_all_need_cp_library = list(set(list(pd_need_copy['each_library'])))
				pd_need_cp_data_pooling_info = pd_cancer_pooling[pd_cancer_pooling['文库编号'].isin(list_all_need_cp_library)]

				email_info = f'/mnt/rhd/cancer_storage/data/product_raw_data/rawdata/log/email_info/email_{self.localtime_for_tag}.log'
				pd_need_cp_data_pooling_info.to_csv(email_info, sep='\t', index=False)
				self.run_log.write(f'\tEmail信息写出: {email_info}\n')
				#
				# 发送DingTalk ：提示下机数据信息
				# '订单号', '批次', '文库名称', '文库编号', '代码', '实验使用芯片', '实验备注'
				pd_need_cp_data_pooling_info_for_dingtalk = pd_need_cp_data_pooling_info[['批次', '订单号', '文库名称', '文库编号', '代码', '实验使用芯片']]
				DingTalkFile = f'/mnt/rhd/cancer_storage/data/product_raw_data/rawdata/log/email_info/ding_{self.localtime_for_tag}.log'
				pd_need_cp_data_pooling_info_for_dingtalk.to_csv(DingTalkFile, sep=',', index=False)
				DingTalkFile = str(DingTalkFile).replace("/mnt/rhd/cancer_storage/data", "/mnt/share05/data")
				self.send_dingding_by_file_info(DingTalkFile, head_info="预告：以下数据已下机正在拷贝，拷贝完成后会有邮件路径提示，请知悉")

				# ******************************************************************************************************
				# 进行拷贝
				# ******************************************************************************************************
				# columns=['each_finsh_json', 'slide_id', 'library', 'start', 'end', 'statues', 'log']
				pd_cp_log = self.deal_cp_path(pd_need_copy)
				self.send_dingding_by_string_info(f"本轮拷贝邮件汇总路径: {email_info.replace('/mnt/rhd/cancer_storage/data', '/mnt/share05/data')}")

				# ******************************************************************************************************
				# 处理过的json文件进行记录，下次扫描将会跳过
				# ******************************************************************************************************
				# pd_cp_filed = pd_cp_log[pd_cp_log['statues'] == '失败']
				# pd_cp_filed_json = set(pd_cp_filed['each_finsh_json'])
				# pd_cp_pass_json = set(need_deal_finsh_json1 | need_deal_finsh_json2) - pd_cp_filed_json
				# open(self.complete_json, 'at').write('\n'.join(pd_cp_pass_json) + '\n')

		self.run_log.write(f'\t本轮拷贝程序休眠：30min\n')
		self.run_log.write(f'本轮拷贝程序结束：{time.strftime("%Y.%m.%d_%H:%M:%S")}\n')

	def find_exist_file(self):
		list_exist_file = []
		# list_exist_file.extend(glob.glob('/mnt/rhd/cancer_storage/data/product_raw_data/*'))
		# list_exist_file.extend(glob.glob('/mnt/rhd/cancer_storage/data/product_raw_data/*/*'))
		# list_exist_file.extend([os.path.basename(x) for x in
		# 						open('/mnt/rhd/cancer_storage/data/product_raw_data/complete.txt').read().strip().split('\n')])
		list_exist_file.extend([os.path.basename(x) for x in open(self.CONFIG['COMPLETE']).read().strip().split('\n')])
		return list(set([os.path.basename(x) for x in list_exist_file]))

	def find_need_cp_list(self, list_new_down_data):
		# pd_new_down_data : 新下机的全部数据信息
		pd_new_down_data = pd.DataFrame(list_new_down_data, columns=['each_finsh_json', 'slide_id', 'each_library', 'file'])
		pd_new_down_data.to_csv('全部未记录过的json识别到的下机数据.txt', sep='\t', index=False)
		# pd_new_down_data_of_cancer_code : 新下机中相关肿瘤产品代码的数据信息
		# 获取肿瘤代码 以及测试备注 的汇总表
		pd_cancer_pooling = self.pooling_info()
		cancer_librarys = pd_cancer_pooling['文库编号']
		pd_new_down_data_of_cancer_code_temp = pd_new_down_data[pd_new_down_data['each_library'].isin(cancer_librarys)]
		pd_new_down_data_of_cancer_code = pd_new_down_data_of_cancer_code_temp.copy()

		# 新下机中相关肿瘤产品代码未进行拷贝的的数据信息
		list_exist_file = self.find_exist_file()
		pd_new_down_data_of_cancer_code['basename'] = [os.path.basename(x) for x in pd_new_down_data_of_cancer_code['file']]
		pd_new_down_data_of_cancer_code_of_have_not_copy = pd_new_down_data_of_cancer_code[~pd_new_down_data_of_cancer_code['basename'].isin(list_exist_file)]
		pd_new_down_data_of_cancer_code_of_have_copy = pd_new_down_data_of_cancer_code[pd_new_down_data_of_cancer_code['basename'].isin(list_exist_file)]
		ready_cp_detail = f'/mnt/rhd/cancer_storage/data/product_raw_data/rawdata/log/cp_log/ready_cp_detail_{time.strftime("%Y%m%d_%H%M%S")}.log'
		if not pd_new_down_data_of_cancer_code_of_have_not_copy.empty:
			pd_new_down_data_of_cancer_code_of_have_not_copy.to_csv(ready_cp_detail, sep='\t', index=False)
		pd_new_down_data_of_cancer_code_of_have_copy.to_csv('pd_new_down_data_of_cancer_code_of_have_copy.txt', sep='\t', index=False)
		return pd_new_down_data_of_cancer_code_of_have_not_copy, pd_cancer_pooling

	def pooling_info(self):
		"""
		重复	代码	批次	文库名称	文库编号		barcode	实验使用芯片	捕获文库名称	data/G	实验设计	优先级	实验备注		上机时间	文库构建使用名称	订单号	对应的DNA编号	外部文库编号查重复	Unnamed: 17
		1		NTb01F	5032	9DC102	DYDFB-703-1	48		NT01T		AFGD01		24.0	A		9CT548	2020-01-02	9DC102	DD19012699		9CT548	1
		1		NTb01F	5032	9DC103	DYDFB-703-2	49		NT01T		AFGD01		9.0		A		9CT549	2020-01-02	9DC103	DD19012699		9CT549	1
		"""
		excel_up_seq_1 = self.CONFIG['excel_up_seq_1']
		excel_up_seq_2 = self.CONFIG['excel_up_seq_2']
		excel_up_seq = self.CONFIG['excel_up_seq']
		text_up_seq = self.CONFIG['text_up_seq']
		text_all_cancer_up_seq = self.CONFIG['text_all_cancer_up_seq']
		cp_shell = f'ssh cancer-master cp {excel_up_seq_1} {excel_up_seq_2}'
		run_code, run_log = subprocess.getstatusoutput(cp_shell)
		if run_code != 0:
			print(f'汇总表拷贝出错：{run_log}')

		pd_pooling = pd.read_excel(excel_up_seq)
		pd_pooling['实验备注'] = pd_pooling['实验备注'].fillna(value='空')
		pd_pooling['实验备注'] = pd_pooling['实验备注'].str.replace('\n', '').replace(' ', '').replace('\t', '').replace('\r', '')
		pd_pooling.to_csv(text_up_seq, sep='\t', index=False)
		remarks = self.CONFIG['实验备注关键词']
		remarks_pattern = '|'.join(remarks)
		pd_pooling = pd_pooling[pd_pooling['代码'].isin(self.cancer_code) | pd_pooling['实验备注'].str.contains(remarks_pattern)]
		pd_pooling.to_csv(text_all_cancer_up_seq, sep='\t', index=False)
		return pd_pooling

	def deal_each_finsh_json(self, raw_path, finsh_json):
		for each_finsh_json in finsh_json:
			dict_finsh = json.load(open(each_finsh_json))
			slide_id = dict_finsh['slide']
			list_library = [str(x) for x in list(dict_finsh['speciesBarcodes'].keys())]
			for each_library in list_library:
				self.deal_each_library(each_finsh_json, raw_path, slide_id, each_library)

	def deal_each_library(self, each_finsh_json, raw_path, slide_id, each_library):
		fastq1 = f"{raw_path}/{slide_id}/L01/{slide_id}_L01_{each_library}_1.fq.gz"
		fastq2 = f"{raw_path}/{slide_id}/L01/{slide_id}_L01_{each_library}_2.fq.gz"
		stat_1 = f"{raw_path}/{slide_id}/L01/{slide_id}_L01_{each_library}_1.fq.fqStat.txt"
		stat_2 = f"{raw_path}/{slide_id}/L01/{slide_id}_L01_{each_library}_2.fq.fqStat.txt"
		report = f"{raw_path}/{slide_id}/L01/{slide_id}_L01_{each_library}.report.html"
		if os.path.isfile(fastq1) and os.path.isfile(fastq2):
			self.list_new_down_data.append([each_finsh_json, slide_id, each_library, fastq1])
			self.list_new_down_data.append([each_finsh_json, slide_id, each_library, fastq2])
			if self.CONFIG['参数']['是否要stat文件']:
				self.list_new_down_data.append([each_finsh_json, slide_id, each_library, stat_1])
				self.list_new_down_data.append([each_finsh_json, slide_id, each_library, stat_2])
			if self.CONFIG['参数']['是否要report文件']:
				self.list_new_down_data.append([each_finsh_json, slide_id, each_library, report])
		else:
			self.run_log.write(f'\t不存在 {fastq1}' + '\n')
			self.run_log.write(f'\t不存在 {fastq2}' + '\n')

	def deal_cp_path(self, pd_need_copy):
		# columns=['each_finsh_json', 'slide_id', 'each_library', 'file']
		list_cp_log = []
		for index, line in pd_need_copy.iterrows():
			each_finsh_json = str(line['each_finsh_json'])
			slide_id = str(line['slide_id'])
			each_library = str(line['each_library'])
			seq_down_file = str(line['file'])

			if self.CONFIG['参数']['开发者模式']:
				cp_shell = f'touch {self.PATH_RAW}/{os.path.basename(seq_down_file)}'
			else:
				# 休眠2分钟
				time.sleep(2)
				# cp_shell = f'rsync -avP {seq_down_file} {self.PATH_RAW}'
				seq_down_file_of_rmt = seq_down_file.replace('/mnt/writefq-01/rawdata', 'writefqdata@10.100.1.251:/mnt/rhd/rawdata/result/OutputFq/upload')
				seq_down_file_of_rmt = seq_down_file_of_rmt.replace('/mnt/writefq-02/rawdata', 'writefqdata@10.100.1.252:/mnt/rhd/rawdata/result/OutputFq/upload')
				cp_shell = f"rsync -avP -e 'ssh -p 10122' {seq_down_file_of_rmt} {self.PATH_RAW}"
				if seq_down_file.endswith('gz'):
					seq_down_file_size = os.path.getsize(seq_down_file)
					if seq_down_file_size < 100 * 1024 * 1024:
						self.send_dingding_by_string_info(f'警告: {seq_down_file} 数据小于100M, 请检查~')

			start_time = time.strftime("%Y.%m.%d_%H:%M:%S")
			run_code, run_log = subprocess.getstatusoutput(cp_shell)
			end_time = time.strftime("%Y.%m.%d_%H:%M:%S")
			if run_code == 0:
				list_cp_log.append([each_finsh_json, slide_id, each_library, start_time, end_time, '成功', '-'])
				open(self.COMPLETE, 'at').write(os.path.basename(seq_down_file) + '\n')
			else:
				list_cp_log.append([each_finsh_json, slide_id, each_library, start_time, end_time, '失败', run_log])
		pd_cp_log = pd.DataFrame(list_cp_log, columns=['each_finsh_json', 'slide_id', 'library', 'start', 'end', 'statues', 'log'])
		pd_cp_log.to_csv(self.cp_log, sep='\t', index=False)
		return pd_cp_log


if __name__ == '__main__':

	while True:
		time.sleep(10)
		cp_object = RawCopy()
		try:
			cp_object.main()
			cp_object.run_log.close()
			time.sleep(60 * 30)
		except IOError as e:
			cp_object.send_dingding_by_string_info(f'报错：{e}, 拷贝程序退出，请及时重启程序！')




