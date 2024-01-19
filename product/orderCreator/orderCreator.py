#!/mnt/share01/tools/miniconda/bin/python
# -*- coding: utf-8 -*-

"""
依据下机样本自动创建钉钉任务。
功能：
1. √ 每隔半小时自动扫描/mnt/share05/data/product_raw_data/rawdata/script/nas_data/外送测序信息汇总-肿瘤全部信息.txt，并判断是否生成新任务
2. √ 批次数据按订单号和批次排序，使用订单号到lims抓取订单信息，生成项目标题
3. 功能详细
    3.1 卡替：NCWm需要看样本备注是否有组织(跟lims沟通加此词条)
    √ 3.2 多组合：一个订单多个组合如卡替实体瘤+免疫组库，NBZ6 wes+panel+NBDI，把免疫组库拆出来单独建任务（因下机时间差太远）
    √ 3.3 缺少样本无法执行的也创建项目（有的以历史订单为对照），并发送提示、写入文件
    3.4 如果lims有订单备注，一并发送钉钉群(跟lims沟通加此词条)
    √ 3.5 如果生产群有单独备注，抓取并添加
    3.6 如果某批次已生成订单信息后又检测到有内容更新，发送钉钉群
4. √ 将项目标题+样本信息发送到钉钉
5. √ 并标注是否是此订单样本第一次下机（比对总文件确定），如果是第一次下机需要钉钉创建项目，否则只需要更新样本信息即可
6. √ 将本批次所有订单及样本详细信息写入文件，将本批次输出信息写入log文件
7. 其它：√ WPDL1替换为PDL1
         √ 过滤NHC遗传订单
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "0.1.0"
__date__ = "2021-10-29 10:16:38"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import base64
import hashlib
import hmac
import json
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

sys.path.append('/mnt/share02/xuefei/mysql/script/')
import insert_db

def get_args():
    parser = argparse.ArgumentParser(description='根据下机批次自动创建钉钉任务。')
    parser.add_argument('--batch', '-b', help='批次，支持6742-6745及逗号分割格式')
    parser.add_argument('--appendix', '-a', help='实验室和MSL沟通后的订单附加信息表路径',
                        default=f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/实验室质检结果反馈.xlsx')
    parser.add_argument('--log', '-l', help='是否输出log，默认不输出', action='store_true')
    parser.add_argument('--dingd', '-d', help='是否发送钉钉，默认不发送', action='store_true')
    parser.add_argument('--test', '-t', help='表明当前在进行测试', action='store_true')
    return parser.parse_args()


class OrderCreator(object):
    def __init__(self, batch_data, appendix_data, log, dingd, newbatch, test):
        self.batch_data = batch_data
        self.newbatch = newbatch
        self.log = log
        self.dingd = dingd
        self.test = test
        self.path = os.path.split(os.path.realpath(sys.argv[0]))[0]
        self.recordfile = f'{self.path}/下机任务创建.txt'
        # self.recordfile = f'{self.path}/下机任务创建.test.txt'  # 测试用
        self.real_batchs = sorted({line.split(',')[0] for line in self.batch_data})
        # self.orders = {'DDN23026194'}   # 测试用DDN23026601
        self.orders = self.get_order()
        self.df = pd.DataFrame(
            columns=['批次', '订单', '任务名', '样本名', '产品代码', '芯片类型', '订单备注', '创建日期'])
        self.appendix_data = appendix_data

    def get_order(self):
        if self.newbatch:  # 如果是新批次，直接取值
            return sorted({line.split(',')[1] for line in self.batch_data})
        else:  # 如果是卡替原批次插队更新，排除掉已创建的订单
            return sorted({line.split(',')[1] for line in self.batch_data if
                           line.split(',')[1] not in open(self.recordfile).read()})

    def creat_order(self):
        split_line = '-' * 80   # '\n' + '测试专用' + '-' * 30
        # split_line = ''
        localtime = time.strftime("%Y年%m月%d日 %H:%M:%S")  # %Y/%m/%d %H:%M:%S
        logtime = time.strftime("%Y_%m_%d_%H%M%S")
        msg = '###################测试信息，请忽略！！！#######################\n' if self.test else ''
        if self.newbatch:
            msg += f"{localtime}\n下机表中更新以下批次：\n{','.join(sorted(self.real_batchs))}\n{split_line}\n"
        else:
            msg += f"{localtime}\n下机表中批次{self.real_batchs[0]}内容有更新：\n{split_line}\n"

        print(msg)

        for order_num in self.orders:
            request = urllib.request.Request(f"https://lims-api.chigene.cn/api/v1/tumor/orders/{order_num}",
                                             headers={"token": "zhiyindongfang_zhongliu_2020"})
            response = urllib.request.urlopen(request)
            dict_response = json.loads(response.read().decode('utf-8'))["data"]
            sample_data = [line for line in self.batch_data if line.split(',')[1] == order_num]
            appendix = self.get_appendix(order_num)
            order = Order(order_num, dict_response['title'], dict_response['family_member_code'],
                          dict_response['products'], sample_data, appendix)
            self.df = self.df.append(
                pd.DataFrame.from_records(
                    order.data,
                    columns=['批次', '订单', '任务名', '样本名', '产品代码', '芯片类型', '订单备注', '创建日期'],
                ),
                ignore_index=True
            )

            for task in order.tasks:
                note = ''
                if task:
                    if subprocess.getoutput(f'grep \"{task["title"]}\" {self.recordfile}'):
                        note = '注：该订单已创建任务，请勿重复创建\n'
                    print(f"{task['title']}\n{task['samples']}\n{note}{split_line}")
                    msg += f"{task['title']}\n{task['samples']}\n{note}{split_line}\n"
                    try:
                        insert_db.insert_db(time.strftime("%Y-%m-%d"), order.order, self.real_batchs, task['title'], '', task['samples'])
                    except BaseException as e:
                        print(f"{order.order} {task['title']}写入数据库错误！请关注@薛飞\n{e}\n")

        if self.dingd:
            dingd = DingDingWebHook()
            dingd.send_text(msg, '18810681046')

        if len(self.real_batchs) == 1:
            log_name = f'{self.path}/log/b{self.real_batchs[0]}_{logtime}_log.txt'
        else:
            log_name = f'{self.path}/log/b{self.real_batchs[0]}_b{self.real_batchs[-1]}_{logtime}_log.txt'
        if self.log:
            with open(log_name, 'w') as w:
                w.write(msg)
            self.df.to_csv(f'{self.recordfile}', mode='a', header=False, sep='\t', index=False)

    def get_appendix(self, order_num):
        def addinfo(row, dic):
            if '是' in row['不做RNA'].astype('str').values:
                dic['title'].add('无RNA')
            if '是' in row['不做PDL1'].astype('str').values:
                dic['title'].add('无PDL1')
            if '是' in row['加急'].astype('str').values:
                dic['title'].add('加急')
            for i in row['监测用对照订单号']:
                if not pd.isna(i):
                    dic['info'].add(f"监测订单，对照在{i}")
            return dic

        def addinfo2(row, dic):
            if '是' in row['风险建库'].astype('str').values:
                dic['info'].add('风险建库')
            if '是' in row['风险杂交'].astype('str').values:
                dic['info'].add('风险杂交')
            if '是' in row['风险上机'].astype('str').values:
                dic['info'].add('风险上机')
            return dic

        append_dict = {'title': set(), 'info': set()}
        msl_df = self.appendix_data['MSL反馈']
        extract_df = self.appendix_data['提取质检反馈']
        library_df = self.appendix_data['文库质检反馈']
        msl_row = msl_df.loc[msl_df['订单号'] == order_num]
        if not msl_row.empty:
            append_dict = addinfo(msl_row, append_dict)
        extract_row = extract_df.loc[extract_df['订单号'] == order_num]
        if not extract_row.empty:
            append_dict = addinfo(extract_row, append_dict)
        library_row = library_df.loc[library_df['订单编号'] == order_num]
        if not library_row.empty:
            append_dict = addinfo2(library_row, append_dict)
        return {'title': '_'.join(sorted(append_dict['title'])),
                'info': '，'.join(sorted(append_dict['info']))}


class Order(object):
    def __init__(self, order, name, id, products, sample_data, appendix):
        self.order = order
        self.name = name
        self.id = id
        self.products = products
        self.data = []  # 储存全部样本信息，用来输出到记录文件
        self.tasks = []
        self.packages = self.trim_packages()
        self.product_in_package = []
        self.sample_data = sample_data
        self.time = time.strftime("%Y/%m/%d")
        self.appendix = appendix
        self.create_tasks()

    def trim_packages(self):
        """去除单细胞等生产流程没有的产品组合代码"""
        packages = {product['package_code'] for product in self.products}
        trim_list = ['NBESR', 'QPK', 'S0301', 'NHC', 'NRS0301', 'ScRNAseq', 'ScRNAseq-V2', 'ScTCRseq', 'SR01', 'SH01']  # 'Si-D7', 'Si-D14', 'Si-D21', 'Si-D28'
        return [package for package in packages if package not in trim_list]

    def trim_products(self):
        """去除不体现在任务标题里的产品代码"""
        trim_list = ['WARV7', 'WHE', 'WMGMT', 'NCP980m1', 'NCP980m2', 'NCP980m3', 'NBESR', 'QPK', 'S0301', 'S0301J', 'NRS0301', 'NRS0301J',
                     'NCP116m1', 'NCP116m2', 'NCP116m3', 'Ncecm1', 'Ncecm2', 'Ncecm3', 'ScRNAseq', 'ScRNAseq-V2', 'ScTCRseq',
                     'MRD01', 'MRD02', 'MRD03', 'MRD04', 'MRD05', 'MRD06', 'MRD07', 'MRD08', 'I0301']
        new_products = []
        for product in self.product_in_package:
            if product not in trim_list:
                new_products.append(product.replace('WPDL1', 'PDL1'))
        return new_products

    def create_tasks(self):
        """分别生成本订单每个产品组合的任务，包括任务标题和所含样本"""
        # 卡替TCR免疫组库'I0101'、'I0102'合并
        if 'I0101' in self.packages and 'I0102' in self.packages:
            sample_in_package = self.add_sample('kt_tcr', 'kt_tcr')
            self.tasks.append(self.create_info('', sample_in_package, 'kt_tcr'))
            self.packages.remove('I0101')  # 把免疫组库从组合里删除以免重复创建
            self.packages.remove('I0102')
        # 卡替BCR免疫组库'I0201'、'I0202'合并
        if 'I0201' in self.packages and 'I0202' in self.packages:
            sample_in_package = self.add_sample('kt_bcr', 'kt_bcr')
            self.tasks.append(self.create_info('', sample_in_package, 'kt_bcr'))
            self.packages.remove('I0201')  # 把免疫组库从组合里删除以免重复创建
            self.packages.remove('I0202')
        # 血肿组合里的'NBDI'免疫组库拆分出来单独创建任务
        if len(self.contain_NBDI()) > 0:
            self.packages.append('NBDI')
            sample_in_package = self.add_sample('NBDI', 'NBDI')
            self.tasks.append(self.create_info('', sample_in_package, 'NBDI'))
            self.packages.remove('NBDI')  # 把NBDI免疫组库从组合里删除以免重复创建
            self.products.remove(self.contain_NBDI()[0])  # 把NBDI免疫组库从产品里删除，否则免疫组库样本会出现在其它组合里
        for package in self.packages:
            print(package)
            sample_in_package = self.add_sample(package, 'normal')
            self.tasks.append(self.create_info(package, sample_in_package, 'normal'))

    def add_sample(self, package, product_type):
        """添加本批次此产品组合的所有样本，如果是免疫组库或有样本备注则保留备注"""
        sample_in_package = []
        if product_type == 'kt_tcr':
            package = 'kt_tcr'
            products = [{'product_code': 'I0101', 'package_code': 'kt_tcr'},
                        {'product_code': 'I0102', 'package_code': 'kt_tcr'}]
        elif product_type == 'kt_bcr':
            package = 'kt_bcr'
            products = [{'product_code': 'I0201', 'package_code': 'kt_bcr'},
                        {'product_code': 'I0202', 'package_code': 'kt_bcr'}]
        elif product_type == 'NBDI':
            products = [{'product_code': 'NBDI', 'package_code': 'NBDI'}]
        else:
            products = self.products

        for product in products:
            if product_type == 'normal' and product['package_code'] == package:
                self.product_in_package.append(product['product_code'])
            for line in self.sample_data:
                # 只确认产品代码会把其他组合的样本也纳入，因此要确认产品代码+组合代码
                if line.split(',')[4] == product['product_code'] and product['package_code'] == package:
                    sample_in_package.append(line)
        return sample_in_package

    def create_info(self, package, sample_in_package, product_type):
        """生成任务信息"""
        self.product_in_package = self.trim_products()
        # print(self.order, self.packages, self.products)

        if product_type == 'normal' and len(self.product_in_package) > 0 and package == self.product_in_package[-1]:
            title_suffix = package
        elif product_type == 'normal' and len(self.product_in_package) == 0:
            # 防止新出现组合代码时报错并提示
            title_suffix = package
            print(f'此组合中无产品，请及时更新产品组合代码！！！{self.order}---{self.packages}')
        elif product_type == 'kt_tcr':
            title_suffix = '免疫组库(AB+GD)'
        elif product_type == 'kt_bcr':
            title_suffix = '免疫组库(H+KL)'
        elif product_type == 'NBDI':
            title_suffix = 'NBDI'
        else:
            title_suffix = f"{package}_{'+'.join(self.product_in_package)}"

        title = f"{self.order}_{self.name}_{title_suffix}"
        if self.appendix['title']:  # 如果有生产群备注信息如无RNA等，添加到标题
            title = f"{title}_{self.appendix['title']}"

        sample_in_package_rmnotes = []
        for sample_data in sorted(sample_in_package):
            # print(sample_data)
            batch, order, sample, library, code, chip, notes = sample_data.split(',')
            self.data.append([batch, order, title, sample, code, chip, notes, self.time])
            if product_type == 'normal' and sample_data[-1] == '空':  # 如果是普通订单且无备注，去掉备注信息
                sample_data = ','.join([batch, order, sample, library, code, chip])
            sample_in_package_rmnotes.append(sample_data)

        if len(sample_in_package_rmnotes) > 0:
            samples = '\n'.join(sample_in_package_rmnotes)
            if self.appendix['info']:  # 如果有生产群备注信息如监测订单、风险建库等，添加到样本信息前
                samples = f"{self.appendix['info']}\n{samples}"
            task = {'title': title, 'samples': samples}
            return task
        else:
            # samples = '多组合订单，本次未下机此组合的样本'
            return ''


    def contain_NBDI(self):
        return [product for product in self.products if product['product_code'] == 'NBDI']


class DingDingWebHook(object):
    def __init__(self, secret=None, url=None):
        if not secret:
            secret = 'SECcd0589445bb546f60b9e471c918548d5775a737d3d3cb25cffeb2dd51247adc3'  # 原始秘钥
        if not url:
            url = 'https://oapi.dingtalk.com/robot/send?access_token=4b25b5853f696cdfd1a0bd0e968f7d1549d30bd40fd484674e95111a6327aaea'  # 无加密的url
        timestamp = round(time.time() * 1000)  # 时间戳
        string_to_sign = f'{timestamp}\n{secret}'.encode('utf-8')  # 时间戳+秘钥，转为utf8
        hmac_code = hmac.new(secret.encode('utf-8'), string_to_sign, digestmod=hashlib.sha256).digest()
        sign = urllib.parse.quote_plus(base64.b64encode(hmac_code))  # 最终签名
        self.webhook_url = f"{url}&timestamp={timestamp}&sign={sign}"  # 最终url，url+时间戳+签名

    def send_text(self, msg, atMobiles, isAtAll=False):
        """发送文本信息"""
        send_data = {
            "msgtype": 'text',
            "text": {"content": msg},
            "at": {"atMobiles": atMobiles, "isAtAll": isAtAll}
        }
        self.send_data(send_data)

    def send_markdown(self, title, text, atMobiles, isAtAll=False):
        """发送markdown信息。注：长信息会截短、自动换行"""
        send_data = {
            "msgtype": "markdown",
            "markdown": {"title": title, "text": text},
            "at": {"atMobiles": atMobiles, "isAtAll": isAtAll}
        }
        self.send_data(send_data)

    def send_data(self, send_data):
        headers = {"Content-Type": "application/json", "charset": "utf-8"}
        send_data = json.dumps(send_data).encode("utf-8")  # 将字典类型数据转化为json格式，并编码为UTF-8格式
        request = urllib.request.Request(url=self.webhook_url, data=send_data, headers=headers)  # 发送请求
        opener = urllib.request.urlopen(request)  # 将请求发回的数据构建成为文件格式
        print(opener.read())  # 打印返回的结果


def format_batchs(batchs):
    """处理参数输入的批次信息，支持','分割和'-'连续的批次"""
    batch_list = []
    if ',' in batchs and '-' in batchs:
        a = batchs.split(',')
        for i in a:
            if '-' in i:
                start, end = i.split('-')
                for l in range(int(start), int(end) + 1):
                    batch_list.append(str(l))
            else:
                batch_list.append(str(i))
    elif ',' in batchs:
        for i in batchs.split(','):
            batch_list.append(str(i))
    elif '-' in batchs:
        start, end = batchs.split('-')
        for l in range(int(start), int(end) + 1):
            batch_list.append(str(l))
    else:
        batch_list.append(batchs)
    return batch_list


def get_batch_data(batchs):
    ngs_file = '/mnt/share05/data/product_raw_data/rawdata/script/nas_data/外送测序信息汇总-肿瘤全部信息.txt'
    codes = '/mnt/share05/data/product_raw_data/rawdata/script/cancer_code.list'
    for b in batchs:
        subprocess.call(
            f"grep -P '\t{b}\t' {ngs_file} | grep -w -f {codes} | grep -E -v 'ZRNDZL|ZRnDZL|NHCeF|HLA-|mNGS-seq|NHC|S0301|NBESR|NRS0301J|S0301J|ScRNAseq|ScRNAseq-V2|SR01|ScTCRseq|SH01|S2H01' >> batch_tmp.txt",
            shell=True)
    data = subprocess.getoutput(f"awk -F '\t' '{{OFS=\",\"}}{{$1=$1;gsub(/,/, \"，\", $12); print $3,$15,$4,$5,$2,$7,$12}}' batch_tmp.txt | "
                                f"sort -k 2,2 -k 1n,1 -t,")
    os.remove('batch_tmp.txt')
    return data


if __name__ == "__main__":
    args = get_args()
    appendix_data = pd.read_excel(args.appendix, header=1, sheet_name=None)
    if args.batch:  # 使用参数手动运行
        batchs = format_batchs(args.batch)
        batch_data = get_batch_data(batchs)
        if batch_data:
            c = OrderCreator(batch_data.split('\n'), appendix_data, args.log, args.dingd, True, args.test)
            c.creat_order()
        else:
            print('该批次数据未下机！')
    else:  # 无参数，使用nohup投递后台每个半小时自动运行
        ngs_file = '/mnt/share05/data/product_raw_data/rawdata/script/nas_data/外送测序信息汇总-肿瘤全部信息.txt'
        # ngs_file = f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/外送测序信息汇总-肿瘤全部信息.test.txt'    # 测试用！！！
        record_path = f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/下机任务创建.txt'
        # record_path = f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/下机任务创建.test.txt'   # 测试用！！！
        file_size = os.stat(ngs_file).st_size
        while True:
            try:
                current_size = os.stat(ngs_file).st_size
                ngs_batch = subprocess.getoutput(f'cut -f3 {ngs_file}|sed "1d"|sort -u|tail -1')
                record_batch = subprocess.getoutput(f'cut -f1 {record_path}|sed "1d"|sort -u|tail -1')
                # print(current_size, file_size, ngs_batch, record_batch, sep='\n')
                if current_size != file_size:
                    if ngs_batch > record_batch:
                        batchs = format_batchs(f'{int(record_batch) + 1}-{int(ngs_batch)}')
                        batch_data = get_batch_data(batchs)
                        if batch_data:  # 新下机批次可能没有肿瘤数据，因此要判断batch_data是否为空
                            # c = OrderCreator(batch_data.split('\n'), appendix_data, True, False, True)   # 测试用
                            c = OrderCreator(batch_data.split('\n'), appendix_data, True, True, True, False)
                            c.creat_order()
                    elif ngs_batch == record_batch:  # 卡替订单原批次插队
                        if current_size < file_size:    # 某批次样本撤销上机
                            pass
                        else:
                            batchs = format_batchs(ngs_batch)
                            batch_data = get_batch_data(batchs)
                            if batch_data:
                                # c = OrderCreator(batch_data.split('\n'), appendix_data, True, False, False)  # 测试用
                                c = OrderCreator(batch_data.split('\n'), appendix_data, True, True, False, False)
                                c.creat_order()
                    file_size = current_size
                time.sleep(60 * 30)
                # time.sleep(5)   # 测试用
            except BaseException as e:
                dingd = DingDingWebHook()
                dingd.send_text(f'报错：{e}.\n自动创建项目程序终止，请及时重启程序！', '18810681046')
                print(f'程序终止：{e}')
                raise
