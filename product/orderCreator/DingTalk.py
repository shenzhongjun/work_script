#!/mnt/share01/tools/bin/python
# -*- coding: UTF-8 -*-

__author__ = "Fei Xue"
__email__ = "xuefei@zhuanhuayixue.org"
__version__ = "1a"
__date__ = "9/3/21 09:04 "
__copyright__ = "Copyright (C) 2021 All rights reserved"

import time
import hmac
import hashlib
import base64
import urllib.parse
import requests
import json
import getpass
import argparse

USER = getpass.getuser()


class DingMsg:

	def __init__(self, webhook, secret):
		"""
		使用签名密钥进行发送消息
		"""

		self.webhook = webhook
		self.secret = secret

		self.contacts = {"崔思佳": "18522926390",
				 "解飞": "15201052659",
				 "李伟翔": "18811302141",
				 "李秀欣": "15081557830",
				 "刘春月": "13020178149",
				 "任恒达": "17346523190",
				 "闫亚丽": "15201497166",
				 "张现仓": "15632245279",
				 "薛飞": "17696073046",
				 "张佳丽": "17749855194",
				 "祝福": "18704140682",
				 "王子豪": "18861963730",
				 "周义杰": "18810681046"}

	def verify(self):
		timestamp = str(round(time.time() * 1000))
		secret_enc = self.secret.encode('utf-8')
		string_to_sign = '{}\n{}'.format(timestamp, self.secret)
		string_to_sign_enc = string_to_sign.encode('utf-8')
		hmac_code = hmac.new(secret_enc, string_to_sign_enc, digestmod=hashlib.sha256).digest()
		sign = urllib.parse.quote_plus(base64.b64encode(hmac_code))
		webhook = f"{self.webhook}&timestamp={timestamp}&sign={sign}"
		return webhook

	def send_msg(self, msg_type, msg_info, atMobiles, isAtAll=False):
		"""
		发送 文本信息
		@return:
		"""
		webhook = self.verify()
		headers = {"Content-Type": "application/json; charset=utf-8"}
		post_data = {"msgtype": msg_type, "text": {"content": msg_info},
					 "at": {"atMobiles": atMobiles, "isAtAll": isAtAll}
					 }

		r = requests.post(webhook, headers=headers, data=json.dumps(post_data))
		print(r.content)

	def sand_markdown(self, title, markdown_info, atMobiles, atUserIds, isAtAll=False):
		"""
		发送 markdown
		@return:
		"""
		webhook = self.verify()
		headers = {"Content-Type": "application/json; charset=utf-8"}
		post_data = {"msgtype": "markdown",
					 "markdown": {"title": title,  "text": markdown_info},
					 "at": {"atMobiles": atMobiles, "atUserIds": atUserIds,  "isAtAll": isAtAll}}

		r = requests.post(webhook, headers=headers, data=json.dumps(post_data))
		print(r.content)

	def receive_msg(self):
		"""
		接收信息 自动化
		暂时未打通
		@return:
		"""
		pass

	def up_file(self):
		"""
		上传文件 自动化
		暂时未打通
		@return:
		"""
		pass


def my_parser():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('--webhook', help="", required=True)
	parser.add_argument('--secret', help="", required=True)
	parser.add_argument('--head_info', help="", default='')
	parser.add_argument('--msg_type', help="", default='text')
	parser.add_argument('--msg_file', help="", default=None)
	parser.add_argument('--msg_info', help="", default=None)
	parser.add_argument('--atMobiles', help="需要@同事的电话号码，多个的话逗号分割", default=[])
	parser.add_argument('--isAtAll', help="是否全部@", default=False)
	return parser.parse_args()


if __name__ == '__main__':
	args = my_parser()
	webhook = args.webhook
	secret = args.secret
	msg_type = args.msg_type
	if args.msg_file:
		msg_info = args.head_info + '\n' + open(args.msg_file).read()
	else:
		msg_info = args.msg_info

	atMobiles = str(args.atMobiles).strip().split(',')
	isAtAll = args.isAtAll

	ding_msg = DingMsg(webhook, secret)
	ding_msg.send_msg(msg_type, msg_info, atMobiles, isAtAll)
