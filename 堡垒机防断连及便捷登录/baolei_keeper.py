#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
堡垒机防断连脚本。
第一次运行请跟随提示进行设置，后续即可直接运行。
如需修改运行条件，可带参数运行或直接修改脚本目录下的baolei_keeper.config文件。
"""

import argparse
import os
import time
import pyautogui as pg


def get_args():
    parser = argparse.ArgumentParser(description='Keep baolei connecting when you leave.')
    parser.add_argument('--loc', help='Location of your ssh software icon, you can get it when first run or set \
    --first parameter. It should be like "800,1000".')
    parser.add_argument('--tabs', help='Number of tabs to keep.')
    parser.add_argument('--time', help='Total time to run baolei keeper, you can set like 24h or 3d.')
    parser.add_argument('--interval', help='Time interval to run baolei keeper, it cannot be greater than 10h, \
    because the baolei machine will be disconnected if there is no operation for 10 hours.')
    parser.add_argument('--first', help='Set first time to run if ssh software icon location is changed.',
                        action='store_true')
    return parser.parse_args()


def keeper(lo, tab, rtime, itime):
    localtime = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f'''{localtime}
开始运行...
SSH软件坐标({lo})，预计运行{rtime}，保持标签页{tab}个，运行间隔{itime}。
''')
    time.sleep(2)
    pg.moveTo(lo)
    pg.click()
    rtime = float(rtime.split('h')[0]) * 3600 if 'h' in rtime else float(rtime.split('d')[0]) * 3600 * 24
    itime = float(itime.split('h')[0]) * 3600
    for i in range(int(rtime / itime)):
        for a in range(int(tab)):
            pg.write(f'sleep {itime} &', interval=0.08)     # 有bug，不支持判断大小写锁定状态
            pg.press('Enter', interval=0.15)
            pg.hotkey('ctrl', 'tab', interval=0.15)
        print(f'完成第{i+1}次运行；')
        time.sleep(itime)
    print('\n设定时间到，退出运行。')


def first_run():
    print('首次运行，请根据以下提示进行操作：\n')
    input('请将鼠标放置在SSH软件图标上，然后按回车键。PS：建议将软件图标固定在任务栏，后续无需再进行此操作。')
    lo = pg.position()
    print(f'你的SSH软件图标位置是{lo}\n')
    tab = set_tab('默认保持连接12个标签页，如果要重新设定此数目，请输入具体数值：\n', 12)
    rtime = set_time('默认保持连接12个小时，如果要重新设定时间，请输入具体时间，格式参照18h或3d：\n', '12h', 0, 1000)
    itime = set_time('默认3小时进行一次模拟操作，如果要重新设定时间间隔，请输入具体时间，格式参照3h且不能超过10小时：\n',
                     '3h', 0, 10, True)
    print('完成设置，后续运行默认按照此设置进行。如需修改运行条件，可带参数运行或直接修改脚本目录下的baolei_keeper.config文件。\n')
    with open(config_path(), 'w') as w:
        info = '\t'.join(['loc_x', 'loc_y', 'tab_num', 'total_time', 'interval_time']) + '\n'
        info += '\t'.join([str(lo[0]), str(lo[1]), str(tab), rtime, itime])
        w.write(info)
    keeper(lo, tab, rtime, itime)


def set_tab(prompt, default):
    while True:
        try:
            tab = input(prompt)
            tab = int(tab) if tab else default
            if tab > 0:
                return tab
            else:
                print('请输入一个大于0的整数。\n')
                continue
        except ValueError as ve:
            print(f"请输入一个整数。error is : {ve}\n")
            continue


def set_time(prompt, default, min, max, *interval):
    while True:
        t = input(prompt)
        t = t if t else default
        if format_time(t, min, max, interval):
            return format_time(t, min, max, interval)
        else:
            continue


def format_time(t, min, max, *interval):
    try:
        if 'h' in t and min < float(t.split('h')[0]) < max:
            return t
        elif 'd' in t and min < float(t.split('d')[0]) < max and not interval:
            return t
        else:
            print('时间格式错误，请重新输入。')
            return None
    except ValueError as ve:
        print(f"请输入一个整数或小数。erroe is {ve}")
        return None


def config_path():
    path = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'baolei_keeper.config')
    return path


def set_parameter(obj):
    with open(config_path()) as conf:
        conf.readline()
        x, y, tab, rtime, itime = conf.readline().strip().split('\t')
    if obj.loc:
        x = obj.loc.split(',')[0]
        y = obj.loc.split(',')[1]
    if obj.tabs:
        tab = obj.tabs
    if obj.time:
        rtime = format_time(obj.time, 0, 1000)
    if obj.interval:
        itime = format_time(obj.interval, 0, 10)
    return x, y, tab, rtime, itime


if __name__ == "__main__":
    args = get_args()
    if not os.path.exists(config_path()) or args.first:
        first_run()
    else:
        loc_x, loc_y, tabs, run_time, interval_time = set_parameter(args)
        keeper(pg.Point(loc_x, loc_y), tabs, run_time, interval_time)
