#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/2/25 16:38
# @Author  : kailiu@acorndx.com
# @Descript:
# @File    : visHtml.py
# @Software: PyCharm

import re

# 颜色设置
misMatchbg = '#00ff00'  # '#97FFFF'#'#006400'   #Lime
misMatchColor = '#0000AA'  # '#FF00FF'#'#FFFF00'   # 黄色
lowqcColor = '#C0C0C0'  # grey
indelColor = '#FF8040'  # Orange
InsertBG = '#EEEE00'
softColor = '#C4C4C4'  # '#6960EC'  #Slate Blue2
softBG = '#8B6914'  # 土色
forwardColor = '#0165FC'  # bright blue
reverseColor = '#EF4026'  # tomato
unequalColor = '#C0D9D9'


htmlhead = '''
<html>
<head>
<meta charset = 'utf-8'>
<title> 变异展示</title>
</head>

<script src="https://cdn.staticfile.org/jquery/1.10.2/jquery.min.js">
</script>
<script>
$(document).ready(function(){{
  $("#hide").click(function(){{
    $("insert").hide();
  }});
  $("#show").click(function(){{
    $("insert").show();
  }});
  $("#hideSoft").click(function(){{
    $("soft").hide();
  }});
  $("#showSoft").click(function(){{
    $("soft").show();
  }});
  $("#showUnequal").click(function(){{
    $("unequal").show();
  }});
  $("#hideUnequal").click(function(){{
    $("unequal").hide();
  }});
}});
</script>
<button id="hide">隐藏插入序列</button> <button id="show">显示插入序列</button>
<button id="hideSoft">隐藏soft序列</button> <button id="showSoft">显示soft序列</button>
<button id="hideUnequal">隐藏变异不绝对相等序列</button> <button id="showUnequal">显示变异不绝对相等序列</button>
<p>
  碱基质量阈值：<font color='red'>{bpQ}</font>
  <br>
  比对质量阈值：<font color='red'>{mapQ}</font>
  <br>
  对于组织样本，变异可信的标准是：至少有4个高质量fragment支持（一个fragment包括read1和read2，两个reads拥有相同的readsID--每条reads倒数第三列信息）
  <br>
  高质量reads：reads上碱基质量大部分大于阈值，尤其是变异位置碱基质量（小于阈值会标记为灰色），同时reads上零散的变异较少（比较干净）
</p>
<pre>
'''

htmltail = '''
</pre>
</html>
'''


def misMatchbase(base, baseQual, cutoff, is_equal):
    color = lowqcColor if (baseQual - 33 < cutoff) else misMatchColor
    if is_equal:
        content = '<span style="color:%s;background-color:%s;">%s</span>' % (color, misMatchbg, base)
    else:
        content = '<unequal><span style="color:%s;background-color:%s;">%s</span></unequal>' % (color, unequalColor, base)
    return content


def indelBase(base, baseQual, cutoff, status='D', is_equal=True):
    if status[0] == 'I':
        color = lowqcColor if (baseQual - 33 < cutoff) else indelColor
        if is_equal:
            bg_color = InsertBG
            type_ = 'insert'
        else:
            bg_color = unequalColor
            type_ = 'unequal'
        content = '<%s><span style="color:%s;background-color:%s;">%s</span></%s>' % (type_, color, bg_color, base, type_)
    else:
        if is_equal:
            content = '<font color = %s>%s</font>' % (indelColor, base)
        else:
            content = '<unequal><span style="color:%s;background-color:%s;">%s</span></unequal>' % (indelColor, unequalColor, base)
    return content


def softBase(base, baseQual, cutoff, is_equal):
    color = lowqcColor if (baseQual - 33 < cutoff) else softColor
    if is_equal:
        bg_color = softBG
        type_ = 'soft'
    else:
        bg_color = unequalColor
        type_ = 'unequal'

    content = '<%s><span style="color:%s;background-color:%s;">%s</span></%s>' % (type_, color, bg_color, base, type_)

    return content


# 相邻的碱基html标签合并
def mergeAdjacentTags(baseList, base):
    if len(baseList) > 0:
        preBase = re.sub(r'(<.+?>)(.).+?(</.+?>)', r'\1\2\3', baseList[-1])
    else:
        preBase = ''
    if base == preBase:
        temp = baseList[-1] + base
        baseList[-1] = re.sub(r'</.+?><.+?>', r'', temp)
    else:
        baseList.append(base)
    return baseList


# 获得序列的基本元素
def baseContent(reads, baseQual, cutoff, is_equal):
    IsReads1 = reads.is_read1
    IsReverse = reads.is_reverse
    if IsReads1:
        char = '-'
    else:
        char = '='
    if IsReverse:
        color = reverseColor
    else:
        color = forwardColor
    if baseQual - 33 < cutoff:
        color = lowqcColor
    if is_equal:
        base = '<font color=%s>%s</font>' % (color, char)
    else:
        base = '<unequal><span style="color:%s;background-color:%s;">%s</span></unequal>' % (color, unequalColor, char)
    return base


# 将cigar转成一个数列展示
def findInsertBase(reads, posInReads):
    cigars = re.findall(r'\d+\w', reads.cigarstring)
    align = []
    for cigar in cigars:
        status = cigar[-1]
        length = int(cigar[:-1])
        align += [status] * length
    return align[posInReads]


def posReference(start, end):
    posRef = [i for i in range(int(start), int(end))]
    return posRef


# 展示reads
def htmlHeader(html,bpQ,mapQ):
    htmlW = open(html, 'w')
    htmlW.write(htmlhead.format(bpQ=bpQ,mapQ=mapQ))
    htmlW.close()


def htmlCoordinate(html, posRef, refFasta):
    htmlW = open(html, 'a')
    content = []
    coordinate = []
    htmlW.write('<font color="#0000CD">')

    preLen = 0
    prePos = 0
    for pos in posRef:
        if pos % 20 == 0:
            preLen = len(str(pos))
            prePos = pos
            # print(length)
            coordinate.append(str(pos))
        elif (pos - prePos < preLen):
            # print(pos,prePos)
            pass
        else:
            # print(pos,prePos,'-')
            coordinate.append(' ')

    for pos in posRef:
        if pos % 20 == 0:
            content.append('|')
        elif pos % 10 == 0:
            content.append(':')
        elif pos % 5 == 0:
            content.append('.')
        else:
            content.append(' ')

    htmlW.write(''.join(coordinate) + '\n')
    htmlW.write(''.join(content) + '\n')
    htmlW.write(refFasta + '\n')
    htmlW.write('</font>')
    htmlW.close()


# 转换变异的展示形式
def mut2string(mut):
    content = '{chrom}:{start}-{end}:{ref}>{allele}'.format(
        chrom=mut[0], start=mut[1], end=mut[2], ref=mut[3], allele=mut[4])
    if len(mut) > 5:
        content += '\t' + mut[5]
    # print(content,mut)
    content = '<font color = "blue">%s</font>' % content
    return content


def htmlMutSite(html, posRef, mut):
    htmlW = open(html, 'a')
    htmlW.write('<font color="#FF0000">')
    content = []
    mutSite = []
    start = mut[1]
    for pos in posRef:
        if pos == int(start):
            content.append('|')
        else:
            content.append(' ')

    length = len(start)
    for pos in posRef:
        if pos == int(start):
            mutSite.append(start)
        elif pos - int(start) < length and pos - int(start) > 0:
            pass
        else:
            mutSite.append(' ')

    htmlW.write(''.join(content) + '\n')
    htmlW.write(''.join(mutSite) + '\n')
    htmlW.write('</font>')
    htmlW.write(mut2string(mut))
    htmlW.close()


def countReads(reads, count):
    count['total'] += 1
    if reads.is_reverse:
        count['reverse'] += 1
        if reads.is_read1:
            count['Rreads1'] += 1
        else:
            count['Rreads2'] += 1
    else:
        count['forward'] += 1
        if reads.is_read1:
            count['Freads1'] += 1
        else:
            count['Freads2'] += 1
    return count


def readsTOtable(mutReads, refReads):
    COUNT = {'type': 'mut',
             'forward': 0,
             'reverse': 0,
             'Freads1': 0,
             'Freads2': 0,
             'Rreads1': 0,
             'Rreads2': 0,
             'total': 0}
    refCOUNT = {'type': 'ref',
                'forward': 0,
                'reverse': 0,
                'Freads1': 0,
                'Freads2': 0,
                'Rreads1': 0,
                'Rreads2': 0,
                'total': 0}
    for reads in mutReads:
        COUNT = countReads(reads, COUNT)
    for reads in refReads:
        refCOUNT = countReads(reads, refCOUNT)

    table = ['''
    <table border = "1" cellspacing = "0" cellpadding = "2">
    <caption style = "text-align:left"> reads 统计 </caption>
    <tr>
    ''']
    for item in sorted(COUNT, reverse=True):
        table.append('<th>%s</th>' % item)
    table.append('</tr>')

    for item in sorted(COUNT, reverse=True):
        table.append('<td>%s</td>' % COUNT[item])
    table.append('</tr>')

    for item in sorted(COUNT, reverse=True):
        table.append('<td>%s</td>' % refCOUNT[item])
    table.append('</tr>')

    table.append('</table>\n')
    # 为EQA 临时修改的程序
    # 可以通过这里条件的调整，给每个变异一个最终的筛选状态
    if COUNT['total'] >= 50:
        status = True
    else:
        status = False
    return ''.join(table), status


# 把reads序列转成html的内容
def reads2string(reads, posRef, bpQ, is_equal):
    start = reads.get_blocks()[0][0]  # 获得reads比对的起始终止位置。
    end = reads.get_blocks()[-1][1]

    refFasta = reads.get_reference_sequence().upper()
    readsFa = reads.query_sequence
    readsQual = reads.qual

    alignment = reads.get_aligned_pairs(with_seq=True)
    content = []
    Store = False

    offSetRef = 0  # 调整参考序列的下标
    offSetReads = 0  # 调整reads的下标，寻找比对的Cigar信息，遇到Delete就加+1，其他不动
    for pos in posRef:
        if pos <= start or pos >= end:
            content.append(' ')
        else:
            if not Store:
                pass
            else:
                continue
            Store = True

            for site in alignment:
                if site[0] == None:
                    # print('This is deletion')
                    temp = indelBase('*', 70, bpQ, 'D', is_equal)
                    content = mergeAdjacentTags(content, temp)
                    offSetRef += 1  # 有一个空位，参考序列快进一位。
                    offSetReads += 1
                elif site[1] == None:
                    # print('This is insertion or soft clip')
                    if findInsertBase(reads, site[0] + offSetReads) == 'I':
                        temp = indelBase(readsFa[site[0]], ord(readsQual[site[0]]), bpQ, 'I', is_equal)
                    else:
                        temp = softBase(readsFa[site[0]], ord(readsQual[site[0]]), bpQ, is_equal)

                    content.append(temp)
                    offSetRef -= 1

                else:
                    # print(site[0],len(refFasta),len(readsFa))
                    if readsFa[site[0]] == refFasta[site[0] + offSetRef]:
                        base = baseContent(reads, ord(readsQual[site[0]]), bpQ, is_equal)
                        content = mergeAdjacentTags(content, base)
                    else:
                        temp = misMatchbase(readsFa[site[0]], ord(readsQual[site[0]]), bpQ, is_equal)
                        content.append(temp)
                        # print(readsFa[site[0]] , refFasta[site[0]])
    if is_equal:
        content.append(reads.query_name + '\t')
        content.append(reads.cigarstring + '\t')
        content.append(str(reads.flag))
    else:
        content.append('<unequal><span >%s\t%s\t%s</span></unequal>' % (reads.query_name, reads.cigarstring, str(reads.flag)))

    # print(content)
    return ''.join(content)


# 变异reads展示
def htmlMutReads(html, mutReads, refReads, bpQ, posRef, is_equals):
    htmlW = open(html, 'a')
    htmlW.write("<font color = '#0000CD'>")
    content, status = readsTOtable(mutReads, refReads)

    htmlW.write(content)
    for index, reads in enumerate(mutReads):
        content = reads2string(reads, posRef, bpQ, is_equals[index])
        htmlW.write(content + '\n')
    htmlW.write("</font>\n")
    htmlW.close()
    return status


def htmlTail(html):
    htmlW = open(html, 'a')
    htmlW.write(htmltail)
    htmlW.close()
