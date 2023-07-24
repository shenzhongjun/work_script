"""
提取bamdst软件的结果，并计算0.2x和0.5x
"""
import pandas as pd
import re
from argparse import ArgumentParser

class pick_qc:
    """
    从bamdst的结果文件挑选质控结果
    """
    def __init__(self, depth, stats, outfile, gzGCqual):
        self.depth = depth
        self.stats = stats
        self.outfile = outfile
        self.gzGCqual = gzGCqual

    def pick_result(self):
        raw_reads, raw_data, map_ratio, mapdata_ratio = "", "", "", ""
        PCR_dup_ratio = ""
        target_data_ratio, target_read_ratio = "", ""
        aver_depth, aver_depth_rmdup = "", ""
        for line in open(self.stats, "r").readlines():
            if re.search("^#", line):
                continue
            line_data = line.strip().split("\t")
            if "Raw Reads" in line:
                raw_reads = line_data
            elif "Raw Data" in line:
                raw_data = line_data
            elif "Fraction of Mapped Reads" in line:
                map_ratio = line_data
            elif "Fraction of Mapped Data" in line:
                mapdata_ratio = line_data
            elif "Fraction of PCR duplicate reads" in line:
                PCR_dup_ratio = line_data
            elif "Fraction of Target Data in mapped data" in line:
                target_data_ratio = line_data
            elif "Fraction of Target Reads in mapped reads" in line:
                target_read_ratio = line_data
            elif "[Target] Average depth" in line_data:
                aver_depth = line_data
            elif "[Target] Average depth(rmdup)" in line_data:
                aver_depth_rmdup = line_data
            else:
                continue
        return [raw_reads, raw_data, map_ratio, mapdata_ratio, PCR_dup_ratio, target_data_ratio, target_read_ratio, aver_depth, aver_depth_rmdup]

    def count_cover(self, aver_depth):
        data = pd.read_csv(self.depth, sep="\t", compression='gzip')
        cov_2 = float(aver_depth) * 0.2
        cov_5 = float(aver_depth) * 0.5
        line, line_2, line_5 = 0, 0, 0
        for index, row in data.iterrows():
            if row["Rmdup depth"] >= cov_2:
                line_2 += 1
            if row["Rmdup depth"] >= cov_5:
                line_5 += 1
            line += 1
        ratio_2 = round(line_2 / line * 100, 2)
        ratio_5 = round(line_5 / line * 100, 2)
        return ratio_2, ratio_5

    def read_gzGCqual(self):
        data = pd.read_csv(self.gzGCqual, sep="\t")
        q20_rate = data["Q20%"].values.tolist()[0]
        q30_rate = data["Q30%"].values.tolist()[0]
        gc_content = data["GC%"].values.tolist()[0]
        return [["Q20(%)", q20_rate], ["Q30(%)", q30_rate], ["GC(%)", gc_content]]

    def deploy(self):
        result = self.pick_result()
        aver_depth = result[-1][-1]
        ratio_2, ratio_5 = self.count_cover(aver_depth)
        result.append(["0.2×_coverage", ratio_2])
        result.append(["0.5×_coverage", ratio_5])
        gz_stats = self.read_gzGCqual()
        result.extend(gz_stats)
        out = open(self.outfile, "w")
        out.write("质控项目\t基本信息\n")
        for r in result:
            r = [str(i) for i in r]
            out.write("\t".join(r) + "\n")
        out.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("--depth", action="store", required=True, help="depth file")
    parser.add_argument("--stats", action="store", required=True, help="stats file")
    parser.add_argument("--outfile", action="store", required=True, help="out file")
    parser.add_argument("--gzGCqual", action="store", required=True, help="gzGCqual stats file")
    argv = parser.parse_args()
    p = pick_qc(argv.depth, argv.stats, argv.outfile, argv.gzGCqual)
    p.deploy()

if __name__ == '__main__':
    main()