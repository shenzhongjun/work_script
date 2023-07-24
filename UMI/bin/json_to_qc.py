#读取json文件，获得质控数据

import json
from argparse import ArgumentParser

class QC:
    def __init__(self, json_content, outfile):
        self.json_content = json_content
        self.outfile = outfile

    def deploy(self):
        total_reads = self.json_content["summary"]["before_filtering"]["total_reads"]
        total_bases = self.json_content["summary"]["before_filtering"]["total_bases"]
        raw_data = round(total_bases / 1000000000, 2)
        passed_reads = self.json_content["summary"]["after_filtering"]["total_reads"]
        passed_bases = self.json_content["summary"]["after_filtering"]["total_bases"]
        clean_data = round(passed_bases / 1000000000, 2)
        gc_content = self.json_content["summary"]["after_filtering"]["gc_content"] * 100
        q20_rate = self.json_content["summary"]["after_filtering"]["q20_rate"] * 100
        q30_rate = self.json_content["summary"]["after_filtering"]["q30_rate"] * 100
        dup_rate = self.json_content["duplication"]["rate"] * 100
        peak_insert_size = self.json_content["insert_size"]["peak"]
        adapter_bases = self.json_content["adapter_cutting"]["adapter_trimmed_bases"]
        adapter_rate = round(float(adapter_bases / total_bases) * 100, 2)

        out = open(self.outfile, "w")
        out.write("reads数(PE)\t原始测序数据总碱基数(Gb)\t质控后测序数据总reads数(PE)\t质控后测序数据总碱基数(Gb)\t接头比率(%)\tGC含量(%)\tQ20(%)\tQ30(%)\tDuplication(%)\tPeak_Insert_Size\n")
        result = [total_reads, raw_data, passed_reads, clean_data, adapter_rate, gc_content, q20_rate, q30_rate, dup_rate, peak_insert_size]
        result = [str(i) for i in result]
        out.write("\t".join(result)+"\n")

def main():
    parser = ArgumentParser()
    parser.add_argument("--input", action="store", required=True, help="input file, json format")
    parser.add_argument("--output", action="store", required=True, help="output file")
    argv = parser.parse_args()

    inf = open(argv.input, "r")
    json_content = json.load(inf)
    t = QC(json_content, argv.output)
    t.deploy()

if __name__ == '__main__':
    main()