#对报告中突变频率低于5%的位点，以及回溯的位点，生成html
import os
from argparse import ArgumentParser
import pandas as pd
import yaml
import subprocess

def vcf(somatic, outfile):
	somatic_infor = pd.read_table(somatic)
	somatic_infor = somatic_infor[somatic_infor['is_delete'].isin([0])]
	somatic_infor['AltDepth'] = ''
	for index, row in somatic_infor.iterrows():
		somatic_infor.loc[index, 'AltDepth'] = int(int(row['Depth']) * float(row['MutFreq']) / 100)
	somatic_infor.index = [i for i in range(len(somatic_infor))]
	out = open(outfile, 'w')
	for i in range(len(somatic_infor['Gene'])):
		if somatic_infor['AltDepth'][i] > 50:
			continue
		if (somatic_infor['MutFreq'][i] > 5 and somatic_infor['Ref'][i] != '-' and somatic_infor['Alt'][i] != "-") or (somatic_infor['MutFreq'][i] > 20 and (somatic_infor['Ref'][i] == '-' or somatic_infor['Alt'][i] == "-")):
			continue
		data = [somatic_infor['Chr'][i], somatic_infor['Start'][i], somatic_infor['End'][i], somatic_infor['Ref'][i], somatic_infor['Alt'][i], somatic_infor['HGVS'][i]]
		data = [str(i) for i in data]
		out.write('\t'.join(data)+'\n')
	out.close()

def main():
	parser = ArgumentParser()
	configFile = os.path.join(os.path.dirname(__file__), 'somaticVision.yaml')
	parser.add_argument("--config", action="store", dest="config", default=configFile, help="config file")
	parser.add_argument("--bam", action="store", dest="bam", required=True, help="final bam")
	parser.add_argument("--somatic", action="store", dest="somatic", required=True, help="somatic target mutation")
	parser.add_argument("--outfile", action="store", dest="outfile", required=True, help="output file, a vcf list")
	o = parser.parse_args()
	f = open(o.config, 'r')
	configDic = yaml.load(f, Loader=yaml.FullLoader)
	retrieveList = configDic["database"]["retrievelist"]
	##########################Run##################
	vcf(o.somatic, o.outfile)
	subprocess.getstatusoutput('cut -f 1-5,10 '+retrieveList+' >> '+o.outfile)

if __name__ == '__main__':
	main()
