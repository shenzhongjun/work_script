"""研发任务，VarSome数据库70个肿瘤突变评级SNV位点的数据手动导出成json后整理成table"""
import json

json_path = 'E:/work/01_development/15_新merge表/VarSome数据库/varsome_70sites.json'
table_path = 'E:/work/01_development/15_新merge表/VarSome数据库/varsome_70sites_result.xls'
with open(json_path, encoding='utf-8') as f, open(table_path, 'w', encoding='utf-8') as w:
    w.write('Gene\thgvsp\tclassifi\tCrtd\tCrtd_explain\tDrug\tDrug_explain\tGerm\tGerm_explain\tPath\tPath_explain\t'
            'Pubs\tPubs_explain\tSoma\tSoma_explain\tFreq\tFreq_explain\tType\tType_explain\tPred\tPred_explain\n')
    for line in f:
        if line.strip().startswith('{'):
            json_content = json.loads(line.strip())
        else:
            json_content = json.loads(line.strip())[0]

        gene = json_content['refseq_transcripts'][0]['items'][0]['gene_symbol']
        hgvsp = json_content['refseq_transcripts'][0]['items'][0]['hgvs_p1']
        amp_annotation = json_content['amp_annotation']
        verdict = amp_annotation['verdict']['tier']

        print(f'{gene}:{hgvsp}:{verdict}')
        contents = []
        contents += [gene, hgvsp, verdict]

        for part in amp_annotation['classifications']:
            print(part)
            user_explain = next(iter(part['user_explain'].values()))[0].replace('\n', ' ')
            contents.append(part['tier'])
            contents.append(user_explain)
        w.write('\t'.join(contents) + '\n')
        # print(amp_annotation)
        # break


