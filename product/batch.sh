# 输入pici 6715,6716,6717查看这三个批次的样本下机情况（按订单号排序）；
# 输入pici 6715,6716,6717 1查看带备注的样本下机情况（按订单号排序）

array=(`echo $1 | tr ',' ' '`)
codes=/mnt/share05/data/product_raw_data/rawdata/script/cancer_code.list
for batch in "${array[@]}"
do
  query grep $batch | grep -w -f $codes | grep -E -v 'ZRNDZL|ZRnDZL|NHCeF|HLA-|mNGS-seq|NHC' >> batch_tmp.txt
done
if [ $2 ]
then
  # 把英文逗号替换成中文逗号，否则后续脚本报错！
  awk -F '\t' '{OFS=","}{$1=$1;gsub(/,/, "，", $12);print $3,$15,$4,$5,$2,$7,$12}' batch_tmp.txt | sort -k 2,2 -k 1n,1 -t,
else
  awk -F '\t' '{OFS=","}{$1=$1;print $3,$15,$4,$5,$2,$7}' batch_tmp.txt | sort -k 2,2 -k 1n,1 -t,
fi
rm batch_tmp.txt
