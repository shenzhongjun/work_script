/mnt/share01/tools/bin/perl /mnt/share02/lixx/01.script/vjob|awk '{print $12}'|cut -f7,8 -d/|sort -u|sed '1d'|tee running_jobs
wc -l running_jobs
rm running_jobs
