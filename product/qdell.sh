if [ $# == 0 ];then
	echo "need a keyword!"
	exit 1
fi

keyword=$1
read -p "Are you really want to qdel jobs of keyword: \"${keyword}\"?
`/mnt/share01/tools/bin/perl /mnt/share02/lixx/01.script/vjob | grep ${keyword} | awk '{print $1,$2,$3,$4,$9,$12}' | tr ' ' '\t'`
y/n": name

if [ "$name" = "y" ];then
	/mnt/share01/tools/bin/perl /mnt/share02/lixx/01.script/vjob | grep ${keyword} | cut -f1 -d. | xargs -i qdel {}
else
	echo "nothing happened."
fi
