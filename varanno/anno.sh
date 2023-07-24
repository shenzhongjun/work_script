FILE=table
INPUT=`pwd`/table
OUTPUT=`pwd`/anno_result.txt

/mnt/share01/tools/bin/perl /mnt/share03/tumor_pipeline/Somatic/DNA/WES_pipe_v4/bin/VarAnnotation_v1.3.pl \
	-f $FILE \
	-i $INPUT \
	-o $OUTPUT
