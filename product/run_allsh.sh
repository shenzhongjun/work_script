if [[ ! -f "*_result.sh" && ! -f "*_lims_done.sh" ]]
then
	ls | xargs -i sh {}
fi
