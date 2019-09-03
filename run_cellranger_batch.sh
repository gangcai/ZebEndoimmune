#!/bin/bash
for sample in EC1dpf EC2dpf EC4dpf EC8dpf
do
	echo "start for "$sample
	date
	bash run_cellranger.sh $sample > "$sample".log 
	echo "end for "$sample
	date
done
