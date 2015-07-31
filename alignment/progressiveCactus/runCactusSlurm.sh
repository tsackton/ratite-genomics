#!/bin/bash 
/n/home12/tsackton/cuff/progressiveCactus/bin/runProgressiveCactus.sh \
	--configFile=ratite_config.xml \
	--stats \
	--logDebug \
	--logFile=ratite_align.log \
	--maxLogFileSize=100000 \
 	--batchSystem=slurm \
        --defaultMemory=6000000000 \
	--jobTime=600 \
	--maxThreads=48 \
	--bigBatchSystem=singleMachine \
	--bigMemoryThreshold=5000000000 \
        --slurm-partition=serial_requeue,general,tsackton \
        --slurm-scriptpath=slurm_scripts \
        --slurm-time=1440 \
	--retryCount=3 \
	ratite.seqFile \
	./ratiteDir \
	ratiteAlign.hal
