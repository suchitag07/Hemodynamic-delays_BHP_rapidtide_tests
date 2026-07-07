#!/bin/bash

for subject in "$@"; do

echo "Processing subject: $subject"
echo "HOSTNAME: $(hostname)" 

subjID=${subject//_/}

mkdir -p "/Outputs/main/rapidtide_3.1.11/timepoint_1/defaultfilt/sub-${subjID}"
chmod -R 775 "/Outputs/main/rapidtide_3.1.11/timepoint_1/defaultfilt/sub-${subjID}"

CMD="apptainer exec \
  -B /Outputs/main/rsfMRI_Preproc_Derivatives/timepoint_1/fmriprep/sub-${subjID}/ses-01:/data_in \
  -B /Outputs/main/rapidtide_3.1.11/timepoint_1/defaultfilt/sub-${subjID}:/data_out \
  /Software/sif/rapidtide_latest_3.1.11.sif rapidtide \
  	/data_in/func/sub-${subjID}_ses-01_task-rest_desc-preproc_bold.nii.gz\
    /data_out/sub-${subjID} \
    --numnull 10000 \
	--filterband lfo \
	--preppass \
	--sharpenregressor \
	--searchrange -5 30 \
	--passes 3 \
	--nofitfilt \
	--simcalcrange 130 -1 \
	--corrmask /data_in/func/sub-${subjID}_ses-01_task-rest_desc-brain_mask.nii.gz \
	--globalmeaninclude /data_in/fMRIPrep_parc/bold_space/sub-${subjID}_ses-01_task-rest_space-bold_desc-aparcaseg.nii.gz:8,47 \
	--refineinclude /data_in/fMRIPrep_parc/bold_space/sub-${subjID}_ses-01_task-rest_space-bold_desc-aparcaseg.nii.gz:8,47 \
	--whitemattermask /data_in/fMRIPrep_parc/bold_space/sub-${subjID}_ses-01_task-rest_space-bold_desc-WM_probseg_bin.nii.gz \
	--csfmask /data_in/fMRIPrep_parc/bold_space/sub-${subjID}_ses-01_task-rest_space-bold_desc-CSF_probseg_bin.nii.gz \
	--motionfile /data_in/func/sub-${subjID}_ses-01_task-rest_desc-confounds-trimmed_timeseries_revised_final.tsv \
	--motpowers 2 \
	--mklthreads 1 \
	--nprocs 1 \
	--numskip 0 \
	--outputlevel max"
        
echo -e "${CMD}\n\n"

eval ${CMD}

done