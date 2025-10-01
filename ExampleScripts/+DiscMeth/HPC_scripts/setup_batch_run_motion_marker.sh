#!/usr/bin/env python

import subprocess

subjectNames = ['Subject_02'] #, 'Subject_03', 'Subject_04', 'Subject_05', 'Subject_06', 'Subject_07', 'Subject_08', 'Subject_09', 'Subject_10', 'Subject_11']
motions = ['straightslowrunning'] # 'curvedslowrunning', 'vcut']
trials = ['trial0025', 'trial0100', 'trial0122']  # Subject_02
        #  ['trial0028', 'trial0075', 'trial0089'],  # Subject_03
        #  ['trial0020', 'trial0063', 'trial0079'],  # Subject_04
        #  ['trial0031', 'trial0061', 'trial0084'],  # Subject_05
        #  ['trial0025', 'trial0076', 'trial0103'],  # Subject_06
        #  ['trial0026', 'trial0076', 'trial0091'],  # Subject_07
        #  ['trial0025', 'trial0061', 'trial0076'],  # Subject_08
        #  ['trial0030', 'trial0117', 'trial0152'],  # Subject_09
        #  ['trial0032', 'trial0079', 'trial0104'],  # Subject_10
        #  ['trial0026', 'trial0065', 'trial0093']]  # Subject_11

DiscAll = ['BE', 'ME', 'LIIIc-2', 'RIIa-2', 'LIIIc-3', 'RIIa-3']
nNodes = ['15', '25', '50'] # , '75', '100', '150', '200', '250', '300', '400']

nSub = len(subjectNames)
for iSub in range(len(subjectNames)):
    subject = subjectNames[iSub]
    for iTrial in range(len(motions)):
        trial = trials[iSub][iTrial]
        motion = motions[iTrial]
        for DM in DiscAll:
            for iNodes in nNodes:
                sbatch_command = """sbatch --export=ALL,SUBJECT={},MOTION={},TRIAL={},DISCMETH={},NNODES={} --job-name={}_{}_{}_{}_{} submit_run_motion_marker_3D.slurm""".format(
                    subject, motion, trial, DM, iNodes, subject, motion, trial, DM, iNodes)
                print(sbatch_command)  # Uncomment this line when testing to view the qsub command
                # Comment the following 3 lines when testing to prevent jobs from being submitted
                exit_status = subprocess.call(sbatch_command, shell=True)
                if exit_status == 1:  # Check to make sure the job submitted
                    print("Job {0} failed to submit".format(sbatch_command))

print("Done submitting jobs!")
