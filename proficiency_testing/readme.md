# Proficiency testing for HTS workshops 

In this series of workshops we have covered a broad introduction to bioinformatics skills, learnt about a variety of different tools and applied these to answer diagnostic questions relavent to each discipline here at PHEL. 

In order to showcase how far we have all come, its time for the much anticipated final project (or proficiency test) 

A final project has been created for each team in order to closely replicate the HTS processes of that team during routine diagnostic work, for example the Entomology team will need to use sequencing data to identify an insect species through sequence analysis. 

#### Objectives

For each project you will: 

* Use the bioinformatics skills you have developed to navigate around NeSI, copy, move and interact with files to create a logically organised, easy to follow, analysis folder.
* Demonstrate your ability to use the NeSI platform to identify relevant bioinformatic tools and load them
* Write `slurm` scripts to deploy jobs to the cluster, with appropriate resource usage, to perform the necessary analysis 
* Evaluate the findings of your analysis, to complete a set of discipline-specific tasks  

The due date for your completed analysis is [TO BE DETERMINED].

---

## Housekeeping

Before beginning your work, please keep in mind the following considerations:

1. As this work is considered training, please ensure that you use the `nesi03181` project code for submitting all `slurm` jobs, not the MPI diagnostic project.
1. There is no expectation of a particualr folder structure or file layout for this assessment, but do your best to keep your work well organised and easy to follow, as you would for any other work.
1. We will be primarily using your `slurm` scripts to assess your profficiency with the tools that you have used, so please make sure that these are easy to locate.
1. During the time for this assessment, you will have **_limited_** ability to consult with the trainers for assistance in scheduled meetings. Please respect our time and use these times for questions and help.
1. Remember that this is an assessment of your ability, which includes the ability to troubleshoot issues. We will provide limited assistance but will not be fixing errors that we find.
1. Do not leave this until the last minute! Some of the jobs you will need to perform may take several hours of compute time so cannot be run on the last day. Similarly, you will need time to assess and verify your findings so make sure that you leave sufficient time to do so.

---

## Contents

Navigate to the final project folder for your discipline to learn more about your specifc assignment. If you do not belong to one of these teams, just pick the discipline or evaluation that most closely represents your own discipline.

1. [Entomology](#entomology)
1. [Mycology and Bacteriology](#mycology-and-bacteriology)
1. [Virology](#virology)

---

## Entomology

**Georgia Breckell**

Damage to a harvest of apples has been found and it is suspected insect damage. 

<img src='../img/prof_testing_crop_damage.png' alt='crop damage' width='400' />

Fragments of a suspect insect were found and the sample has been sent for molecular analysis.


<img src='../img/prof_testing_bug_leg.png' alt='Insect leg' width='350' />

For this project you will need to use the sequencing data provided here `/nesi/project/nesi03181/phel/proficiency_testing/Entomology_Test`  to identify the insect and determine if this is likely to be responsible for the damage observed. 

Notes for the assesment: 

- You do not know how these reads were produced (ie any chances for contamination or which Illumina platform was used), therefore ensure you investiagate your data and incorperate quality control at various stages of your analysis. 

---

## Mycology and Bacteriology

**Luciano Rigano**

In `/nesi/project/nesi03181/phel/proficiency_testing/Mycology_Bacteriology_Test/` you will find two directories:

1. `dataset_A_ONT/`
1. `dataset_B_illumina/`

These datasets correspond to genome sequencing data from two different bacterial organisms obtained using ONT and Illumina sequencing technologies, as indicated in the file names.

Using your preferred approach, identify the species that organism A and organism B belong to.

---

## Virology

**David Waite**

In the `/nesi/project/nesi03181/phel/proficiency_testing/Virology_Test/` folder, you will find three files, corresponding to Illumina MiSeq and MinION sequencing runs from a symtomatic *Vitis* sample. Make a copy of these files in your own directory, then use your knowledge of HTS analysis tools and workflows to identify which virus(es) are present in the sample.
  
Notes for the assessment:

1. While both of these libraries were generated from the same nucleic acid extraction, remember that there is a significant difference in the sequencing depth of each platform.
1. There are multiple viruses in the library, not all of which are diagnostically relevant. For this assessment you will need to:
   1. Identify the likely viral species in the sample.
   1. Perform additional analyses to confirm the likelihood that each detection is valid.
1. Validation can be performed through whichever means you feel is necessary, which may include your knowledge of virology (i.e. known information of host range or risk potential for each virus). However, you should be including sequence data to support your conclusions.

---
