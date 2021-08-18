# Level 3.2 revision

This revision material consists of three short exercises covering the main points of the second part of Level 3 training. For each of these exercises you will need to work from your `3_Assembly_mapping/` directory. **_However_** - keep a close eye on the names of your own folders during these exercises. Each person created their folders and therefore there are several variations in names that we have seen in the workshops, such as;

* `3_Assembly_mapping/`
* `3_Assembly-mapping/`
* `3_assembly_mapping/`

It does not matter which you use, but make sure that your commands are specific to your own folders and not those of the trainers. If you find that your `3_Assembly_mapping/` folder has become fragmented across different names, it might be worth consolidating your data into a single location before proceeding.

When you have completed all exercises, email the requested material and answers to the trainers.

> Remember, if you are struggling with any exercises:
>
> 1. Every exercise here was covered in the training material. Refer to this for hints.
> 1. It is not cheating to use Google to find help on how various commands work.
> 1. You are welcome to contact the trainers for advice on a particular exercise, but please attempt the first two options before resorting to this.

---

## Exercise 1 - Genome assembly with SPAdes

Select one of the three *M. bovis* samples to assemble with `SPAdes`. Create a `slurm` script to assemble the MiSeq data, creating an output folder in your equivalent of the `3_Assembly-mapping/` folder. For this assembly, you must set the following parameters:

1. Assembly will take place in **isolate** mode.
1. Manually specify *k*-mer sizes for assembly. It is up to you which values to choose, but you must make the selection rather than let `SPAdes` pick automatically.
1. Set the number of threads to 20.
   1. *Note - you will have to account for this in your `slurm` batch file as well.

You must also provide an appropriate amount of memory and run time for the command to finish. Below as some approximations of the resources needed to assemble each genome.

|Sample|Run time (hh:mm:ss)|Memory (GB)|
|:---|:---:|:---:|
|Mb1|01:00:00|40 GB|
|Mb152|00:10:00|20 GB|
|Mb168|01:30:00|60 GB|

---

## Exercise 2 - ...

...

---

## Exercise 3 - ...

...

---
