# Introduction to NeSI HPC Platforms

* Teaching: 20 minutes
* Exercises: 10 minutes

#### Objectives

* Write and submit a `BLAST` job to the NeSI compute cluster using the `slurm` system
* Understand the compatibility considerations that must be made when loading modules

#### Keypoints

* NeSI uses a login node (`lander`) to submit large computational tasks to the `mahuika` and `maui` clusters.
* Computational resources are requested and provisioned using the `slurm` management system.
* Pre-installed software can be accessed using the `module load` command.

---

## Contents

1. [Introduction to NeSI HPC Platforms](#introduction-to-nesi-hpc-platforms)
1. [Searching and load software with the `module` system](#searching-and-load-software-with-the-module-system)
1. [Considerations when working with modules](#considerations-when-working-with-modules)
---

## Introduction to NeSI HPC Platforms

This section will be delivered through the a presentation. The slides are available [here](Introduction_to_NeSI_HPCplatforms.pdf).

---

## Searching and load software with the `module` system

In the `slurm` example we just worked through, we were provided the software names and versions to a piece of software (`BLASTn`) and a database which were both already installed on NeSI. But what do we do if we don't know what is available, or what version it is?

Before we worry about loading software, there are a few handy commands and resources for finding existing software on NeSI. The most user-friendly option is the NeSI [Supported Applications](https://support.nesi.org.nz/hc/en-gb/sections/360000040076-Supported-Applications) web page which provides an up to date list of everything available on NeSI. Each piece of software lists the versions installed and links to the software documentation, and every entry is tagged with some handy keywords to enable quick searching.

Alternatively, if we are already logged into NeSI then we can search from the command line to find software relevant to us. To view our currently loaded software modules, we can use the `module list` command.

If you run this from a fresh log in from the JupyterHub portal you will see an output similar to:

```bash
$ module list

Currently Loaded Modules:
  1) NeSI                        (S)  11) gimkl/2020a                        21) libxml2/2.9.10-GCCcore-9.2.0               31) Tk/8.6.10-GCCcore-9.2.0
  2) slurm                            12) bzip2/1.0.8-GCCcore-9.2.0          22) libxslt/1.1.34-GCCcore-9.2.0               32) Python/3.8.2-gimkl-2020a
  3) GCCcore/9.2.0                    13) XZ/5.2.4-GCCcore-9.2.0             23) cURL/7.64.0-GCCcore-9.2.0                  33) nodejs/14.16.1-GCCcore-9.2.0
  4) zlib/1.2.11-GCCcore-9.2.0        14) libpng/1.6.37-GCCcore-9.2.0        24) PCRE/8.43-GCCcore-9.2.0                    34) git/2.23.3
  5) binutils/2.32-GCCcore-9.2.0      15) freetype/2.10.1-GCCcore-9.2.0      25) netCDF/4.7.3-gimpi-2020a                   35) JupyterLab/.2021.5.0-gimkl-2020a-3.0.15 (H)
  6) GCC/9.2.0                        16) Szip/2.1.1-GCCcore-9.2.0           26) SQLite/3.31.1-GCCcore-9.2.0                36) craype-broadwell
  7) libpmi                           17) HDF5/1.10.5-gimpi-2020a            27) METIS/5.1.0-GCCcore-9.2.0                  37) craype-network-infiniband
  8) impi/2019.6.166-GCC-9.2.0        18) libjpeg-turbo/2.0.2-GCCcore-9.2.0  28) tbb/2019_U9-GCCcore-9.2.0
  9) gimpi/2020a                      19) ncurses/6.1-GCCcore-9.2.0          29) SuiteSparse/5.6.0-gimkl-2020a-METIS-5.1.0
 10) imkl/2020.0.166-gimpi-2020a      20) libreadline/8.0-GCCcore-9.2.0      30) Tcl/8.6.10-GCCcore-9.2.0

  Where:
   H:  Hidden Module
```

Alternatively, if you have logged in through a basic `ssh` connection, you will see a much smaller output:

```bash
$ module list

Currently Loaded Modules:
  1) slurm   2) NeSI (S)
```

If we want to see what is available to load, there are two options. The first is to simply report a list of every software module available on NeSI:

```bash
$ module avail

$ module avail blast
```

In the first instance, we see everything available on NeSI. In the second, we see everything available on NeSI with the keyword 'blast' in the name. There is also a more thorough search option:

```bash
$ module spider blast
```

The difference here is that `avail` searches the module names, and `spider` searches their description and other information for the keyword. Both are useful, and in both cases the search is case insensitive. For example, above we used the lowercase spelling of `BLAST` but in the results we have a mixture of cases (`BLAST/2.3.0`, `RMBlast/2.6.0-gimkl-2017a`, `samblaster/0.1.24-gimkl-2017a`).

When we want to go and load a module, the `module load` command **_is case sensitive_** so we must use the exact result from `module avail` or `module spider`.

```bash
$ module load blast/2.3.0
Lmod has detected the following error:  The following module(s) are unknown: "blast/2.3.0"

Please check the spelling or version number. Also try "module spider ..."
It is also possible your cache file is out-of-date try:
  $ module --ignore-cache load "blast/2.3.0"

$ module load BLAST/2.3.0
```

When we load a module, if there is no feedback from the command prompt then the load was successful. We can now access our tool from the command line like we can for the native tools like `grep`, `ls`, and `cp`.

---

## Considerations when working with modules

Software modules are used for a number of reasons. The first is that whenever a session starts (i.e. you log into NeSI) all required tools must be found and loaded by the operating system. With the number of tools available on NeSI this would be prohibitive.

The second reason is that software versions change through time, as new features are added or bugs are fixed. For the `samtools` software, which we will use later in this training program, there are currently 7 versions of the software installed into NeSI (0.1.18, 0.1.19, 1.3.1, 1.8, 1.9, 1.10, 1.12). All of these are executed through the `samtools` command, so if all were simultaneously available NeSI would not know which one to use. In most situations we would want to be working with the most recent version of the software, and if we use the `module load` command without providing a version number then this is what will load, but sometimes there are reasons to use an older version.

```bash
$ module avail blast
BLAST/2.3.0
BLAST/2.6.0-gimkl-2017a
BLAST/2.6.0-gimkl-2018b
BLAST/2.9.0-gimkl-2018
BLAST/2.10.0-GCC-9.2.0  (D)

$ module load BLAST
blastn -version
blastn: 2.10.0+
```

It is always better to load a specific version as you then have a record of which version was used should you ever need to revisit your work. This helps to avoid dependency clashes, which is the final consideration when loading modules.

When we write software, it is very rare to write the entire program from scratch. There are a wealth of publicly available resources which can be used when developing a tool which saves the developer from writing every single line of code they need to achieve their intent. This saves time and leads to more stable and robust code but can be a problem if we are trying to use tools with incompatible dependencies. An example of this issue can be found with the popular `python` scripting language. When the language was updated from version `2.7` to version `3.0` in 2006, there were fundamental changes to how the language was used so that code written in version `2.7` could not run in version `3.0+`. In order to save the global `python` community from having to rebuild all their work, version `2.7` was maintained and updated in parallel to the newer releases of `3.0+`. Unfortunately, because people didn't *need* to swap to the newer version, a situation arose in which even in 2019 tools were still being developed in the older version and users needed to implement complicated work arounds to run both `python2.7` and `python3` code in parallel as pipelines often consisted of a mixture of tools written in both versions.

For a common example of clashing dependencies, we will try to load two modules at the same time on NeSI.

> ### Exercise
>
> Try to load the `BLAST/2.6.0-gimkl-2017a` and `SAMtools/1.8-gimkl-2018b` modules in the same NeSI session. Once you have entered both load commands, you can confirm the tools are working by running them from the command line with a `--help` parameter to view the manual for the tool:
> ```bash 
> $ blast --help
> $ samtools --help
> ```
> 
> If both tools are loaded you will be able to see the help menus for both tools. If you want to unload your current tools and load in a different order, you can use the following command to reset the loaded modules.
> 
> ```bash
> $ module purge
> ```
> 
> <details>
> <summary>Solution</summary>
>
> This cannot be done. The first module can be loaded without issue, but the second one will always produce an error on load. For example;
>
> ```bash
> $ module load BLAST/2.6.0-gimkl-2017a
> $ module load SAMtools/1.8-gimkl-2018b
> Lmod has detected the following error:  These module(s) exist but cannot be loaded as requested: "gimkl/2018b"
>   Try: "module spider gimkl/2018b" to see how to load the module(s).
> ```
>
> ```bash
> $ module load SAMtools/1.8-gimkl-2018b
> $ module load BLAST/2.6.0-gimkl-2017a
>  Lmod has detected the following error:  These module(s) exist but cannot be loaded as requested: "GCCcore/5.4.0"
>   Try: "module spider GCCcore/5.4.0" to see how to load the module(s).
> ```
> </details>

This issue occurs because these two different tools were created using a different set of development tools, and both sets of code cannot be active in parallel.

This is one of the key considerations we must keep in mind when working with NeSI. Ideally, our `slurm` scripts will be minimalistic and only load a single module for each job, as in order to keep our resource usage minimal and efficient we will use resource requests tailored for the specific job. When we **_need_** to perform multiple commands in a single script we must make sure that all modules can be loaded together, or make use of the `module purge` command to isolate the software at each step of the script.

The NeSI module nomenclature works in the manner of `SOFTWARE/VERSION-BUILD`, for the `BLASTn` example above we can see this as:

```
BLAST/2.6.0-gimkl-2017a
|-----|-----|----------
Tool  Ver.  BUILD
```

If tools have the same build description then they should be able to load together. There are also cases where tools with different build types can still load together. For example:

```bash
$ module load BLAST/2.10.0-GCC-9.2.0
$ module load SAMtools/0.1.19-GCCcore-9.2.0
```

Works without issue, as the two build environments are compatible. **Always test that your modules are compatible by loading them from your JupyterHub terminal before loading them in a script.**

If you absolutely must use incompatible versions of tools in the same script, separate them using the `module purge` command to unload the first version prior to loading the second.

```bash
$ module load BLAST/2.6.0-gimkl-2017a
# Do some BLAST work...

$ module purge
$ module load SAMtools/1.8-gimkl-2018b
# Do some samtools work...
```

---

[Next lesson](04-working-with-illumina-data.md)
