# 1.4 - The NeSI computing environment

## Overview

!!! clock "time"

    * Teaching: 30 minutes

!!! circle-info "Objectives and Key points"

    #### Objectives
    
    * Understand the differences between the `project/` and `nobackup/` directories.
    * Understand the key commands when working with the module system.
    * Understand the key commands when working with `slurm`.

    #### Keypoints
    
    * NeSI provide two data directories for performing bioinformatic work. Data should be kept in the location most appropriate for the task.
    * The module system is used to load and unload software for specific computing tasks.
    * Software versioning can be controlled via the module system.
    * Jobs are submitted to the cluster via the `slurm`, where they run asynchronously from the users session.
    
---

## Persistant and temporary storage locations

When working on NeSI there are two locations we can use to store our data; the `project/` directory (where we work for these training workshops) and the `nobackup/` location. These file locations are similar in terms of how they are accessed:

!!! terminal "code"

    ```bash
    cd /nesi/project/nesi03181/

    cd /nesi/nobackup/nesi03181/
    ```

There are important differences between them. Data on the `project/` directory is *persistent*, backed up, and a billed resource. Anything you place in the `project/` directory will stay there until it is deleted.

In constrast, files written to the `nobackup/` directory do not last forever. Files that have not been modified (a proxy for 'used') in the last three months are automatically removed from the system. The benefit of this directory though is that we are not charged for our usage of the `nobackup/` directory.

There are many tools we work with which create a lot of intermediate files as they run, and realisatically we do not need to keep everything that produce. In particular, assembly tools such as `SPAdes` and `Canu` create a lot of temporary files as they move through rounds of data cleaning, assembly, and refinement but in the end only a handful of the outputs are actually useful for the end user.

When working with large projects, it is best practice to keep your raw data in the `project/` directory but perform your day to day work in the `nobackup/` side. Only copy the critical outputs back to the `project/` side.

---

## Search and load software with the `module` system

Up until this point we have worked exclusively with software that is avilable by default in `bash` environments. For biology-related tasks though, these tools alone are not sufficient and we must use specialised software for performing analyses.

While we are free to install new software into NeSI there are many common bioinformatics tools already available on the platform, we just need to know how to access them. The most user-friendly option for finding pre-installed software is the NeSI [Supported Applications](https://support.nesi.org.nz/hc/en-gb/sections/360000040076-Supported-Applications) web page. This provides an up to date list of everything available on NeSI. Each piece of software lists the versions installed and links to the software documentation, and every entry is tagged with some handy keywords to enable quick searching.

Alternatively, if we are already logged into NeSI then we can search from the command line to find software relevant to us. To view our currently loaded software modules, we can use the `module list` command.

!!! terminal "code"

    ```bash
    module list
    ```

??? success "Output"

    ```
    Currently Loaded Modules:
    1) XALT/minimal                      15) imkl-FFTW/2022.0.2-gimpi-2022a      29) libxml2/2.9.10-GCCcore-11.3.0
    2) NeSI                         (S)  16) gimkl/2022a                         30) libxslt/1.1.34-GCCcore-11.3.0
    3) slurm                             17) nodejs/16.15.1-GCCcore-11.3.0       31) cURL/7.83.1-GCCcore-11.3.0
    4) GCCcore/11.3.0                    18) git/2.23.3                          32) netCDF/4.8.1-gimpi-2022a
    5) zlib/1.2.11-GCCcore-11.3.0        19) ZeroMQ/4.3.4-GCCcore-11.3.0         33) SQLite/3.36.0-GCCcore-11.3.0
    6) binutils/2.38-GCCcore-11.3.0      20) bzip2/1.0.8-GCCcore-11.3.0          34) Tcl/8.6.10-GCCcore-11.3.0
    7) GCC/11.3.0                        21) XZ/5.2.5-GCCcore-11.3.0             35) Tk/8.6.10-GCCcore-11.3.0
    8) libpmi/2-slurm                    22) libpng/1.6.37-GCCcore-11.3.0        36) OpenSSL/1.1.1k-GCCcore-11.3.0
    9) numactl/2.0.14-GCC-11.3.0         23) freetype/2.11.1-GCCcore-11.3.0      37) Python/3.10.5-gimkl-2022a
    10) UCX/1.12.1-GCC-11.3.0             24) Szip/2.1.1-GCCcore-11.3.0           38) JupyterLab/.2023.1.0-gimkl-2022a-3.5.3 (H)
    11) impi/2021.5.1-GCC-11.3.0          25) HDF5/1.12.2-gimpi-2022a             39) craype-broadwell
    12) AlwaysIntelMKL/1.0                26) libjpeg-turbo/2.1.3-GCCcore-11.3.0  40) craype-network-infiniband
    13) imkl/2022.0.2                     27) ncurses/6.2-GCCcore-11.3.0
    14) gimpi/2022a                       28) libreadline/8.1-GCCcore-11.3.0

    Where:
    H:  Hidden Module
    ```

If we want to see what additional software is available to load, there are two options. The first is to simply report a list of every software module available on NeSI:

!!! terminal "code"

    ```bash
    # See all tools...
    module avail

    # See all tools that match a given keyword...
    module avail blast
    ```

??? success "Output (blast version)"

    ```
    ----------------------------------------------- /opt/nesi/CS400_centos7_bdw/modules/all -----------------------------------------------
    BLAST/2.3.0                BLAST/2.12.0-GCC-9.2.0         BLASTDB/2023-04                  samblaster/0.1.24-gimkl-2017a
    BLAST/2.6.0-gimkl-2017a    BLAST/2.13.0-GCC-11.3.0 (D)    BLASTDB/2023-07           (D)    samblaster/0.1.26-GCC-9.2.0   (D)
    BLAST/2.6.0-gimkl-2018b    BLASTDB/2022-07                RMBlast/2.6.0-gimkl-2017a
    BLAST/2.9.0-gimkl-2018b    BLASTDB/2022-10                RMBlast/2.9.0-GCC-7.4.0
    BLAST/2.10.0-GCC-9.2.0     BLASTDB/2023-01                RMBlast/2.10.0-GCC-9.2.0  (D)

    Where:
    D:  Default Module

    Use "module spider" to find all possible modules.
    Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".
    ```

In the first instance, we see everything available on NeSI. In the second, we see everything available on NeSI with the keyword 'blast' in the name. There is also a more thorough search option:

!!! terminal "code"

    ```bash
    module spider blast
    ```

??? success "Output"

    ```
    -----------------------------------------------------------------------------------------------------------------------------------
    BLAST:
    -----------------------------------------------------------------------------------------------------------------------------------
        Description:
        Basic Local Alignment Search Tool, or BLAST, is an algorithm for comparing primary biological sequence information, such as
        the amino-acid sequences of different proteins or the nucleotides of DNA sequences. 

        Versions:
            BLAST/2.3.0
            BLAST/2.6.0-gimkl-2017a
            BLAST/2.6.0-gimkl-2018b
            BLAST/2.9.0-gimkl-2018b
            BLAST/2.10.0-GCC-9.2.0
            BLAST/2.12.0-GCC-9.2.0
            BLAST/2.13.0-GCC-11.3.0

    -----------------------------------------------------------------------------------------------------------------------------------
    For detailed information about a specific "BLAST" module (including how to load the modules) use the module's full name.
    For example:

        $ module spider BLAST/2.9.0-gimkl-2018b
    -----------------------------------------------------------------------------------------------------------------------------------

    -----------------------------------------------------------------------------------------------------------------------------------
    BLASTDB:
    -----------------------------------------------------------------------------------------------------------------------------------
        Description:
        BLAST databases downloaded from NCBI.

        Versions:
            BLASTDB/2022-07
            BLASTDB/2022-10
            BLASTDB/2023-01
            BLASTDB/2023-04
            BLASTDB/2023-07

    -----------------------------------------------------------------------------------------------------------------------------------
    For detailed information about a specific "BLASTDB" module (including how to load the modules) use the module's full name.
    For example:

        $ module spider BLASTDB/2023-07
    -----------------------------------------------------------------------------------------------------------------------------------

    -----------------------------------------------------------------------------------------------------------------------------------
    RMBlast:
    -----------------------------------------------------------------------------------------------------------------------------------
        Description:
        RMBlast supports RepeatMasker searches by adding a few necessary features to the stock NCBI blastn program. These include:
        Support for custom matrices ( without KA-Statistics ). Support for cross_match-like complexity adjusted scoring. Cross_match
        is Phil Green's seeded smith-waterman search algorithm. Support for cross_match-like masklevel filtering.. 

        Versions:
            RMBlast/2.6.0-gimkl-2017a
            RMBlast/2.9.0-GCC-7.4.0
            RMBlast/2.10.0-GCC-9.2.0

    -----------------------------------------------------------------------------------------------------------------------------------
    For detailed information about a specific "RMBlast" module (including how to load the modules) use the module's full name.
    For example:

        $ module spider RMBlast/2.9.0-GCC-7.4.0
    -----------------------------------------------------------------------------------------------------------------------------------

    -----------------------------------------------------------------------------------------------------------------------------------
    samblaster:
    -----------------------------------------------------------------------------------------------------------------------------------
        Description:
        samblaster is a fast and flexible program for marking duplicates in read-id grouped paired-end SAM files. It can also
        optionally output discordant read pairs and/or split read mappings to separate SAM files, and/or unmapped/clipped reads to a
        separate FASTQ file. When marking duplicates, samblaster will require approximately 20MB of memory per 1M read pairs. 

        Versions:
            samblaster/0.1.24-gimkl-2017a
            samblaster/0.1.26-GCC-9.2.0

    -----------------------------------------------------------------------------------------------------------------------------------
    For detailed information about a specific "samblaster" module (including how to load the modules) use the module's full name.
    For example:

        $ module spider samblaster/0.1.26-GCC-9.2.0
    -----------------------------------------------------------------------------------------------------------------------------------
    ```

The difference here is that `avail` searches the module names, and `spider` searches their description and other information for the keyword. Both are useful, and in both cases the search is case insensitive. For example, above we used the lowercase spelling of `BLAST` but in the results we have a mixture of cases (`BLAST/2.3.0`, `RMBlast/2.6.0-gimkl-2017a`, `samblaster/0.1.24-gimkl-2017a`).

When we want to go and load a module, the `module load` command **_is case sensitive_** so we must use the exact result from `module avail` or `module spider`.

!!! terminal "code"

    ```bash
    module load blast/2.3.0
    ```

??? failure "Output"

    ```
    Lmod has detected the following error:  The following module(s) are unknown: "blast/2.3.0"

    Please check the spelling or version number. Also try "module spider ..."
    It is also possible your cache file is out-of-date try:
    $ module --ignore-cache load "blast/2.3.0"
    ```

!!! terminal "code"

    ```bash
    module load BLAST/2.3.0
    ```

When we load a module, if there is no feedback from the command prompt then the load was successful. We can now access our tool from the command line like we can for the native tools like `grep`, `ls`, and `cp`.

---

## Considerations when working with modules

Software modules are used for a number of reasons. The first is that whenever a session starts (i.e. you log into NeSI) all required tools must be found and loaded by the operating system. With the number of tools available on NeSI this would be prohibitive.

The second reason is that software versions change through time, as new features are added or bugs are fixed. For the `samtools` software, which we will use later in this training program, there are currently 10 versions of the software installed into NeSI (0.1.18, 0.1.19, 1.3.1, 1.8, 1.9, 1.10, 1.12, 1.13, 1.15, 1.16). All of these are executed through the `samtools` command, so if all were simultaneously available NeSI would not know which one to use.

In most situations we would want to be working with the most recent version of the software, and if we use the `module load` command without providing a version number then this is what will load, but sometimes there are reasons to use an older version.

!!! terminal "code"

    ```bash
    module avail samtools
    ```

??? success "Output"

    ```
    ----------------------------------------------- /opt/nesi/CS400_centos7_bdw/modules/all -----------------------------------------------
    SAMtools/0.1.18-gimkl-2017a      SAMtools/0.1.19-GCCcore-11.3.0    SAMtools/1.8-gimkl-2018b    SAMtools/1.12-GCC-9.2.0
    SAMtools/0.1.18-gimkl-2018b      SAMtools/0.1.19-gimkl-2017a       SAMtools/1.9-GCC-7.4.0      SAMtools/1.13-GCC-9.2.0
    SAMtools/0.1.19-GCCcore-7.4.0    SAMtools/1.3.1-gimkl-2017a        SAMtools/1.10-GCC-9.2.0     SAMtools/1.15.1-GCC-11.3.0
    SAMtools/0.1.19-GCCcore-9.2.0    SAMtools/1.8-gimkl-2017a          SAMtools/1.11-GCC-7.4.0     SAMtools/1.16.1-GCC-11.3.0 (D)

    Where:
    D:  Default Module

    Use "module spider" to find all possible modules.
    Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".
    ```

!!! note "Be specific in your versioning"

    It is always better to load a specific version as you then have a record of which version was used should you ever need to revisit your work. This helps to avoid dependency clashes, which is the final consideration when loading modules.

---

## Conflicting dependencies between software

One major benefit of working with the module system is that it allows us to avoid conflicts between tools which were written in different programming environments. When we write software, it is very rare to write the entire program from scratch. There are a wealth of publicly available resources which can be used when developing a tool which saves the developer from writing every single line of code they need to achieve their intent. This saves time and leads to more stable and robust code but can be a problem if we are trying to use tools with incompatible dependencies.

For a common example of clashing dependencies, we can try to load the `BLAST/2.6.0-gimkl-2017a` and `SAMtools/1.8-gimkl-2018b` modules in the same NeSI session.

!!! terminal "code"

    ```bash
    module load BLAST/2.6.0-gimkl-2017a
    module load 
    ```

??? failure "Output"

    ```
    Lmod has detected the following error:  These module(s) exist but cannot be loaded as requested: "gimkl/2018b"
    Try: "module spider gimkl/2018b" to see how to load the module(s).
    ```

This cannot be done. The first module can be loaded without issue, but the second one will always produce an error on load. For example;

This issue occurs because these two different tools were created using a different set of development tools, and both sets of code cannot be active in parallel.

This is one of the key considerations we must keep in mind when working with NeSI.

Ideally, our sessions will be minimalistic and only load a single module for each job, as in order to keep our resource usage minimal and efficient we will use resource requests tailored for the specific job. When we **_need_** to perform multiple commands in a single script we must make sure that all modules can be loaded together, or make use of the `module purge` command to isolate the software at each step of the script.

!!! terminal "code"

    ```bash
    module load BLAST/2.6.0-gimkl-2017a
    # Do some BLAST work...

    module purge
    module load SAMtools/1.8-gimkl-2018b
    # Do some samtools work...
    ```

The NeSI module nomenclature works in the manner of `SOFTWARE/VERSION-BUILD`, for the `BLASTn` example above we can see this as:

!!! info ""

    ```
    BLAST/2.6.0-gimkl-2017a
    |-----|-----|----------
    Tool  Ver.  BUILD
    ```

If tools have the same build description then they should be able to load together, but **_always test that your modules are compatible by loading them from your JupyterHub terminal before loading them in a script_.**

---
