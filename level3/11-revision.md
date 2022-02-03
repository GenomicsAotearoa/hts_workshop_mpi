# Revision

* Teaching: 30 minutes
* Exercises: 30 minutes

#### Objectives

* Revisit the basics of working in a HPC environment
* Refamiliarise with the project directory
* Practice with simple command line tools
* Submit simple `slurm` scripts to NeSI

---

## Contents

1. [Returning to NeSI](#returning-to-nesi)
1. [xxx](#xxx)
1. [xxx](#xxx)

---

## Returning to NeSI

It has been quite a few months since the last training workshop, so the main point of the training today is to make sure everyone is up to speed on the main considerations when working with NeSI and high-performance computing. To begin, log into NeSI through the [JupyterHub web portal](https://jupyter.nesi.org.nz/) and start a session. Make sure that you use the training project (nesi03181) rather than the usual project.

Once you have logged in, open a terminal and navigate to **your** directory.

```bash
$ cd /nesi/project/nesi03181/phel/YOUR_USER_NAME/
```

We are going to do a few quick exercises to refresh on the basics of command line operations.

You should still have a folder named `shell_data/` in your directory. This is where we did our initial practice introduction to the command line and common operations. Ignoring any folders you created during your own working, there were three folders provided to you at the start of this training workshop.

> ### Exercise
>
> There are three folders in your `shell_data/` directory. One of these is hidden, the other two are not. List these folders, and report back with their names and how long ago they were created.
>
> <details>
> <summary>Solution</summary>
>
> ```bash
> $ cd shell_data/
> # ls -a -l
> drwxrws---+  2 dsen018        nesi03181 4096 Jun  2  2021 .hidden
> drwxrws---+  3 dsen018        nesi03181 4096 Jul 11  2021 sra_metadata
> drwxrws---+  3 dsen018        nesi03181 4096 Jun 23  2021 untrimmed_fastq
> ```
> </details>

---
