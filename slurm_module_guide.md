## Common commands - slurm

#### Submit a `slurm` batch job

```bash
$ sbatch job_script.sh
```

#### Monitor your current jobs

```bash
$ squeue -u USERNAME
```

#### Cancel a queued or running job

```bash
$ scancel JOB_ID
```

#### View the status of completed jobs

```bash
# A particular job
$ sacct -j JOB_ID

# All jobs run since the date YYYY-MM-DD
$ sacct -S YYYY-MM-DD

# All jobs run before the date YYYY-MM-DD
$ sacct -E YYYY-MM-DD
```

>**Note:** There are many more ways to fine tune the output of `sacct`. Refer to the [documentation](https://slurm.schedmd.com/sacct.html) for more detailed information.

#### View the efficiency statistics for a completed job

>**Note:** This command *can* be run for a job that is currently in progress, but the values will not be accurate for the full run (i.e. they are the values of the job to date, not the full job).

```bash
$ seff JOB_ID

$ nn_seff JOB_ID
```

---

## Common commands - module

#### Load a specific program

>**Note:** All modules on NeSI have version and toolchain/environment suffixes. If none is specified, the default version for the tool is loaded. The default version can be seen with the `module avail` command.

```bash
$ module load MY_TOOL
```

#### View available modules

```bash
# View all modules
$ module avail

# View all modules which match the keyword in their name
$ module avail KEYWORD

# View all modules which match the keyword in their name or description
$ module spider KEYWORD
```

#### Unload all current modules

```bash
$ module purge
```

#### Swap a currently loaded module for a different one

```bash
$ module switch CURRENT_MODULE DESIRED_MODULE
```
