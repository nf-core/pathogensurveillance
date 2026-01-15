## 2025-07-16

I tried to set the `maxForks` based on `workflow.executor` but it seems that `workflow.executor` is not available when the configuration files are read and `maxForks` cannot be set using a closure for deferred execution.

## 2025-07-25

On the OSU CQLS Slum HPC, I was getting the following error from multiqc:

.command.sh: line 11: 41 Illegal instruction (core dumped) multiqc --force --outdir all_multiqc --config multiqc_config.yml .

This happens for all tested data on multiqc 1.30 and 1.29, but not 1.28, so reverting to 1.28 for now.
