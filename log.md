## 2025-07-16

I tried to set the `maxForks` based on `workflow.executor` but it seems that `workflow.executor` is not available when the configuration files are read and `maxForks` cannot be set using a closure for deferred execution.
