# ECPred Release Notes

## v1.3.5 (2025-12-29)
- Fix: Correct `lib/EC` path joining using `Paths.get(args[2], "lib", "EC")` in `ECPred.java`.
- Fix: FASTA whitespace removal regex corrected to `"\\s+"`.
- Rebuilt `ECPred.jar` with the above fixes and recent code quality improvements.

## v1.3.4 (2025-10-13)
- Perf: Optimize CPU utilization during EC predictions.
- Use full available CPUs for class-level processing; BLAST threads set to `userCpuCount/6`.
- Dynamic thread pool management via `ParallelExecutor`.

## v1.3.3 (2025-10-13)
- Fix: Robust handling of short protein sequences.
- Generate zero/default feature vectors for unprocessable proteins to avoid index mismatches.
- Prevent `ArrayIndexOutOfBoundsException` in hierarchical predictions; ensure output includes all inputs.

## v1.3.2 (2025-10-13)
- Maintenance: JAR rebuild and cleanup; removed outdated docs.

## v1.3.1
- Fix: Correct non-enzyme confidence calculation.

## v1.3.0
- Feature: In-JVM SVMlight classifier; major speed-up and parallel-safety improvements.

## Notes
- Usage supports an optional 6th argument for `threads` to control CPU usage.
- Run from the parent directory where `lib/` exists (i.e., pass `libraryDir` such that `libraryDir/lib/EC` and `libraryDir/subclasses` are accessible).
- Tested on Linux with Java 17.
