import java.util.concurrent.*;
import java.util.List;
import java.util.ArrayList;
import java.util.function.Supplier;

/**
 * Utility class for managing parallel execution in ECPred
 * Provides thread pools for different levels of parallelization:
 * - Protein level
 * - Method level (BLAST, SPMAP, PEPSTATS)
 * - EC class level
 */
public class ParallelExecutor {
    
    // Singleton instance
    private static ParallelExecutor instance;
    
    // Thread pools for different levels
    private final ExecutorService proteinLevelExecutor;
    private final ExecutorService methodLevelExecutor;
    private final ExecutorService ecClassLevelExecutor;
    
    // Configuration
    private static int userCpuCount = Runtime.getRuntime().availableProcessors();
    private static int proteinLevelThreads = userCpuCount;
    private static final int METHOD_LEVEL_THREADS = 3; // BLAST, SPMAP, PEPSTATS
    private static int ecClassLevelThreads = Math.max(1, userCpuCount / 2);
    
    private ParallelExecutor() {
        // Set thread counts from current userCpuCount
        proteinLevelThreads = userCpuCount;
    ecClassLevelThreads = Math.max(1, userCpuCount / 2);

        // Create thread pools with appropriate sizes
        this.proteinLevelExecutor = Executors.newFixedThreadPool(
            proteinLevelThreads,
            new ThreadFactory() {
                private int counter = 0;
                public Thread newThread(Runnable r) {
                    Thread t = new Thread(() -> {
                        r.run();
                    }, "ProteinWorker-" + counter++);
                    t.setDaemon(false);
                    return t;
                }
            }
        );

        this.methodLevelExecutor = Executors.newFixedThreadPool(
            METHOD_LEVEL_THREADS,
            new ThreadFactory() {
                private int counter = 0;
                public Thread newThread(Runnable r) {
                    Thread t = new Thread(() -> {
                        r.run();
                    }, "MethodWorker-" + counter++);
                    t.setDaemon(false);
                    return t;
                }
            }
        );

        this.ecClassLevelExecutor = Executors.newFixedThreadPool(
            ecClassLevelThreads,
            new ThreadFactory() {
                private int counter = 0;
                public Thread newThread(Runnable r) {
                    Thread t = new Thread(() -> {
                        r.run();
                    }, "ECClassWorker-" + counter++);
                    t.setDaemon(false);
                    return t;
                }
            }
        );
    }
    
    /**
     * Get singleton instance
     */
    public static synchronized ParallelExecutor getInstance() {
        if (instance == null) {
            instance = new ParallelExecutor();
        }
        return instance;
    }
    
    /**
     * Execute tasks in parallel at protein level
     */
    public <T> List<T> executeProteinLevel(List<Callable<T>> tasks) throws InterruptedException, ExecutionException {
        List<Future<T>> futures = proteinLevelExecutor.invokeAll(tasks);
        List<T> results = new ArrayList<>();
        for (Future<T> future : futures) {
            results.add(future.get());
        }
        return results;
    }
    
    /**
     * Execute tasks in parallel at method level (BLAST, SPMAP, PEPSTATS)
     */
    public <T> List<T> executeMethodLevel(List<Callable<T>> tasks) throws InterruptedException, ExecutionException {
        List<Future<T>> futures = methodLevelExecutor.invokeAll(tasks);
        List<T> results = new ArrayList<>();
        for (Future<T> future : futures) {
            results.add(future.get());
        }
        return results;
    }
    
    /**
     * Execute tasks in parallel at EC class level
     */
    public <T> List<T> executeECClassLevel(List<Callable<T>> tasks) throws InterruptedException, ExecutionException {
        List<Future<T>> futures = ecClassLevelExecutor.invokeAll(tasks);
        List<T> results = new ArrayList<>();
        for (Future<T> future : futures) {
            results.add(future.get());
        }
        return results;
    }
    
    /**
     * Submit single task at protein level
     */
    public <T> CompletableFuture<T> submitProteinTask(Supplier<T> task) {
        return CompletableFuture.supplyAsync(task, proteinLevelExecutor);
    }
    
    /**
     * Submit single task at method level
     */
    public <T> CompletableFuture<T> submitMethodTask(Supplier<T> task) {
        return CompletableFuture.supplyAsync(task, methodLevelExecutor);
    }
    
    /**
     * Submit single task at EC class level
     */
    public <T> CompletableFuture<T> submitECClassTask(Supplier<T> task) {
        return CompletableFuture.supplyAsync(task, ecClassLevelExecutor);
    }
    
    /**
     * Wait for all CompletableFutures to complete
     */
    public static <T> CompletableFuture<List<T>> allOf(List<CompletableFuture<T>> futures) {
        return CompletableFuture.allOf(futures.toArray(new CompletableFuture[0]))
            .thenApply(v -> {
                List<T> results = new ArrayList<>();
                for (CompletableFuture<T> future : futures) {
                    results.add(future.join());
                }
                return results;
            });
    }
    
    /**
     * Shutdown all executors gracefully
     */
    public void shutdown() {
        shutdownExecutor(proteinLevelExecutor, "Protein Level");
        shutdownExecutor(methodLevelExecutor, "Method Level");
        shutdownExecutor(ecClassLevelExecutor, "EC Class Level");
    }
    
    private void shutdownExecutor(ExecutorService executor, String name) {
        executor.shutdown();
        try {
            if (!executor.awaitTermination(60, TimeUnit.SECONDS)) {
                System.err.println("Warning: " + name + " executor did not terminate in time");
                executor.shutdownNow();
                if (!executor.awaitTermination(60, TimeUnit.SECONDS)) {
                    System.err.println("Error: " + name + " executor did not terminate");
                }
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }
    
    /**
     * Get executor for protein level tasks
     */
    public ExecutorService getProteinLevelExecutor() {
        return proteinLevelExecutor;
    }
    
    /**
     * Get executor for method level tasks
     */
    public ExecutorService getMethodLevelExecutor() {
        return methodLevelExecutor;
    }
    
    /**
     * Get executor for EC class level tasks
     */
    public ExecutorService getECClassLevelExecutor() {
        return ecClassLevelExecutor;
    }
    
    /**
     * Set user-provided CPU count for thread pool sizing
     * Must be called before getInstance() to take effect.
     */
    public static void setUserCpuCount(int cpuCount) {
        userCpuCount = cpuCount;
        proteinLevelThreads = userCpuCount;
        ecClassLevelThreads = Math.max(1, userCpuCount / 2);
    }
}
