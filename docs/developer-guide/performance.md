# Performance Guide

High performance is a critical design goal of CATChem. This guide provides an overview of the performance features of the model, as well as best practices and strategies for writing efficient code and optimizing performance.

## Core Performance Features

CATChem's architecture is designed from the ground up for high performance on modern HPC systems.

-   **[Column Virtualization](core/column-virtualization.md)**: This is the most important performance feature of CATChem. By treating the 3D model domain as a collection of independent 1D columns, we achieve excellent cache utilization and data locality, which is key to performance on modern CPUs.

-   **Process Modularity**: The modular, process-based architecture minimizes overhead and allows for clean separation of concerns, which in turn makes it easier to optimize individual processes.

-   **Efficient Memory Management**: CATChem is designed to minimize memory allocations at runtime. We use pre-allocated buffers and a "zero-copy" approach where possible (e.g., for meteorological data in virtual columns) to reduce memory bandwidth bottlenecks.

## Profiling and Identifying Bottlenecks

**"The First Rule of Program Optimization: Don't do it. The Second Rule of Program Optimization (for experts only!): Don't do it yet." - Michael A. Jackson**

Before you attempt to optimize any code, you must first **profile** it to identify the actual bottlenecks. Time spent optimizing code that is not a bottleneck is time wasted.

### Built-in Profiling

CATChem includes a simple, built-in profiler that can be enabled in the configuration file. This profiler provides basic timing information for the main processes.

```yaml
# In your CATChem_config.yml
profiling:
  enabled: true
  level: detailed
```

This will produce a summary of the time spent in each major process, which can be a good starting point for identifying hotspots.

### External Profiling Tools

For more detailed analysis, you should use a dedicated profiling tool. Some recommended tools are:

-   **Intel VTune Profiler**: An excellent, all-around profiler for identifying hotspots, memory bandwidth issues, and more.
-   **NVIDIA Nsight Systems**: For profiling on NVIDIA GPUs.
-   **GNU gprof**: A simple, widely available profiler.

**Example with gprof:**

1.  Compile with profiling support:
    ```bash
    export FCFLAGS="-pg -O2"
    make
    ```
2.  Run the model:
    ```bash
    ./catchem_driver
    ```
3.  Analyze the output:
    ```bash
    gprof catchem_driver gmon.out > profile.txt
    ```

When analyzing profiling results, look for the functions or loops where the most time is being spent. These are the "hotspots" that are the best candidates for optimization.

## Optimization Strategies

### Compiler Optimization

Using the right compiler flags is the easiest way to get a significant performance boost. The recommended flags provide a good balance of performance and numerical stability.

-   **Intel Fortran (`ifort`)**: `-O3 -xHOST -ipo -no-prec-div -fp-model fast=2`
-   **GNU Fortran (`gfortran`)**: `-O3 -march=native -ffast-math -funroll-loops`

These flags enable aggressive optimization, vectorization, and inter-procedural optimization. Be aware that flags like `-ffast-math` can affect numerical precision, so it is important to validate your results after changing compiler flags.

### Algorithmic and Code-Level Optimization

-   **Follow the [Coding Standards](coding-standards.md)**: The coding standards, especially the guidelines on loop ordering and memory management, are designed to promote good performance.
-   **Vectorization**: Write loops in a way that the compiler can easily vectorize them. This means avoiding complex control flow inside loops and ensuring that the inner-most loop has a sufficient number of iterations.
-   **Data Structures**: Use data structures that are cache-friendly. For example, arrays of structures are often less efficient than structures of arrays.

## Parallel Performance

### OpenMP

CATChem uses OpenMP for on-node parallelism. The most common parallelization strategy is to parallelize the loop over columns.

```fortran
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do i = 1, num_columns
  ! ... process column i ...
end do
!$OMP END PARALLEL DO
```

-   **`SCHEDULE(DYNAMIC)`**: This is often a good choice for column-based processing, as the computational cost of each column can vary. Dynamic scheduling helps to balance the load across threads.
-   **Thread Affinity**: For best performance, it is important to bind threads to specific cores. This can be done with environment variables:
    ```bash
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    ```

### MPI

MPI is used for distributed-memory parallelism across multiple nodes. CATChem uses a 2D horizontal domain decomposition. For good MPI scaling, it is important to minimize communication and balance the workload across MPI ranks.

## Performance Best Practices

-   **Profile First**: Don't guess where the bottlenecks are.
-   **Optimize Hotspots**: Focus your optimization efforts on the parts of the code that are actually slow.
-   **Validate Your Changes**: Always verify that your optimizations have not introduced bugs or changed the scientific results of the model.
-   **Benchmark**: Measure the performance before and after your changes to quantify the improvement.
-   **Consider the Whole System**: A change that improves performance on one system may degrade it on another. Test on multiple platforms if possible.