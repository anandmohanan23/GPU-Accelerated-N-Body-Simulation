# GPU-Accelerated N-Body Simulation with the Barnes-Hut Algorithm
## 👨‍💻 Team Member: Anand Mohanan (CB.AI.U4AID23103) B-17

## 📌 PROJECT OVERVIEW:

This project tackles the challenge of simulating the gravitational interactions between a large number of particles (like stars in a galaxy).

**The Core Problem**: Directly calculating the gravitational force between every pair of N bodies is an O(N 
2
 ) problem, which becomes incredibly slow as the number of bodies increases.

**The Solution**: We implement the Barnes-Hut algorithm, a clever approximation technique that reduces the computational complexity to a much more manageable O(NlogN). This is achieved by grouping distant particles together and treating them as a single, combined mass.

**The Goal**: The project's main objective is to demonstrate the power of GPU acceleration. We first build a standard sequential version in C++ to run on the CPU. Then, we identify the most computationally intensive part—the force calculation—and offload it to the GPU using CUDA.

**The Outcome**: By benchmarking both the CPU and GPU versions with an increasing number of bodies, we will analyze and visualize the performance gains, calculating the speedup to quantify the massive advantage of parallel processing on the GPU for this type of scientific simulation.
___
## ⚙️FEATURES:

**Problem Solved**: The project addresses the computational infeasibility of simulating the gravitational interaction of a large number of bodies using a naive O(N^2)algorithm. It uses the Barnes-Hut approximation algorithm to reduce the complexity to O(NlogN).

**Methodology**: The project involves implementing the Barnes-Hut algorithm first on a CPU sequentially and then parallelizing the most computationally intensive parts on a GPU using CUDA.


**Key Deliverable**: The final deliverable is a comparative analysis of the performance and speedup of the CUDA version over the sequential CPU implementation, for varying numbers of bodies (N).
___
