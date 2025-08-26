#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "body.h"

// Constants
const double G = 6.67430e-11;
const double THETA = 0.5;
const double DT = 0.01;

// --- CPU-Side Tree Representation for GPU ---
// We flatten the tree into arrays for efficient GPU transfer and access.
struct GPUNode {
    double com_x, com_y; // Center of Mass
    double mass;
    double width;
    int is_leaf;
    int child_index; // Index of the first child in the nodes array. -1 if no children.
                     // Children are stored consecutively: [child_index, child_index+1, child_index+2, child_index+3]
};
std::vector<GPUNode> gpu_nodes; // This vector will hold the flattened tree

// --- CUDA Kernel ---
// This kernel calculates the force on each body in parallel.
__global__ void calculateForcesKernel(Body* bodies, GPUNode* nodes, int num_nodes, int num_bodies, double theta) {
    // Each thread handles one body
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_bodies) return;

    Body* current_body = &bodies[idx];
    current_body->acc_x = 0;
    current_body->acc_y = 0;

    // We traverse the tree for this body.
    // We use a stack to avoid recursion. Alternatively, we could use a while loop.
    int stack[64]; // Stack for node indices to process. 64 is a safe depth.
    int stack_ptr = 0;
    stack[stack_ptr++] = 0; // Push the root node (index 0)

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr]; // Pop the next node
        GPUNode node = nodes[node_idx];

        double dx = node.com_x - current_body->pos_x;
        double dy = node.com_y - current_body->pos_y;
        double r_squared = dx*dx + dy*dy;
        // Avoid division by zero and self-interaction
        if (r_squared < 1e-10) continue;

        double r = sqrt(r_squared);

        // Check Barnes-Hut condition: (s/d < θ)
        if ((node.width / r < theta) || node.is_leaf) {
            // If condition met or it's a leaf, calculate force from this node's COM
            double force_mag = G * current_body->mass * node.mass / r_squared;
            current_body->acc_x += force_mag * dx / r;
            current_body->acc_y += force_mag * dy / r;
        } else {
            // Otherwise, push all children onto the stack to process later
            if (node.child_index != -1) {
                for (int i = 0; i < 4; ++i) {
                    stack[stack_ptr++] = node.child_index + i;
                }
            }
        }
    }
}

// --- CPU Functions ---
void generateRandomBodies(std::vector<Body>& bodies, int N, double region_size) {
    // ... (Same as sequential version) ...
}
// TODO: You need a function to build the CPU tree (like in the sequential version)
// and then a NEW function to convert that tree into the flattened `gpu_nodes` vector.
// This is a non-trivial step requiring a tree traversal.

int main() {
    const int N = 10000;
    const double REGION_SIZE = 1000.0;
    const int STEPS = 10;
    const int BLOCK_SIZE = 256;

    std::vector<Body> bodies;
    generateRandomBodies(bodies, N, REGION_SIZE);

    // GPU Pointers
    Body* d_bodies = nullptr;
    GPUNode* d_nodes = nullptr;

    // Allocate GPU memory for bodies
    cudaMalloc(&d_bodies, N * sizeof(Body));
    cudaMemcpy(d_bodies, bodies.data(), N * sizeof(Body), cudaMemcpyHostToDevice);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < STEPS; ++step) {
        // 1. & 2. ON CPU: Build Tree and convert it to flattened `gpu_nodes` vector
        // TODO: This is your core algorithmic work from Phase 1, plus the flattening step.
        // buildTreeAndFlatten(bodies, REGION_SIZE); 
        std::cerr << "Tree building and flattening not implemented yet!" << std::endl;

        // 3. Copy the flattened tree to GPU
        cudaMalloc(&d_nodes, gpu_nodes.size() * sizeof(GPUNode));
        cudaMemcpy(d_nodes, gpu_nodes.data(), gpu_nodes.size() * sizeof(GPUNode), cudaMemcpyHostToDevice);

        // 4. Launch Kernel to calculate forces on GPU
        int grid_size = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
        calculateForcesKernel<<<grid_size, BLOCK_SIZE>>>(d_bodies, d_nodes, gpu_nodes.size(), N, THETA);
        cudaDeviceSynchronize(); // Wait for kernel to finish

        // 5. Copy updated bodies back to CPU to update positions
        cudaMemcpy(bodies.data(), d_bodies, N * sizeof(Body), cudaMemcpyDeviceToHost);

        // 6. Update Positions and Velocities (on CPU)
        for (Body& b : bodies) {
            b.vel_x += b.acc_x * DT;
            b.vel_y += b.acc_y * DT;
            b.pos_x += b.vel_x * DT;
            b.pos_y += b.vel_y * DT;
        }

        // 7. Copy the updated positions back to GPU for the next iteration
        cudaMemcpy(d_bodies, bodies.data(), N * sizeof(Body), cudaMemcpyHostToDevice);

        // 8. Free the GPU memory for the tree for this step
        cudaFree(d_nodes);
        gpu_nodes.clear(); // Clear the CPU-side tree vector

        if (step % 10 == 0) std::cout << "Step " << step << std::endl;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "CUDA Simulation took " << elapsed.count() << " seconds." << std::endl;

    // Cleanup
    cudaFree(d_bodies);
    return 0;
}