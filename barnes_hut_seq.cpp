#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "body.h"

// Constants
const double G = 6.67430e-11; // Gravitational constant
const double THETA = 0.5;     // Barnes-Hut threshold parameter
const double DT = 0.01;       // Time step

class TreeNode {
public:
    double center_x, center_y; // Center of mass
    double total_mass;
    double width;              // Width of the square region
    bool is_leaf;
    Body* body;                // Only valid for leaf nodes
    TreeNode* children[4];     // Quadrants: NW, NE, SW, SE

    TreeNode(double x, double y, double w) : center_x(x), center_y(y), width(w), total_mass(0), is_leaf(true), body(nullptr) {
        for (int i = 0; i < 4; ++i) children[i] = nullptr;
    }

    ~TreeNode() {
        for (int i = 0; i < 4; ++i) delete children[i];
    }

    // TODO: Implement this function
    void insertBody(Body* newBody) {
        // 1. If the node is empty (no body), store the newBody here.
        // 2. If it's a leaf node (already has a body), subdivide the node into 4 quadrants,
        //    re-insert the old body, then insert the new body.
        // 3. If it's an internal node, update the center of mass and total mass,
        //    then find the appropriate quadrant for the new body and insert it there.
        std::cerr << "TreeNode::insertBody not implemented yet!" << std::endl;
    }

    // TODO: Implement this function
    void computeCenterOfMass() {
        // A recursive function that calculates the center of mass and total mass
        // for this node and all its children.
        // center_x = (sum of mass_i * x_i) / total_mass
        // center_y = (sum of mass_i * y_i) / total_mass
        std::cerr << "TreeNode::computeCenterOfMass not implemented yet!" << std::endl;
    }

    // TODO: Implement this function
    void calculateForce(Body* targetBody, double theta) {
        // Calculate the distance (r) between targetBody and this node's center of mass.
        // If (width / r < theta) OR the node is a leaf: 
        //    Use the center of mass to calculate the force on targetBody and add it.
        // Else:
        //    Recursively call calculateForce on all children.
        std::cerr << "TreeNode::calculateForce not implemented yet!" << std::endl;
    }
};

void generateRandomBodies(std::vector<Body>& bodies, int N, double region_size) {
    bodies.clear();
    bodies.reserve(N);
    for (int i = 0; i < N; ++i) {
        double x = static_cast<double>(rand()) / RAND_MAX * region_size - region_size/2.0;
        double y = static_cast<double>(rand()) / RAND_MAX * region_size - region_size/2.0;
        double mass = static_cast<double>(rand()) / RAND_MAX * 1e20; // Random mass
        bodies.emplace_back(x, y, 0, 0, mass);
    }
}

int main() {
    const int N = 1000;         // Number of bodies
    const double REGION_SIZE = 1000.0; // Size of the simulation space
    const int STEPS = 10;       // Number of simulation steps to run

    std::vector<Body> bodies;
    generateRandomBodies(bodies, N, REGION_SIZE);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < STEPS; ++step) {
        // 1. Build the Tree
        TreeNode root(0, 0, REGION_SIZE); // Root node centered at (0,0)
        for (Body& b : bodies) {
            root.insertBody(&b);
        }
        root.computeCenterOfMass();

        // 2. Calculate Forces for each body using the tree
        for (Body& b : bodies) {
            b.acc_x = 0;
            b.acc_y = 0;
            root.calculateForce(&b, THETA);
        }

        // 3. Update Positions and Velocities (Simple Euler Integration)
        for (Body& b : bodies) {
            b.vel_x += b.acc_x * DT;
            b.vel_y += b.acc_y * DT;
            b.pos_x += b.vel_x * DT;
            b.pos_y += b.vel_y * DT;
        }

        // (Optional) Output positions every X steps for visualization
        if (step % 10 == 0) {
            std::cout << "Step " << step << std::endl;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Sequential Simulation took " << elapsed.count() << " seconds." << std::endl;

    return 0;
}