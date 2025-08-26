#ifndef BODY_H
#define BODY_H

struct Body {
    double pos_x, pos_y; // Position (2D for simplicity)
    double vel_x, vel_y; // Velocity
    double acc_x, acc_y; // Acceleration
    double mass;

    // Constructor for easy initialization
    __host__ __device__
    Body(double x = 0, double y = 0, double vx = 0, double vy = 0, double m = 0)
        : pos_x(x), pos_y(y), vel_x(vx), vel_y(vy), acc_x(0), acc_y(0), mass(m) {}
};

#endif