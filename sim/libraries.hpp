#ifndef FLUID_LIBRARIES_H
#define FLUID_LIBRARIES_H


#endif //FLUID_LIBRARIES_H
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <numbers>
#include <span>
#include <algorithm>

constexpr double bmax_coord_x = 0.065;
constexpr double bmax_coord_y = 0.1;
constexpr double bmax_coord_z = 0.065;
constexpr double bmin_coord_x = -0.065;
constexpr double bmin_coord_y = -0.08;
constexpr double bmin_coord_z = -0.065;

static constexpr double const_r = 1.695;
static constexpr double global_density = 1000;
static constexpr double stiff_pressure = 3.0;
static constexpr double stiff_collision = 30000;   //stiffness collision
static constexpr double damping = 128.0;
static constexpr double viscosity = 0.4;
static constexpr double part_size = 0.0002; //Particle Size
static constexpr double time_step = 0.001;
const double distance_minimum = pow(10, -10);

const double gravity = 9.8;