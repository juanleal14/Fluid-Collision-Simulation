#ifndef FLUID_OPTIMIZATION_GRID_H
#define FLUID_OPTIMIZATION_GRID_H

#endif //FLUID_OPTIMIZATION_GRID_H

#include "libraries.h"
#include "block.cpp"


//std::vector<double> bmax = {0.065, 0.1, 0.065};
constexpr double bmax_coord_x = 0.065;
constexpr double bmax_coord_y = 0.1;
constexpr double bmax_coord_z = 0.065;
constexpr double bmin_coord_x = -0.065;
constexpr double bmin_coord_y = -0.08;
constexpr double bmin_coord_z = -0.065;
const double distance_minimum = pow(10, -10);

//class Bmax;
//class Bmax bmax();
//std::vector<double> bmin = {-0.065, -0.08, -0.065};
//class Bmin;
//class Bmin bmin();
//struct GridSize;
//struct Grid;
class vect3;
class GridSize;
class Grid;
class Initial_Values;

std::vector<Block> get_contiguous_blocks(int current_block, GridSize gsize);
Grid gridCreation(GridSize gridSize);
int find_block(Particle particle,GridSize gridSize);
Grid grid_initialization(Initial_Values &initialValues,std::vector<Particle> &particles);

void particle_collision_with_X_axis(Grid &grid);
void particle_collision_with_Y_axis(Grid &grid);
void particle_collision_with_Z_axis(Grid &grid);

void X_boundary_interaction(Grid &grid);
void Y_boundary_interaction(Grid &grid);
void Z_boundary_interaction(Grid &grid);

void particle_collision(Grid &grid);
void boundary_collision(Grid &grid);