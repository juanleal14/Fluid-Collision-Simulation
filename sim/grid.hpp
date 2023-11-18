#ifndef FLUID_OPTIMIZATION_GRID_H
#define FLUID_OPTIMIZATION_GRID_H
#endif

#include "block.cpp"


//std::vector<double> bmax = {0.065, 0.1, 0.065};


//class Bmax;
//class Bmax bmax();
//std::vector<double> bmin = {-0.065, -0.08, -0.065};
//class Bmin;
//class Bmin bmin();
//struct GridSize;
//struct Grid;
class GridSize;
class Grid;

Grid initialize_grid(std::ifstream &file,Initial_Values &initialValues,int &counter);
Grid initial_read(const std::string& file_address,Initial_Values &initialValues);

std::vector<Block> get_contiguous_blocks(int current_block, GridSize gsize);
Grid gridCreation(GridSize gridSize);
int find_block(Particle particle,GridSize gridSize);

void particle_collision_with_X_axis(Grid &grid);
void particle_collision_with_Y_axis(Grid &grid);
void particle_collision_with_Z_axis(Grid &grid);

void X_boundary_interaction(Grid &grid);
void Y_boundary_interaction(Grid &grid);
void Z_boundary_interaction(Grid &grid);

void particle_collision(Grid &grid);
void boundary_collision(Grid &grid);;


