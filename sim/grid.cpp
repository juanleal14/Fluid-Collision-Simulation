#include "grid.hpp"


class GridSize {
public:
    // Use member initializer list to initialize member variables
    GridSize(int x = 0, int y = 0, int z = 0, double sx = 0.0, double sy = 0.0, double sz = 0.0): {}

    // Getter methods for member variables
    [[nodiscard]] int getNumX() const { return num_blocks.x(); }
    [[nodiscard]] int getNumY() const { return num_blocks.y(); }
    [[nodiscard]] int getNumZ() const { return num_blocks.z(); }
    [[nodiscard]] double getSizeX() const { return sizes.x(); }
    [[nodiscard]] double getSizeY() const { return sizes.y(); }
    [[nodiscard]] double getSizeZ() const { return sizes.z(); }

    // Setter methods to modify member variables
    void setNumX(int x) { num_x = x; }
    void setNumY(int y) { num_y = y; }
    void setNumZ(int z) { num_z = z; }
    void setSizeX(double sx) { size_x = sx; }
    void setSizeY(double sy) { size_y = sy; }
    void setSizeZ(double sz) { size_z = sz; }

private:
    Vect3 num_blocks;
    Vect3 sizes;
};

class Grid {
    public:
        // Use member initializer list for member variable initialization
        Grid() : size(), block(size.getNumX(), std::vector<int>(size.getNumY(), 0)) {}

        // Public members for direct access
        GridSize size;
        std::vector<Block> blocks;
};

/*
class Bmin {
public:
    // Constructor to initialize constant features
    Bmin(double x, double y, double z) : bmin_x(x), bmin_y(y), bmin_z(z) {}

    // Getter methods to access individual features
    [[nodiscard]] double x() const { return bmin_x; }
    [[nodiscard]] double y() const { return bmin_y; }
    [[nodiscard]] double z() const { return bmin_z; }

private:
    double bmin_x = bmin_coord_x;
    double bmin_y = bmin_coord_y;
    double bmin_z = bmin_coord_z;
};*/
const Vect3 bmax(bmax_coord_x,bmax_coord_y,bmax_coord_z);
const Vect3 bmin(bmin_coord_x,bmin_coord_y,bmin_coord_z);



Grid gridCreation(GridSize gridSize){
    Grid grid;
    int total_iteration = (gridSize.getNumX())*(gridSize.getNumY())*(gridSize.getNumZ());
    for (int loop_x = 0; loop_x < total_iteration; loop_x++){
        std::vector <Particle> block;
        grid.blocks.push_back(block);
    }
    grid.size = gridSize;
    std::cout<<"Before : "<< grid.blocks.size()<<'\n';
    return grid;
}

std::vector<Block> get_contiguous_blocks(int current_block, Grid &grid){
    std::vector<Block> contiguous_blocks;
    const int block_z = current_block / (grid.size.getNumY() * grid.size.getNumX());
    const int block_y = (current_block - block_z * (grid.size.getNumY() * grid.size.getNumX())) / grid.size.getNumX();
    const int block_x = current_block - (block_z * grid.size.getNumY() * grid.size.getNumX()) - (block_y * grid.size.getNumX());

    for (int loop_x = -1; loop_x < 2; loop_x++){
        if ((block_x + loop_x >= 0) && (block_x + loop_x <= grid.size.getNumX() - 1)) {
            for (int loop_y = -1; loop_y < 2; loop_y++) {
                if ((block_y + loop_y >= 0) && (block_y + loop_y <= grid.size.getNumY() - 1)) {
                    for (int loop_z = -1; loop_z < 2; loop_z++) {
                        if ((block_z + loop_z >= 0) && (block_z + loop_z <= grid.size.getNumZ() - 1)){
                            const int computed_block = (block_x + loop_x) + (block_y + loop_y) * grid.size.getNumX() + (block_z + loop_z) * grid.size.getNumY() * grid.size.getNumX();
                            contiguous_blocks.push_back(grid.blocks[computed_block]);
                        }
                    }
                }
            }
        }
    }

    return contiguous_blocks;
}


int find_block(Particle particle,GridSize gridSize){
    int block_x = floor((particle.pos.x() - bmin.x())/gridSize.getSizeX());
    int block_y = floor((particle.pos.y() - bmin.y())/gridSize.getSizeY());
    int block_z = floor((particle.pos.z() - bmin.z())/gridSize.getSizeZ());
    if (block_x < 0){
        block_x = 0;
    } else if (block_x >= gridSize.getNumX()-1) {
        block_x = gridSize.getNumX()-1;
    }if (block_y < 0){
        block_y = 0;
    } else if (block_y >= gridSize.getNumY()-1) {
        block_y = gridSize.getNumY()-1;
    }if (block_z < 0){
        block_z = 0;
    } else if (block_z >= gridSize.getNumZ()-1) {
        block_z = gridSize.getNumZ()-1;
    }
    //cout << "This is the x block " << block_x << ", y block " << block_y << ", z block " << block_z;
    const int num_block = block_x + block_y*gridSize.getNumX() + block_z*gridSize.getNumY()*gridSize.getNumX();
    return num_block;
}

Grid grid_initialization(Initial_Values &initialValues,Block &particles){
    //Bmax bmax();
    const double boxx = bmax.getFeature1() - bmin.getFeature1();
    const double boxy = bmax.getFeature2() - bmin.getFeature2();
    const double boxz = bmax.getFeature3() - bmin.getFeature3();
    GridSize gridSize;
    gridSize.setNumX(floor(boxx/initialValues.getH()));
    gridSize.setNumY(floor(boxy/initialValues.getH()));
    gridSize.setNumZ(floor(boxz/initialValues.getH()));
    gridSize.setSizeX(boxx/gridSize.getNumX());
    gridSize.setSizeY(boxy/gridSize.getNumY());
    gridSize.setSizeZ(boxz/gridSize.getNumZ());
    const std::vector<std::vector <int>> blocks = gridCreation(gridSize);
    Grid grid;
    grid.size = gridSize;
    grid.blocks = blocks(grid.size.getNumX(), std::vector<int>(grid.size.getNumY(), 0);
    std::cout<<"After : "<<grid.blocks.size()<<'\n';
    for (int i = 0; i < particles.size(); i++){
        const int index = find_block(particles[i],grid.size);
        //cout<<i<<' ';
        grid.blocks[index].push_back(i);
    }
    return grid;
}

void particle_collision_with_Z_axis(Grid &grid) {

    double z_param = 0;
    double increment = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumX()*grid.size.getNumY(); loop_i++) {   //pared Z_0 //the number of particles in the x axis is num_y * num_x twice x min and xmax
        for (auto loop_j: grid.blocks[loop_i]) {
            z_param = loop_j.pos.z() + loop_j.hv.z() * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (z_param - bmin.z());                                //  part_size − (z − zmin)
            if (increment > distance_minimum) {                            //  az + (cs · ∆z − damping · vz)
            loop_j.acceleration.set_z(loop_j.acceleration.z() + (stiff_collision * increment - damping * loop_j.v.z()));
            }
        }
    }
    for (int loop_i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i++){
        for (auto loop_j : grid.blocks[loop_i]) {                                       //pared Z_max
            z_param = loop_j.pos.z() + loop_j.hv.z() * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (bmax.z() - z_param);                                //  part_size − (zmax − x)
            if (increment > distance_minimum) {                             //  az − (cs · ∆z + damping · vz)
                loop_j.acceleration.set_z(loop_j.acceleration.z() - (stiff_collision * increment + damping * loop_j.v.z()));
            }
        }
    }
}

void particle_collision_with_Y_axis(Grid &grid) {

    double y_param = 0;
    double increment = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumZ(); loop_i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (int loop_j = 0, loop_k=grid.size.getNumX()*(grid.size.getNumY()-1); loop_j < grid.size.getNumX(); loop_j++,loop_k++) {//pared Y_0
            for (auto loop_l : grid.blocks[ loop_j + loop_i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_x * grid.size.num_y;
                y_param = loop_l.pos.y()  + loop_l.hv.y() * time_step;   //  y = py + hvy · ∆t
                increment = part_size - (y_param - bmin.y());                        //  part_size − (y − ymin)grid.blocks[block_index]
                if (increment > distance_minimum) {                         //  ay + (stiff_collision · ∆y − damping · vy)
                    loop_l.acceleration.set_y(loop_l.acceleration.y() + (stiff_collision * increment - damping * loop_l.v.y()));
                }
            }                                                                //  pared Y_max
            for (auto loop_l : grid.blocks[loop_k + loop_i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_x * grid.size.num_y;
                y_param = loop_l.pos.y() + loop_l.hv.y() * time_step;    //  y = py + hvy · ∆t
                increment = part_size - (bmax.y() - y_param);                        //  part_size − (ymay − y) i
                if (increment > distance_minimum) {                       //  ay − (stiff_collision · ∆y + damping · vy)
                    loop_l.acceleration.set_y(loop_l.acceleration.y() - (stiff_collision * increment + damping * loop_l.v.y()));
                }
            }
        }
    }
}

void particle_collision_with_X_axis(std::vector<Particle> &particles, Grid &grid) {

    double x_param = 0 ;
    double increment = 0;
    for (int loop_i = 0,loop_j =grid.size.getNumX()-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i+=grid.size.getNumX(), loop_j+=grid.size.getNumX()){
        for (auto loop_l : grid.blocks[loop_i]) {                                       //pared X_0
            x_param = loop_l.pos.x() + loop_l.hv.x() * time_step;        //  x = px + hvx · ∆t
            increment = part_size - (x_param - bmin.x());                             //  part_size − (x − xmin)
            if (increment > distance_minimum) {                              //ax + (stiff_collision · ∆x − damping · vx)
                loop_l.acceleration.set_x(loop_l.acceleration.x() + (stiff_collision * increment - damping * loop_l.v.x()));
            }
        }
        for (auto loop_l : grid.blocks[loop_j]) {                                       //pared X_max
            x_param = loop_l.pos.x() + loop_l.hv.x() * time_step;        //x = px + hvx · ∆t
            increment = part_size - (bmax.x()- x_param);                             //part_size − (xmax− x)
            if (increment > distance_minimum) {                             //ax − (stiff_collision · ∆x + damping · vz)
                loop_l.acceleration.set_x(loop_l.acceleration.x() - (stiff_collision * increment + damping * loop_l.v.x()));
            }
        }
    }
}


void Z_boundary_interaction(std::vector<Particle> &particles, Grid &grid) {

    double distance_z = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumX()*grid.size.getNumY(); loop_i++) {   //pared x_0 //the number of particles in the x axis is num_y * num_z twice x min and xmax
        for (auto loop_j: grid.blocks[loop_i]) {
            distance_z = loop_j.pos.z() - bmin.z();
            if (distance_z < 0) {
                loop_j.pos.set_y(bmin.z() - distance_z);
                loop_j.v.set_z(-loop_j.v.z());
                loop_j.hv.set_z(-loop_j.hv.z());
            }
        }
    }
    for (int loop_i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i++){    //pared x_max
        for (auto loop_j : grid.blocks[loop_i]) {
            distance_z = bmax.z() - loop_j.pos.z();
            if (distance_z < 0) {
                loop_j.pos.set_z(bmax.z() + distance_z);
                loop_j.v.set_z(loop_j.v.z());
                loop_j.hv.set_z(-loop_j.hv.z());
            }
        }
    }
}

void Y_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_y=0;
    double test = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumZ(); loop_i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (int loop_j = 0, loop_k=grid.size.getNumX()*(grid.size.getNumY()-1); loop_j < grid.size.getNumX(); loop_j++,loop_k++) {//pared Y_0
            for (auto loop_l : grid.blocks[ loop_j + loop_i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_z * grid.size.num_y;
                distance_y = loop_l.pos.y() - bmin.y();
                if (distance_y < 0) {
                    loop_l.pos.set_y(bmin.y() + distance_y);
                    loop_l.v.set_y(loop_l.v.y());
                    loop_l.hv.set_y(-loop_l.hv.y());
                }
            }                                                                          //pared Y_max
            for (auto loop_n : grid.blocks[loop_k + loop_i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_z * grid.size.num_y;
                distance_y = bmax.y() - loop_n.pos.y();
                if (distance_y < 0) {
                    loop_n.pos.set_y(bmax.y() + distance_y);
                    loop_n.v.set_y(-loop_n.v.y());
                    loop_n.hv.set_y(-loop_n.hv.y());
                }
            }
        }
    }
}

void X_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_x=0;
    for (int loop_i = 0,loop_j =grid.size.getNumX()-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i+=grid.size.getNumX(), loop_j+=grid.size.getNumX()){
        for (auto loop_l : grid.blocks[loop_i]) {
            distance_x = loop_l.pos.x()-bmin.x();
            if (distance_x < 0){
                loop_l.pos.set_x(bmin.x() - distance_x);
                loop_l.v.set_x(loop_l.v.x());
                loop_l.hv.set_x(-loop_l.hv.x());
            }
        }//pared x_max
        for (auto loop_l : grid.blocks[loop_j]) {
            distance_x = bmax.x() - loop_l.pos.x();
            if (distance_x < 0) {
                loop_l.pos.set_x(bmax.x() + distance_x);
                loop_l.v.set_x(-loop_l.v.x());
                loop_l.hv.set_x(-loop_l.hv.x());
            }
        }
    }
}


void particle_collision(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &particles){
    particle_collision_with_Z_axis(particles,grid,particles);
    particle_collision_with_Y_axis(particles,grid,particles);
    particle_collision_with_X_axis(particles,grid,particles);
}
void boundary_collision(std::vector<Particle> &particles, Grid &grid){
    Z_boundary_interaction(particles,grid);
    Y_boundary_interaction(particles,grid);
    X_boundary_interaction(particles,grid);
}
