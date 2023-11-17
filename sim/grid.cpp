#include "grid.hpp"


class GridSize {
public:
    // Use member initializer list to initialize member variables
    GridSize(int x = 0, int y = 0, int z = 0, double sx = 0.0, double sy = 0.0, double sz = 0.0): num_x(x), num_y(y), num_z(z), size_x(sx), size_y(sy), size_z(sz) {}

    // Getter methods for member variables
    [[nodiscard]] int getNumX() const { return num_x; }
    [[nodiscard]] int getNumY() const { return num_y; }
    [[nodiscard]] int getNumZ() const { return num_z; }
    [[nodiscard]] double getSizeX() const { return size_x; }
    [[nodiscard]] double getSizeY() const { return size_y; }
    [[nodiscard]] double getSizeZ() const { return size_z; }

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
        Grid() : size(), blocks(size.getNumX(), std::vector<int>(size.getNumY(), 0)) {}

        // Public members for direct access
        GridSize size;
        std::vector<std::vector<Particle>> blocks;
    };

/*
class Bmin {
public:
    // Constructor to initialize constant features
    Bmin(double x, double y, double z) : bmin_x(x), bmin_y(y), bmin_z(z) {}

    // Getter methods to access individual features
    [[nodiscard]] double getFeature1() const { return bmin_x; }
    [[nodiscard]] double getFeature2() const { return bmin_y; }
    [[nodiscard]] double getFeature3() const { return bmin_z; }

private:
    double bmin_x = bmin_coord_x;
    double bmin_y = bmin_coord_y;
    double bmin_z = bmin_coord_z;
};*/
const Vect3 bmax(bmax_coord_x,bmax_coord_y,bmax_coord_z);
const Vect3 bmin(bmin_coord_x,bmin_coord_y,bmin_coord_z);
std::vector <std::vector<int>> gridCreation(GridSize gridSize){
    std::vector<std::vector<int>> blocks;
    for (int x = 0; x < (gridSize.getNumX())*(gridSize.getNumY())*(gridSize.getNumZ()); x++){
        const std::vector <int> new_vector;
        blocks.push_back(new_vector);
    }
    std::cout<<"Before : "<< blocks.size()<<'\n';
    return blocks;
}

std::vector<int> get_contiguous_blocks(int current_block, GridSize gsize){
    std::vector<int> contiguous_blocks;
    const int block_z = current_block / (gsize.getNumY() * gsize.getNumX());
    const int block_y = (current_block - block_z * (gsize.getNumY() * gsize.getNumX())) / gsize.getNumX();
    const int block_x = current_block - (block_z * gsize.getNumY() * gsize.getNumX()) - (block_y * gsize.getNumX());

    for (int x = -1; x < 2; x++){
        if ((block_x + x >= 0) && (block_x + x <= gsize.getNumX() - 1)) {
            for (int y = -1; y < 2; y++) {
                if ((block_y + y >= 0) && (block_y + y <= gsize.getNumY() - 1)) {
                    for (int z = -1; z < 2; z++) {
                        if ((block_z + z >= 0) && (block_z + z <= gsize.getNumZ() - 1)){
                            const int computed_block = (block_x + x) + (block_y + y) * gsize.getNumX() + (block_z + z) * gsize.getNumY() * gsize.getNumX();
                            contiguous_blocks.push_back(computed_block);
                        }
                    }
                }
            }
        }
    }

    return contiguous_blocks;
}


int find_block(Particle particle,GridSize gridSize){
    int block_x = floor((particle.px - bmin.getFeature1())/gridSize.getSizeX());
    int block_y = floor((particle.py - bmin.getFeature2())/gridSize.getSizeY());
    int block_z = floor((particle.pz - bmin.getFeature3())/gridSize.getSizeZ());
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

Grid grid_initialization(Initial_Values &initialValues,std::vector<Particle> &particles){
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
    grid.blocks = blocks;
    std::cout<<"After : "<<grid.blocks.size()<<'\n';
    for (int i = 0; i < particles.size(); i++){
        const int index = find_block(particles[i],grid.size);
        //cout<<i<<' ';
        grid.blocks[index].push_back(i);
    }
    return grid;
}

void particle_collision_with_Z_axis( Grid &grid) {

    double z_param = 0;
    double increment = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumX()*grid.size.getNumY(); loop_i++) {   //pared Z_0 //the number of particles in the x axis is num_y * num_x twice x min and xmax
        for (const int loop_j: grid.blocks[loop_i]) {
            z_param = particles[loop_j].pos().z() + particles[loop_j].hv().z() * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (z_param - bmin.z());                                //  part_size − (z − zmin)
            if (increment > distance_minimum) {                            //  az + (cs · ∆z − damping · vz)
            particles[loop_j].acc().z() = particles[loop_j].acc().z() + (stiff_collision * increment - damping * particles[loop_j].vel().z());
            }
        }
    }
    for (int loop_i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i++){
        for (const int loop_j : grid.blocks[loop_i]) {                                       //pared Z_max
            z_param = particles[loop_j].pos().z() + particles[loop_j].hv().z() * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (bmax.z() - z_param);                                //  part_size − (zmax − x)
            if (increment > distance_minimum) {                             //  az − (cs · ∆z + damping · vz)
                particles[loop_j].acc().z() = particles[loop_j].acc().z() - (stiff_collision * increment + damping * particles[loop_j].vel().z());
            }
        }
    }
}

void particle_collision_with_Y_axis(std::vector<Particle> &particles , Grid &grid, std::vector <Acceleration> &particles) {

    double y_param = 0;
    double increment = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumZ(); loop_i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (long loop_j = 0, loop_k=grid.size.getNumX()*(grid.size.getNumY()-1); loop_j < grid.size.getNumX(); loop_j++,loop_k++) {//pared Y_0
            for (int loop_l : grid.blocks[ loop_j + loop_i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_x * grid.size.num_y;
                y_param = particles[loop_l].pos().y()  + particles[loop_l].hv().y() * time_step;   //  y = py + hvy · ∆t
                increment = part_size - (y_param - bmin.y());                        //  part_size − (y − ymin)grid.blocks[block_index]
                if (increment > distance_minimum) {                         //  ay + (stiff_collision · ∆y − damping · vy)
                    particles[loop_l].acc().y() = particles[loop_l].acc().y() + (stiff_collision * increment - damping * particles[loop_l].vel().y());
                }
            }                                                                //  pared Y_max
            for (int loop_l : grid.blocks[loop_k + loop_i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_x * grid.size.num_y;
                y_param = particles[loop_l].pos().y() + particles[loop_l].hv().y() * time_step;    //  y = py + hvy · ∆t
                increment = part_size - (bmax.y() - y_param);                        //  part_size − (ymay − y) i
                if (increment > distance_minimum) {                       //  ay − (stiff_collision · ∆y + damping · vy)
                    particles[loop_l].acc().y() = particles[loop_l].acc().y() - (stiff_collision * increment + damping * particles[loop_l].vel().y());
                }
            }
        }
    }
}

void particle_collision_with_X_axis(std::vector<Particle> &particles, Grid &grid) {

    double x_param = 0 ;
    double increment = 0;
    for (int loop_i = 0,loop_j =grid.size.getNumX()-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i+=grid.size.getNumX(), loop_j+=grid.size.getNumX()){
        for (const int loop_l : grid.blocks[loop_i]) {                                       //pared X_0
            x_param = particles[loop_l].pos().x() + particles[loop_l].hv().x() * time_step;        //  x = px + hvx · ∆t
            increment = part_size - (x_param - bmin.x());                             //  part_size − (x − xmin)
            if (increment > distance_minimum) {                              //ax + (stiff_collision · ∆x − damping · vx)
                particles[loop_l].acc().x() = particles[loop_l].acc().x() + (stiff_collision * increment - damping * particles[loop_l].vel().x());
            }
        }
        for (const int loop_l : grid.blocks[loop_j]) {                                       //pared X_max
            x_param = particles[loop_l].pos().x() + particles[loop_l].hv().x() * time_step;        //x = px + hvx · ∆t
            increment = part_size - (bmax.x()- x_param);                             //part_size − (xmax− x)
            if (increment > distance_minimum) {                             //ax − (stiff_collision · ∆x + damping · vz)
                particles[loop_l].acc().x() = particles[loop_l].acc().x() - (stiff_collision * increment + damping * particles[loop_l].vel().x());
            }
        }
    }
}


void Z_boundary_interaction(std::vector<Particle> &particles, Grid &grid) {

    double distance_z = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumX()*grid.size.getNumY(); loop_i++) {   //pared x_0 //the number of particles in the x axis is num_y * num_z twice x min and xmax
        for (int loop_j: grid.blocks[loop_i]) {
            distance_z = particles[loop_j].pos().z() - bmin.z();
            if (distance_z < 0) {
                particles[loop_j].pos().y() = bmin.z() - distance_z;
                particles[loop_j].vel().z() = -particles[loop_j].vel().z();
                particles[loop_j].hv().z() = -particles[loop_j].hv().z();
            }
        }
    }
    for (int loop_i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i++){    //pared x_max
        for (int loop_j : grid.blocks[loop_i]) {
            distance_z = bmax.z() - particles[loop_j].pos().z();
            if (distance_z < 0) {
                particles[loop_j].pos().z() = bmax.z() + distance_z;
                particles[loop_j].vel().z() = -particles[loop_j].vel().z();
                particles[loop_j].hv().z() = -particles[loop_j].hv().z();
            }
        }
    }
}

void Y_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_y=0;
    for (int loop_i = 0; loop_i < grid.size.getNumZ(); loop_i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (int loop_j = 0, loop_k=grid.size.getNumX()*(grid.size.getNumY()-1); loop_j < grid.size.getNumX(); loop_j++,loop_k++) {//pared Y_0
            for (int loop_l : grid.blocks[ loop_j + loop_i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_z * grid.size.num_y;
                distance_y = particles[loop_l].pos().y() - bmin.y();
                if (distance_y < 0) {
                    particles[loop_l].pos().y() = bmin.y() + distance_y;
                    particles[loop_l].vel().y() = -particles[loop_l].vel().y();
                    particles[loop_l].hv().y() = -particles[loop_l].hv().y();
                }
            }                                                                          //pared Y_max
            for (int loop_n : grid.blocks[loop_k + loop_i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_z * grid.size.num_y;
                distance_y = bmax.y() - particles[loop_n].pos().y();
                if (distance_y < 0) {
                    particles[loop_n].pos().y() = bmax.y() + distance_y;
                    particles[loop_n].vel().y() = -particles[loop_n].vel().y();
                    particles[loop_n].hv().y() = -particles[loop_n].hv().y();
                }
            }
        }
    }
}

void X_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_x=0;
    for (int loop_i = 0,loop_j =grid.size.getNumX()-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i+=grid.size.getNumX(), loop_j+=grid.size.getNumX()){
        for (int loop_l : grid.blocks[loop_i]) {
            distance_x = particles[loop_j].pos().x()-bmin.x();
            if (distance_x < 0){
                particles[loop_j].pos().x() = bmin.x() - distance_x;
                particles[loop_j].vel().x() = -particles[loop_j].vel().x();
                particles[loop_j].hv().x() = -particles[loop_j].hv().x();
            }
        }//pared x_max
        for (int loop_l : grid.blocks[loop_j]) {
            distance_x = bmax.x() - particles[loop_j].pos().x();
            if (distance_x < 0) {
                particles[loop_j].pos().x() = bmax.x() + distance_x;
                particles[loop_j].vel().x() = -particles[loop_j].vel().x();
                particles[loop_j].hv().x() = -particles[loop_j].hv().x();
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
