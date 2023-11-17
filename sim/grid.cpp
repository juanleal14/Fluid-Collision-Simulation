#include "grid.hpp"


class GridSize {
public:
    // Use member initializer list to initialize member variables
    GridSize(int x = 0, int y = 0, int z = 0, double sx = 0.0, double sy = 0.0, double sz = 0.0)
            : num_x(x), num_y(y), num_z(z), size_x(sx), size_y(sy), size_z(sz) {}

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
    int num_x;
    int num_y;
    int num_z;
    double size_x;
    double size_y;
    double size_z;
};


    class Grid {
    public:
        // Use member initializer list for member variable initialization
        Grid() : size(), blocks(size.getNumX(), std::vector<int>(size.getNumY(), 0)) {}

        // Public members for direct access
        GridSize size;
        std::vector<std::vector<int>> blocks;
    };

class vect3 {
public:
    // Constructor to initialize constant features
    vect3(const double x, const double y, const double z) : x(x), y(y), z(z) {}

    // Getter methods to access individual features
    [[nodiscard]] double getFeature1() const { return x; }
    [[nodiscard]] double getFeature2() const { return y; }
    [[nodiscard]] double getFeature3() const { return z; }

private:
    double x;
    double y;
    double z;
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
vect3 bmax(bmax_coord_x,bmax_coord_y,bmax_coord_z);
vect3 bmin(bmin_coord_x,bmin_coord_y,bmin_coord_z);
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

void particle_collision_with_Z_axis(std::vector<Particle> &particles,  Grid &grid, std::vector <Acceleration> &accelerations) {

    double z_param = 0;
    double increment = 0;
    for (int i = 0; i < grid.size.getNumX()*grid.size.getNumY(); i++) {   //pared Z_0 //the number of particles in the x axis is num_y * num_x twice x min and xmax
        for (int j: grid.blocks[i]) {
            z_param = particles[j].pz + particles[j].hvz * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (z_param - bmin.getFeature3());                                //  part_size − (z − zmin)
            if (increment > pow(10, -10)) {                            //  az + (cs · ∆z − damping · vz)
            accelerations[j].az = accelerations[j].az + (stiff_collision * increment - damping * particles[j].vz);
            }
        }
    }
    for (int i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); i++){
        for (const int j : grid.blocks[i]) {                                       //pared Z_max
            z_param = particles[j].pz + particles[j].hvz * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (bmax.getFeature3() - z_param);                                //  part_size − (zmax − x)
            if (increment > pow(10, -10)) {                             //  az − (cs · ∆z + damping · vz)
                accelerations[j].az = accelerations[j].az - (stiff_collision * increment + damping * particles[j].vz);
            }
        }
    }
}

void particle_collision_with_Y_axis(std::vector<Particle> &particles , Grid &grid, std::vector <Acceleration> &accelerations) {

    double y_param = 0;
    double increment = 0;
    for (int i = 0; i < grid.size.getNumZ(); i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (long j = 0, k=grid.size.getNumX()*(grid.size.getNumY()-1); j < grid.size.getNumX(); j++,k++) {//pared Y_0
            for (int l : grid.blocks[ j + i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_x * grid.size.num_y;
                y_param = particles[l].py  + particles[l].hvy * time_step;   //  y = py + hvy · ∆t
                increment = part_size - (y_param - bmin.getFeature2());                        //  part_size − (y − ymin)grid.blocks[block_index]
                if (increment > pow(10, -10)) {                         //  ay + (stiff_collision · ∆y − damping · vy)
                    accelerations[l].ay = accelerations[l].ay + (stiff_collision * increment - damping * particles[l].vy);
                }
            }                                                                //  pared Y_max
            for (int l : grid.blocks[k + i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_x * grid.size.num_y;
                y_param = particles[l].py + particles[l].hvy * time_step;    //  y = py + hvy · ∆t
                increment = part_size - (bmax.getFeature2() - y_param);                        //  part_size − (ymay − y) i
                if (increment > pow(10, -10)) {                       //  ay − (stiff_collision · ∆y + damping · vy)
                    accelerations[l].ay = accelerations[l].ay - (stiff_collision * increment + damping * particles[l].vy);
                }
            }
        }
    }
}

void particle_collision_with_X_axis(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations) {

    double x_param = 0 ;
    double increment = 0;
    for (int i = 0,j =grid.size.getNumX()-1; i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); i+=grid.size.getNumX(), j+=grid.size.getNumX()){
        for (int l : grid.blocks[i]) {                                       //pared X_0
            x_param = particles[l].px + particles[l].hvx * time_step;        //  x = px + hvx · ∆t
            increment = part_size - (x_param - bmin.getFeature1());                             //  part_size − (x − xmin)
            if (increment > pow(10, -10)) {                              //ax + (stiff_collision · ∆x − damping · vx)
                accelerations[l].ax = accelerations[l].ax + (stiff_collision * increment - damping * particles[l].vx);
            }
        }
        for (int l : grid.blocks[j]) {                                       //pared X_max
            x_param = particles[l].px + particles[l].hvx * time_step;        //x = px + hvx · ∆t
            increment = part_size - (bmax.getFeature1()- x_param);                             //part_size − (xmax− x)
            if (increment > pow(10, -10)) {                             //ax − (stiff_collision · ∆x + damping · vz)
                accelerations[l].ax = accelerations[l].ax - (stiff_collision * increment + damping * particles[l].vx);
            }
        }
    }
}


void Z_boundary_interaction(std::vector<Particle> &particles, Grid &grid) {

    double distance_z = 0;
    for (int i = 0; i < grid.size.getNumX()*grid.size.getNumY(); i++) {   //pared x_0 //the number of particles in the x axis is num_y * num_z twice x min and xmax
        for (int j: grid.blocks[i]) {
            distance_z = particles[j].pz - bmin.getFeature3();
            if (distance_z < 0) {
                particles[j].py = bmin.getFeature3() - distance_z;
                particles[j].vz = -particles[j].vz;
                particles[j].hvz = -particles[j].hvz;
            }
        }
    }
    for (int i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); i++){    //pared x_max
        for (int j : grid.blocks[i]) {
            distance_z = bmax.getFeature3() - particles[j].pz;
            if (distance_z < 0) {
                particles[j].pz = bmax.getFeature3() + distance_z;
                particles[j].vz = -particles[j].vz;
                particles[j].hvz = -particles[j].hvz;
            }
        }
    }
}

void Y_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_y=0;
    for (int i = 0; i < grid.size.getNumZ(); i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (int j = 0, k=grid.size.getNumX()*(grid.size.getNumY()-1); j < grid.size.getNumX(); j++,k++) {//pared Y_0
            for (int l : grid.blocks[ j + i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_z * grid.size.num_y;
                distance_y = particles[l].py - bmin.getFeature2();
                if (distance_y < 0) {
                    particles[l].py = bmin.getFeature2() + distance_y;
                    particles[l].vy = -particles[l].vy;
                    particles[l].hvy = -particles[l].hvy;
                }
            }                                                                          //pared Y_max
            for (int n : grid.blocks[k + i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_z * grid.size.num_y;
                distance_y = bmax.getFeature2() - particles[n].py;
                if (distance_y < 0) {
                    particles[n].py = bmax.getFeature2() + distance_y;
                    particles[n].vy = -particles[n].vy;
                    particles[n].hvy = -particles[n].hvy;
                }
            }
        }
    }
}

void X_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_x=0;
    for (int i = 0,j =grid.size.getNumX()-1; i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); i+=grid.size.getNumX(), j+=grid.size.getNumX()){
        for (int l : grid.blocks[i]) {
            distance_x = particles[j].px-bmin.getFeature1();
            if (distance_x < 0){
                particles[j].px = bmin.getFeature1() - distance_x;
                particles[j].vx = -particles[j].vx;
                particles[j].hvx = -particles[j].hvx;
            }
        }//pared x_max
        for (int l : grid.blocks[j]) {
            distance_x = bmax.getFeature1() - particles[j].px;
            if (distance_x < 0) {
                particles[j].px = bmax.getFeature1() + distance_x;
                particles[j].vx = -particles[j].vx;
                particles[j].hvx = -particles[j].hvx;
            }
        }
    }
}


void particle_collision(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations){
    particle_collision_with_Z_axis(particles,grid,accelerations);
    particle_collision_with_Y_axis(particles,grid,accelerations);
    particle_collision_with_X_axis(particles,grid,accelerations);
}
void boundary_collision(std::vector<Particle> &particles, Grid &grid){
    Z_boundary_interaction(particles,grid);
    Y_boundary_interaction(particles,grid);
    X_boundary_interaction(particles,grid);
}
