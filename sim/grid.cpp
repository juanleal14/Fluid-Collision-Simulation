#include "grid.hpp"


class GridSize {
public:
    // Use member initializer list to initialize member variables
    GridSize(): num_blocks(Vect3<int>(0,0,0)), sizes(Vect3<double>(0,0,0)) {}
    GridSize(Vect3<int> n, Vect3<double> s): num_blocks(n), sizes(s){}

    // Getter methods for member variables
    [[nodiscard]] int getNumX() const { return num_blocks.x(); }
    [[nodiscard]] int getNumY() const { return num_blocks.y(); }
    [[nodiscard]] int getNumZ() const { return num_blocks.z(); }
    [[nodiscard]] double getSizeX() const { return sizes.x(); }
    [[nodiscard]] double getSizeY() const { return sizes.y(); }
    [[nodiscard]] double getSizeZ() const { return sizes.z(); }

    // Setter methods to modify member variables
    void setNumX(int x) { num_blocks.set_x(x); }
    void setNumY(int y) { num_blocks.set_y(y); }
    void setNumZ(int z) { num_blocks.set_z(z); }
    void setSizeX(double sx) { sizes.set_x(sx); }
    void setSizeY(double sy) { sizes.set_y(sy); }
    void setSizeZ(double sz) { sizes.set_z(sz); }


private:
    Vect3 <int> num_blocks; //REVISAR, NO ESTOY SEGURO(JUAN)
    Vect3<double> sizes;
};

class Grid {
    public:
        // Use member initializer list for member variable initialization
        Grid()= default;
        Grid(GridSize gsize): size(gsize){
            int const total_iteration = (gsize.getNumX())*(gsize.getNumY())*(gsize.getNumZ());
            size = gsize;
            blocks = std::vector<Block>(total_iteration);
        }
        Block& operator[](int i){return blocks[i];}

        void add_particle(Particle p){
            int const block_num = find_block(p,size);
            blocks[block_num].push_back(p);
        }

        // Public members for direct access
        GridSize size;
        std::vector<Block> blocks;


};


///gridCreation implementada dentro de clase Grid mediante default constructor
/*
Grid gridCreation(GridSize gridSize){
    Grid grid;
    int const total_iteration = (gridSize.getNumX())*(gridSize.getNumY())*(gridSize.getNumZ());
    for (int loop_x = 0; loop_x < total_iteration; loop_x++){
        Block const block;
        grid.blocks.push_back(block);
    }
    grid.size = gridSize;
    std::cout<<"Before : "<< grid.blocks.size()<<'\n';
    return grid;
}*/
Grid initialize_grid(std::ifstream &file,Initial_Values &initialValues,int &counter){ /// Crear Grid con el size correcto y añadir todas las partículas del archivo

    const double boxx = bmax_coord_x - bmin_coord_x;
    const double boxy = bmax_coord_y - bmin_coord_y;
    const double boxz = bmax_coord_z - bmin_coord_z;
    GridSize gridSize;
    gridSize.setNumX(floor(boxx/initialValues.getH()));
    gridSize.setNumY(floor(boxy/initialValues.getH()));
    gridSize.setNumZ(floor(boxz/initialValues.getH()));
    gridSize.setSizeX(boxx/gridSize.getNumX());
    gridSize.setSizeY(boxy/gridSize.getNumY());
    gridSize.setSizeZ(boxz/gridSize.getNumZ());
    Grid grid(gridSize) ;

    while(file.peek()!=EOF){
        grid.add_particle(read_particle(file,counter));
        counter ++;
    }
    std::cout<<"Total particles read from file = "<<counter<<'\n';
    return grid;
}


Grid initial_read(const std::string& file_address,Initial_Values &initialValues){
    //Read file
    std::ifstream file(file_address, std::ios::binary);
    int counter = 0;
    if (!file.is_open()) { //Check error opening
        std::cerr<<"Error: Cannot open " << file_address <<" for reading";
        exit (-3);
    }
    initialValues = read_general_info(file);//call to a function to read parameters
    std::cout<<"Lee bien ppm np";
    Grid grid = initialize_grid(file,initialValues,counter);
    if(counter == 0){
        std::cout<< "Error : Invalid number of particles: " << counter <<".";
    }
    else if(counter != initialValues.getNp()){
        std::cout<<"Error : Number of particles mismatch. Header " << initialValues.getNp() << " Found " << counter <<".";
    }
    file.close();
    return grid;
}

std::vector<int> get_contiguous_blocks(int current_block, Grid &grid){
    std::vector<int> contiguous_blocks;
    const int block_z = current_block / (grid.size.getNumY() * grid.size.getNumX());
    const int block_y = (current_block - block_z * (grid.size.getNumY() * grid.size.getNumX())) / grid.size.getNumX();
    const int block_x = current_block - (block_z * grid.size.getNumY() * grid.size.getNumX()) - (block_y * grid.size.getNumX());

    for (int loop_z = 0; loop_z < 2; loop_z++){
        if ((block_z + loop_z >= 0) && (block_z + loop_z <= grid.size.getNumZ() - 1)) {
            for (int loop_y = 0; loop_y < 2; loop_y++) {
                if ((block_y + loop_y >= 0) && (block_y + loop_y <= grid.size.getNumY() - 1)) {
                    for (int loop_x = 0; loop_x < 2; loop_x++) {
                        if ((block_x + loop_x >= 0) && (block_x + loop_x <= grid.size.getNumX() - 1)){
                            const int computed_block = (block_x + loop_x) + (block_y + loop_y) * grid.size.getNumX() + (block_z + loop_z) * grid.size.getNumY() * grid.size.getNumX();
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
    int block_x = floor((particle.pos.x() - bmin_coord_x)/gridSize.getSizeX());
    int block_y = floor((particle.pos.y() - bmin_coord_y)/gridSize.getSizeY());
    int block_z = floor((particle.pos.z() - bmin_coord_z)/gridSize.getSizeZ());
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


void particle_collision_with_Z_axis(Grid &grid) {

    double z_param = 0;
    double increment = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumX()*grid.size.getNumY(); loop_i++) {   //pared Z_0 //the number of particles in the x axis is num_y * num_x twice x min and xmax
        for (auto loop_j: grid.blocks[loop_i]) {
            z_param = loop_j.pos.z() + loop_j.hv.z() * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (z_param - bmin_coord_z);                                //  part_size − (z − zmin)
            if (increment > distance_minimum) {                            //  az + (cs · ∆z − damping · vz)
            //loop_j.acceleration.set_z(loop_j.acceleration.z() + (stiff_collision * increment - damping * loop_j.v.z()));
                loop_j.acceleration[2] = loop_j.acceleration[2] + (stiff_collision * increment - damping * loop_j.v[2]);
            }
        }
    }
    for (int loop_i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() - grid.size.getNumX()*grid.size.getNumY())-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i++){
        for (auto loop_j : grid.blocks[loop_i]) {                                       //pared Z_max
            z_param = loop_j.pos.z() + loop_j.hv.z() * time_step;        //  z = pz + hvz · ∆t
            increment = part_size - (bmax_coord_z - z_param);                                //  part_size − (zmax − x)
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
                increment = part_size - (y_param - bmin_coord_y);                        //  part_size − (y − ymin)grid.blocks[block_index]
                if (increment > distance_minimum) {                         //  ay + (stiff_collision · ∆y − damping · vy)
                    loop_l.acceleration.set_y(loop_l.acceleration.y() + (stiff_collision * increment - damping * loop_l.v.y()));
                }
            }                                                                //  pared Y_max
            for (auto loop_l : grid.blocks[loop_k + loop_i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_x * grid.size.num_y;
                y_param = loop_l.pos.y() + loop_l.hv.y() * time_step;    //  y = py + hvy · ∆t
                increment = part_size - (bmax_coord_y - y_param);                        //  part_size − (ymay − y) i
                if (increment > distance_minimum) {                       //  ay − (stiff_collision · ∆y + damping · vy)
                    loop_l.acceleration.set_y(loop_l.acceleration.y() - (stiff_collision * increment + damping * loop_l.v.y()));
                }
            }
        }
    }
}

void particle_collision_with_X_axis(Grid &grid) {

    double x_param = 0 ;
    double increment = 0;
    for (int loop_i = 0,loop_j =grid.size.getNumX()-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i+=grid.size.getNumX(), loop_j+=grid.size.getNumX()){
        for (auto loop_l : grid.blocks[loop_i]) {                                       //pared X_0
            x_param = loop_l.pos.x() + loop_l.hv.x() * time_step;        //  x = px + hvx · ∆t
            increment = part_size - (x_param - bmin_coord_x);                             //  part_size − (x − xmin)
            if (increment > distance_minimum) {                              //ax + (stiff_collision · ∆x − damping · vx)
                loop_l.acceleration.set_x(loop_l.acceleration.x() + (stiff_collision * increment - damping * loop_l.v.x()));
            }
        }
        for (auto loop_l : grid.blocks[loop_j]) {                                       //pared X_max
            x_param = loop_l.pos.x() + loop_l.hv.x() * time_step;        //x = px + hvx · ∆t
            increment = part_size - (bmax_coord_x- x_param);                             //part_size − (xmax− x)
            if (increment > distance_minimum) {                             //ax − (stiff_collision · ∆x + damping · vz)
                loop_l.acceleration.set_x(loop_l.acceleration.x() - (stiff_collision * increment + damping * loop_l.v.x()));
            }
        }
    }
}


void Z_boundary_interaction(Grid &grid) {

    double distance_z = 0;
    for (int loop_i = 0; loop_i < grid.size.getNumX()*grid.size.getNumY(); loop_i++) {   //pared x_0 //the number of particles in the x axis is num_y * num_z twice x min and xmax
        for (auto loop_j: grid.blocks[loop_i]) {
            distance_z = loop_j.pos.z() - bmin_coord_z;
            if (distance_z < 0) {
                loop_j.pos.set_y(bmin_coord_z - distance_z);
                loop_j.v.set_z(-loop_j.v.z());
                loop_j.hv.set_z(-loop_j.hv.z());
            }
        }
    }
    for (int loop_i = (grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX() -
            grid.size.getNumX()*grid.size.getNumY())-1; loop_i < grid.size.getNumZ()*
            grid.size.getNumY()*grid.size.getNumX(); loop_i++){    //pared x_max
        for (auto loop_j : grid.blocks[loop_i]) {
            distance_z = bmax_coord_z - loop_j.pos.z();
            if (distance_z < 0) {
                loop_j.pos.set_z(bmax_coord_z + distance_z);
                loop_j.v.set_z(loop_j.v.z());
                loop_j.hv.set_z(-loop_j.hv.z());
            }
        }
    }
}

void Y_boundary_interaction(Grid &grid) {

    double distance_y=0;
    for (int loop_i = 0; loop_i < grid.size.getNumZ(); loop_i++){     //the number of particles in the y axis is num_x * num_z, twice y min and xmax
        for (int loop_j = 0, loop_k=grid.size.getNumX()*(grid.size.getNumY()-1); loop_j < grid.size.getNumX(); loop_j++,loop_k++) {//pared Y_0
            for (auto loop_l : grid.blocks[loop_j + loop_i * grid.size.getNumX() * grid.size.getNumY()]) {//  block_index = j + i * grid.size.num_z * grid.size.num_y;
                distance_y = loop_l.pos.y() - bmin_coord_y;
                if (distance_y < 0) {
                    loop_l.pos.set_y(bmin_coord_y + distance_y);
                    loop_l.v.set_y(loop_l.v.y());
                    loop_l.hv.set_y(-loop_l.hv.y());
                }
            }                                                                          //pared Y_max
            for (auto loop_n : grid.blocks[loop_k + loop_i * grid.size.getNumX() * grid.size.getNumY()]) { //  block_index = k + i * grid.size.num_z * grid.size.num_y;
                distance_y = bmax_coord_y - loop_n.pos.y();
                if (distance_y < 0) {
                    loop_n.pos.set_y(bmax_coord_y + distance_y);
                    loop_n.v.set_y(-loop_n.v.y());
                    loop_n.hv.set_y(-loop_n.hv.y());
                }
            }
        }
    }
}

void X_boundary_interaction(Grid &grid) {

    double distance_x=0;
    for (int loop_i = 0,loop_j =grid.size.getNumX()-1; loop_i < grid.size.getNumZ()*grid.size.getNumY()*grid.size.getNumX(); loop_i+=grid.size.getNumX(), loop_j+=grid.size.getNumX()){
        for (auto loop_l : grid.blocks[loop_i]) {
            distance_x = loop_l.pos.x()-bmin_coord_x;
            if (distance_x < 0){
                loop_l.pos.set_x(bmin_coord_x - distance_x);
                loop_l.v.set_x(loop_l.v.x());
                loop_l.hv.set_x(-loop_l.hv.x());
            }
        }//pared x_max
        for (auto loop_l : grid.blocks[loop_j]) {
            distance_x = bmax_coord_x - loop_l.pos.x();
            if (distance_x < 0) {
                loop_l.pos.set_x(bmax_coord_x + distance_x);
                loop_l.v.set_x(-loop_l.v.x());
                loop_l.hv.set_x(-loop_l.hv.x());
            }
        }
    }
}

void particle_collision(Grid &grid){
    particle_collision_with_Z_axis(grid);
    particle_collision_with_Y_axis(grid);
    particle_collision_with_X_axis(grid);
}
void boundary_collision(Grid &grid){
    Z_boundary_interaction(grid);
    Y_boundary_interaction(grid);
    X_boundary_interaction(grid);
}
//
Vect3<int> belongs_to_boundary_part(Particle particle, GridSize &gridSize) {
    Vect3<int> belongings(0,0,0) ;
    const int block_x = floor((particle.pos.x() - bmin_coord_x)/gridSize.getNumX());
    if (block_x == 0) {
        belongings[0]=-1;}
    else if (block_x == gridSize.getNumX()-1){
        belongings[0]=1;}
    const int block_y = floor((particle.pos.y() - bmin_coord_y)/gridSize.getSizeY());
    if (block_y == 0) {
        belongings[1]=-1;}
    else if (block_y == gridSize.getNumY()-1){
        belongings[1]=1;}
    const int block_z = floor((particle.pos.z() - bmin_coord_z)/gridSize.getSizeZ());
    if (block_z == 0) {
        belongings[2]=-1;}
    else if (block_z == gridSize.getSizeZ()-1){
        belongings[2]=1;}
    return belongings;
}

Vect3<int> belongs_to_boundary_block(int block_index, GridSize &gridSize) {
    Vect3<int> belongings(0,0,0) ;
    if (block_index % gridSize.getNumX()== 0  ) {
        belongings[0]=-1;}
    else if (block_index % gridSize.getNumX() == gridSize.getNumX()-1){
        belongings[0]=1;}
    if (block_index % gridSize.getNumY()== 0) {
        belongings[1]=-1;}
    else if (block_index % gridSize.getNumY() == gridSize.getNumY()-1){
        belongings[1]=1;}
    if (block_index <gridSize.getSizeY()*gridSize.getSizeX()) {
        belongings[2]=-1;}
    else if (block_index >(gridSize.getNumX()*gridSize.getNumY()*(gridSize.getNumZ()-1))){
        belongings[2]=1;
    }

    return belongings;
}


void new_boundary_collision(Vect3<int> belongings, Particle &particle) {
    double cord_param = 0;
    int wall = 0;
    double increment = 0;
    Vect3<double> bmin(bmin_coord_x, bmin_coord_y, bmin_coord_z);
    Vect3<double> bmax(bmax_coord_x, bmax_coord_y, bmax_coord_z);
    for (int loop_i = 0; loop_i < 3; loop_i++) {
        wall = belongings[loop_i];
        cord_param = particle.pos[loop_i] + particle.hv[loop_i] * time_step;
        if (wall == -1) { //WALL_MIN
            increment = part_size - (cord_param - bmin[loop_i]);
            if (increment > distance_minimum) {
                particle.acceleration[loop_i] = particle.acceleration[loop_i] + (stiff_collision * increment - damping * particle.v[loop_i]);
            }
        }
        if (wall == 1) { //WALL_MAX
            increment = part_size - (bmax[loop_i] - cord_param);
            if (increment > distance_minimum){
                particle.acceleration[loop_i]=particle.acceleration[loop_i] - (stiff_collision * increment + damping * particle.v[loop_i]);
            }
        }
    }
}
void general_boundary_collision(Vect3<int> belongings, Particle particle){
    Vect3<int> zero(0,0,0);
    if (belongings == zero) {
        return;
    }
    new_boundary_collision(belongings, particle);
}



void new_boundary_interaction(Vect3<int> belongings, Particle &particle) {
    int wall = 0;
    double distance = 0;
    double posible_pos=0;
    Vect3<double> bmin(bmin_coord_x, bmin_coord_y, bmin_coord_z);
    Vect3<double> bmax(bmax_coord_x, bmax_coord_y, bmax_coord_z);
    for (int loop_i = 0; loop_i < 3; loop_i++) {
        wall = belongings[loop_i];
        if (wall== -1) { //WALL_MIN
            distance = particle.pos[loop_i] - bmin[loop_i];
            posible_pos=bmin[loop_i]-distance;
        }
        else if (wall == 1) { //WALL_MAX
            distance = bmax[loop_i] - particle.pos[loop_i];
            posible_pos=distance+bmax[loop_i];
        }

        if ( distance < 0) {
            particle.pos[loop_i]=posible_pos;
            particle.v[loop_i]=-particle.v[loop_i] ;
            particle.hv[loop_i]=-particle.hv[loop_i] ;
        }
    }
}
void general_boundary_interaction(Vect3<int> belongings, Particle particle){
    Vect3<int> zero(0,0,0);
    if (belongings == zero) {
        return;
    }
    new_boundary_interaction(belongings, particle);
}