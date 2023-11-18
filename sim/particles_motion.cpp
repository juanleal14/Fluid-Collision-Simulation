#include "particles_motion.hpp"



void load_trace(std::string trz, Grid &grid, std::vector<Particle> &particles, Initial_Values &i_v){
    std::ifstream file(trz, std::ios::binary);
    if (!file.is_open()) { //Check error opening
        std::cout<<"Error: Cannot open trace file: " << trz <<" for reading";
        exit (-1);
    }
    int num_blocks;
    file.read(reinterpret_cast<char*>(&num_blocks), sizeof(int));//NOLINT
    std::vector<std::vector<Particle>> blocks(num_blocks);
    grid.blocks = blocks;
    particles = std::vector<Particle> (i_v.getNp());


    long particles_in_block = 0;
    long part_id = 0;
    std::cout<<"Grid.size = "<<grid.blocks.size()<<" total particles = "<<particles.size()<<'\n';
    for (int i = 0; i<num_blocks;i++){
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));//NOLINT
        for (int p = 0; p<particles_in_block;p++){//NOLINT
            file.read(reinterpret_cast<char*>(&part_id), sizeof(long));//NOLINT
            grid[loop_i].push_back(part_id);
            file.read(reinterpret_cast<char*>(&particles[part_id].pos().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].pos().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].pos().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].hv().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].hv().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].hv().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].v.x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].v.y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].v.z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&densities[part_id]), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&read_value), sizeof(double));//NOLINT
            accelerations[part_id].set_acc_x(read_value);
            file.read(reinterpret_cast<char*>(&read_value), sizeof(double));//NOLINT
            accelerations[part_id].set_acc_y(read_value);
            file.read(reinterpret_cast<char*>(&read_value), sizeof(double));//NOLINT
            accelerations[part_id].set_acc_z(read_value);
        }
    }
    std::cout<<"\nTrace loaded\n";
}


void compare_accelerations(Particle &p1, Particle &p2, long id){
    if (p1.acceleration.x()!=p2.acceleration.x()){
        std::cout<<"id = "<<id<<" "<<"Accelerations x() differ, a1.x() = "<<p1.acceleration.x()<<" a2.x() = "<<p2.acceleration.x()<<'\n';
        //exit(-1);
    }
    if (p1.acceleration.y()!=p2.acceleration.y()){
        std::cout<<"id = "<<id<<" "<<"Accelerations ay differ, a1.acceleration.y() = "<<p1.acceleration.y()<<" a2.acceleration.y() = "<<p2.acceleration.y()<<'\n';
        //exit(-1);
    }
    if (p1.acceleration.z()!=p2.acceleration.z()){
        std::cout<<"id = "<<id<<" "<<"Accelerations az differ, a1.az = "<<p1.acceleration.z()<<" a2.az = "<<p2.acceleration.z()<<'\n';
        //exit(-1);
    }
}

void compare_particle(Particle &p1, Particle &p2,long id){
    if (p1.pos.x() != p2.pos.x()){
        std::cout<<"id = "<<id<<" "<<"Particles x pos differ, p1.pos.x() = "<<p1.pos.x()<<" p2.pos.x() = "<<p2.pos.x()<<'\n';
        //exit(-1);
    }if (p1.pos.y() != p2.pos.y()) {
        std::cout << "id = " << id << " " << "Particles y pos differ, p1.pos.y() = " << p1.pos.y() << " p2.pos.y() = "
                  << p2.pos.y() << '\n';
        //exit(-1);
    }if (p1.pos.z() != p2.pos.z()){
        std::cout<<"id = "<<id<<" "<<"Particles z pos differ, p1.pos.z() = "<<p1.pos.z()<<" p2.pos.z() = "<<p2.pos.z()<<'\n';
        //exit(-1);
    }if (p1.hv.x() != p2.hv.x()){
        std::cout<<"id = "<<id<<" "<<"Particles hvx pos differ, p1.hv.x() = "<<p1.hv.x()<<" p2.hv.x() = "<<p2.hv.x()<<'\n';
        //exit(-1);
    }if (p1.hv.y() != p2.hv.y()){
        std::cout<<"id = "<<id<<" "<<"Particles hvy pos differ, p1.hv.y() = "<<p1.hv.y()<<" p2.hv.y() = "<<p2.hv.y()<<'\n';
        //exit(-1);
    }if (p1.hv.z() != p2.hv.z()){
        std::cout<<"id = "<<id<<" "<<"Particles hvz pos differ, p1.hv.z() = "<<p1.hv.z()<<" p2.hv.z() = "<<p2.hv.z()<<'\n';
        //exit(-1);
    }if (p1.v.x() != p2.v.x()){
        std::cout<<"id = "<<id<<" "<<"Particles vx pos differ, p1.v.x() = "<<p1.v.x()<<" p2.v.x() = "<<p2.v.x()<<'\n';
        //exit(-1);
    }if (p1.v.y() != p2.v.y()) {
        std::cout << "id = " << id << " " << "Particles vy pos differ, p1.v.y() = " << p1.v.y() << " p2.v.y() = "
                  << p2.v.y() << '\n';
        //exit(-1);
    }if (p1.v.z() != p2.v.z()){
        std::cout << "id = " << id << " " << "Particles vz pos differ, p1.v.z() = " << p1.v.z() << " p2.v.z() = "
                  << p2.v.z() << '\n';
        //exit(-1);
    }
    }

void find_elem(long id, Block &block){
    ///Only works if both vectors have same size
    int found = 0;
    for (auto i : block){
        if(id == i.id){
            found = 1;
        }
    }
    if (found==0){
        std::cout<<"Id "<<id<<" not found in grid block\n";
        exit(-1);
    }
}

void check_trace(std::string trz, Grid &grid){
    std::ifstream file(trz, std::ios::binary);
    if (!file.is_open()) { //Check error opening
        std::cout<<"Error: Cannot open trace file: " << trz <<" for reading";
        exit (-1);
    }

    int num_blocks;
    file.read(reinterpret_cast<char*>(&num_blocks), sizeof(int));//NOLINT

    if (num_blocks != grid.blocks.size()){
        std::cout<<"Number of blocks differ from trace:\n"<<"trz_blocks = "<<num_blocks<<'\t'<<"grid_blocks = "<<grid.blocks.size()<<'\n';
        //cout<<"Pruebita"<<(grid.size.nx-1)*(grid.size.ny-1)*(grid.size.nz-1)<<'\n';
        exit(-1);
    }

    long particles_in_block;
    long id;
    Particle part;
    double d;
    Vect3<double> a;
    double write_value;

    for (auto current_block: grid.blocks) {
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));//NOLINT
        std::cout<<"Entering block: "<<" particles in block = "<<particles_in_block<<'\n';

        if(current_block.size()!=particles_in_block){
            std::cout<<"Number of particles for block "<<" mismatch: "<<"grid["<<i<<"].size() = "<<grid[i].size()<<" particles in block = "<<particles_in_block<<'\n';
            exit(-1);
        }
        for (auto current_particle : current_block){
            file.read(reinterpret_cast<char*>(&id), sizeof(long));//NOLINT
            //cout<<"Particle "<<id<<" in block["<<i<<"] : ";
            find_elem(id,current_block);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.pos.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.pos.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value),sizeof(double));//NOLINT
            current_particle.pos.set_z(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.hv.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.hv.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.hv.set_z(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.v.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.v.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.v.set_z(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.density = write_value;
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.acceleration.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.acceleration.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            current_particle.acceleration.set_z(write_value);

            compare_particle(particles[id],part,id);

            if (densities[id]!=d){
                std::cout<<"Densities for particle "<<id<<" differ, d = "<<d<<" densities["<<id<<"] = "<<densities[id]<<" ; Difference = "<<d-densities[id]<<'\n';
                //exit(-1);
            }

            compare_accelerations(accelerations[id],a,id);

        }
    }


    std::cout<<"\nTrace is equal to current state of the simulation\n";

    file.close();
}

void write_to_file(const std::string& output_file_address,Grid grid, Initial_Values &initialValues){
    //Write to the file all the new values
    std::ofstream output_file(output_file_address, std::ios::binary);
    if (!output_file.is_open()) { //Check error opening
        std::cout<<"Error: Cannot open " << output_file_address <<" for writing";
        exit (-4);
    }
    float part_coord_x = 0;
    float part_coord_y = 0;
    float part_coord_z = 0;
    float part_vel_x = 0;
    float part_vel_y = 0;
    float part_vel_z = 0;
    float hvx = 0;
    float hvy = 0;
    float hvz = 0;
    float read_value = initialValues.getPpm();
    output_file.write(reinterpret_cast<char*>(&read_value), sizeof(float));//NOLINT
    int read_value2 = initialValues.getNp();
    output_file.write(reinterpret_cast<char*>(&read_value2), sizeof(int));//NOLINT
    std::cout<<"\n"<<initialValues.getPpm()<<"HOLA"<<initialValues.getNp()<<"\n";
    //create a loop to write all the particles

    for (auto current_block : grid.blocks){
        for (int loop_i = 0; loop_i <= current_block.size(); loop_i++) {//NOLINT
            part_coord_x = static_cast<float>(current_block[loop_i].pos.x());
            part_coord_y = static_cast<float>(current_block[loop_i].pos.y());
            part_coord_z = static_cast<float>(current_block[loop_i].pos.z());
            hvx = static_cast<float>(current_block[loop_i].hv.x());
            hvy = static_cast<float>(current_block[loop_i].hv.y());
            hvz = static_cast<float>(current_block[loop_i].hv.z());
            part_vel_x = static_cast<float>(current_block[loop_i].v.x());
            part_vel_y = static_cast<float>(current_block[loop_i].v.y());
            part_vel_z = static_cast<float>(current_block[loop_i].v.z());
            output_file.write(reinterpret_cast<const char *>(&part_coord_x), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&part_coord_y), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&part_coord_z), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&hvx), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&hvy), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&hvz), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&part_vel_x), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&part_vel_y), sizeof(float));//NOLINT
            output_file.write(reinterpret_cast<const char *>(&part_vel_z), sizeof(float));//NOLINT
            //output_file.write("\n", sizeof(char ));
            std::cout << "\n" <<current_block[loop_i].pos.x() << "\n";
        }
    }
    output_file.close();
}
void particles_motion(Grid &grid) {
    for (auto current_block: grid.blocks) {
        for (auto particle : current_block) {
            double move_x = particle.hv.x() * time_step +
                            particle.acceleration.x() * pow(time_step, 2);
            double move_y = particle.hv.y() * time_step +
                            particle.acceleration.y() * pow(time_step, 2);
            double move_z = particle.hv.z() * time_step +
                            particle.acceleration.z() * pow(time_step, 2);
            //optimization
            int old_block = find_block(particle, grid.size);
            for (auto blocks: grid.blocks) {
                particle.pos.set(particle.pos.x()+move_x,particle.pos.y()+move_y,particle.pos.z()+move_z);
                particle.v.set(particle.hv.x() + (particle.acceleration.x() * time_step)*.5,particle.hv.y() + (particle.acceleration.y() * time_step)*.5,particle.hv.z() + (particle.acceleration.z() * time_step)*.5);
                particle.hv.set(particle.hv.x() + particle.acceleration.x() * time_step,particle.hv.y() + particle.acceleration.y() * time_step,particle.hv.z() + particle.acceleration.z() * time_step);
                /*int new_block = find_block(particle, grid.size);
                if (old_block != new_block) {
                   grid[[new_block].push_back(particle);
                    Particle part_x = grid[old_block].begin();
                    while (part_x != particle) { part_x++; }
                    grid[old_block].erase(part_x);
                }*/
            }
        }
    }
}
/*
void particles_motion(Grid &grid){
    for (auto current_block: grid.blocks) {
        for (int loop_i = 0; loop_i <= current_block.size(); loop_i++) {
        double move_x = loop_i.hv.x() * time_step + accelerations[i].ax * pow(time_step, 2);
        double move_y = particles[i].hvy * time_step + accelerations[i].ay * pow(time_step, 2);
        double move_z = particles[i].hvz * time_step + accelerations[i].az * pow(time_step, 2);

        int old_block = find_block(particles[i], grid.size);

        particles[i].px += move_x;
        particles[i].py += move_y;
        particles[i].pz += move_z;
        particles[i].vx = particles[i].hvx + (accelerations[i].ax * time_step) / 2;
        particles[i].vy = particles[i].hvy + (accelerations[i].ay * time_step) / 2;
        particles[i].vz = particles[i].hvz + (accelerations[i].az * time_step) / 2;
        particles[i].hvx = particles[i].hvx + accelerations[i].ax * time_step;
        particles[i].hvy = particles[i].hvy + accelerations[i].ay * time_step;
        particles[i].hvz = particles[i].hvz + accelerations[i].az * time_step;

        int new_block = find_block(particles[i], grid.size);

        if (old_block != new_block) {
            grid[new_block].push_back(i);

            auto x = grid[old_block].begin();
            while (*x != i) { x++; }

            grid[old_block].erase(x);

        }
    }
    }
}*/
void densities_increase(Grid &grid, Initial_Values initialValues){ /// Cambiar p por part, porque ya hay una varibale gloabl p
    std::vector<Block> contiguous_blocks;
;    for (int i_current_b = 0; i_current_b < grid.blocks.size();i_current_b++){ ///Go through all blocks
        contiguous_blocks = get_contiguous_blocks(i_current_b,grid.size); ///Get contiguous blocks to current block
        for (auto particle_current : grid[i_current_b]){ ///Go through each particle of the current block
            for (const auto& current_cont_block : contiguous_blocks){ ///Traverse the contiguous blocks
                for (auto particle_cont_current : current_cont_block){ /// Go through each particle in the contiguous block
                    if !(particle_current==(particle_cont_current)) { /// Check particle_i != particle_j
                        if (distance_squared(particle_current, particle_cont_current) < (pow(initialValues.getH(), 2))) {
                            grid[i_current_b] += pow(pow(initialValues.getH(), 2) - distance_squared(particle_current, particle_cont_current), 3);
                        }
                    }
                }
            }
        }
    }
}
/*
void densities_increase(std::vector<Particle> &particles, Grid &grid, std::vector<double> &densities){ /// Cambiar p por part, porque ya hay una varibale gloabl p
    std::vector<int> contiguous_blocks;
    int particle_i_index;
    Particle pi;Particle pj;
    int c_block_index;
    int particle_j_index;
    for (int i = 0; i < grid.blocks.size();i++){ ///Go through all blocks
        contiguous_blocks = get_contiguous_blocks(i,grid.size); ///Get contiguous blocks to current block
        for (int p = 0; p < grid[i].size(); p++){ ///Go through each particle of the current block
            particle_i_index = grid[i][p];
            pi = particles[particle_i_index];
            for (int b = 0; b < contiguous_blocks.size(); b++){ ///Traverse the contiguous blocks
                c_block_index = contiguous_blocks[b]; /// Get the index of the contiguous block to traverse
                for (int j = 0; j < grid[c_block_index].size();j++){ /// Go through each particle in the contiguous block
                    particle_j_index = grid[c_block_index][j];
                    pj = particles[particle_j_index];
                    if (particle_i_index != particle_j_index) { /// Check pi != pj
                        if (distance_squared(pi, pj) < (pow(h, 2))) {
                            densities[particle_i_index] += pow(pow(h, 2) - distance_squared(pi, pj), 3);
                        }
                    }
                }
            }
        }
    }
}

*/
void densities_transform(Grid &grid,Initial_Values initialValues){
    for (const auto& current_block: grid.blocks){
        for (auto particle : current_block){
            particle.density = (particle.density + pow(initialValues.getH(),6))* (315*initialValues.getM())/(64*std::numbers::pi* pow(initialValues.getH(),9));
        }
    }
}

void acceleration_transfer(std::vector<Particle> &particles, Grid &grid, Initial_Values initialValues){
    for (int loop_i = 0; loop_i<=grid.blocks.size();loop_i++){ ///Go through all blocks
        std::vector<Block> const contiguous_blocks = get_contiguous_blocks(loop_i,grid.size); ///Get contiguous blocks to current block
        for (int part = 0; part < grid[loop_i].size(); part++){ ///Go through each particle of the current block
            Particle particle_i = grid[loop_i][part];
            for (int loop_b = 0; loop_b < contiguous_blocks.size(); loop_b++){ ///Traverse the contiguous blocks
                int c_block_index = contiguous_blocks[loop_b]; /// Get the index of the contiguous block to traverse
                for (int loop_j = 0; loop_j < grid[c_block_index].size();loop_j++){ /// Go through each particle in the contiguous block
                    int particle_j_index = grid[c_block_index][loop_j];
                    Particle particle_j = particles[particle_j_index];
                    if (particle_i_index != particle_j_index) { /// Check particle_i != particle_j
                        double dist_squared = distance_squared(particle_i, particle_j);
                        if (dist_squared < (pow(h, 2))) {
                            double distij = sqrt(std::max(dist_squared, pow(10,-12))); /// In these 4 lines calculate distij as stated in project and update accelerations
                            grid[i].acceleration.x() += ((particle_i.pos.x() - particle_j.pos.x()) * (15 / (std::numbers::particle_i*pow(initialValues.getH(),6))) * (3 * initialValues.getM() * stiff_pressure/2) * pow(initialValues.getH()-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*global_density) + (particle_j.v.x() - particle_i.v.x()) * (45/(std::numbers::particle_i*pow(initialValues.getH(),6)) ) * viscosity * initialValues.getM()) / (densities[particle_i_index] * densities[particle_j_index]);
                            grid[i].acceleration.y() += ((particle_i.pos.y() - particle_j.pos.y()) * (15 / (std::numbers::particle_i*pow(initialValues.getH(),6))) * (3 * initialValues.getM() * stiff_pressure/2) * pow(initialValues.getH()-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*global_density) + (particle_j.v.y() - particle_i.v.y()) * (45/(std::numbers::particle_i*pow(initialValues.getH(),6)) ) * viscosity * initialValues.getM()) / (densities[particle_i_index] * densities[particle_j_index]);
                            grid[i].acceleration.y() += ((particle_i.pos.z() - particle_j.pos.z()) * (15 / (std::numbers::particle_i*pow(initialValues.getH(),6))) * (3 * initialValues.getM() * stiff_pressure/2) * pow(initialValues.getH()-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*global_density) + (particle_j.v.z() - particle_i.v.z()) * (45/(std::numbers::particle_i*pow(initialValues.getH(),6)) ) * viscosity * initialValues.getM()) / (densities[particle_i_index] * densities[particle_j_index]);
                        }
                    }
                }
            }
        }
    }
}

void accelerations_computation(Grid &grid){
    //initialization of densities and accelerations
    for (auto current_block: grid.blocks) {
        for (auto loop_i : current_block) {
        //for (int loop_i = 0; loop_i < current_block.size(); loop_i++) {
            loop_i.density = 0;
            loop_i.acceleration.set_x(0);
            loop_i.acceleration.set_y(-gravity);
            loop_i.acceleration.set_z(0);
            //accelerations.push_back(a);
        }
    }
    //check_trace("../trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    densities_increase(grid);
    //check_trace("./trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    densities_transform(grid);
    acceleration_transfer(grid);

}

void simulate(int nsteps, Grid &grid){
    for(int i=0;i < nsteps; i++) {
        //std::vector<double> densities;
        //Stages of Simulation
        //Stage 2: Accelerations computation
        accelerations_computation(grid);
        //Stage 3: Particle Collisions
        particle_collision(grid);
        //stage 4: Particles motion
        particles_motion(grid);
        //Stage 5: Boundary collisions
        boundary_collision(grid);
    }
}
