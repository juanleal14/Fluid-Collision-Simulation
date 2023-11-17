#include "particles_motion.hpp"
#include <numbers>


void load_trace(std::string trz, Grid &grid, std::vector<Particle> &particles, std::vector<double> &densities, std::vector<Acceleration> &accelerations, Initial_Values &i_v){
    std::ifstream file(trz, std::ios::binary);
    if (!file.is_open()) { //Check error opening
        std::cout<<"Error: Cannot open trace file: " << trz <<" for reading";
        exit (-1);
    }
    int num_blocks;
    file.read(reinterpret_cast<char*>(&num_blocks), sizeof(int));//NOLINT
    std::vector<std::vector<int>> blocks(num_blocks);
    grid.blocks = blocks;
    particles = std::vector<Particle> (i_v.getNp());
    densities = std::vector<double> (i_v.getNp());
    accelerations = std::vector<Acceleration> (i_v.getNp());
    long particles_in_block = 0;
    long part_id = 0;
    std::cout<<"Grid.size = "<<grid.blocks.size()<<" total particles = "<<particles.size()<<'\n';
    for (int i = 0; i<num_blocks;i++){
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));//NOLINT
        for (int p = 0; p<particles_in_block;p++){//NOLINT
            file.read(reinterpret_cast<char*>(&part_id), sizeof(long));//NOLINT
            grid.blocks[loop_i].push_back(part_id);
            file.read(reinterpret_cast<char*>(&particles[part_id].pos().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].pos().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].pos().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].hv().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].hv().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].hv().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].vel().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].vel().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&particles[part_id].vel().z()), sizeof(double));//NOLINT
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


void compare_accelerations(Acceleration &a1, Acceleration &a2, long id){
    if (a1.x()!=a2.x()){
        std::cout<<"id = "<<id<<" "<<"Accelerations x() differ, a1.x() = "<<a1.x()<<" a2.x() = "<<a2.x()<<'\n';
        //exit(-1);
    }
    if (a1.ay!=a2.ay){
        std::cout<<"id = "<<id<<" "<<"Accelerations ay differ, a1.ay = "<<a1.ay<<" a2.ay = "<<a2.ay<<'\n';
        //exit(-1);
    }
    if (a1.az!=a2.az){
        std::cout<<"id = "<<id<<" "<<"Accelerations az differ, a1.az = "<<a1.az<<" a2.az = "<<a2.az<<'\n';
        //exit(-1);
    }
}

void compare_particle(Particle &p1, Particle &p2,long id){
    if (p1.pos().x() != p2.pos().x()){
        std::cout<<"id = "<<id<<" "<<"Particles x pos differ, p1.pos().x() = "<<p1.pos().x()<<" p2.pos().x() = "<<p2.pos().x()<<'\n';
        //exit(-1);
    }
    if (p1.pos().y() != p2.pos().y()){
        std::cout<<"id = "<<id<<" "<<"Particles y pos differ, p1.pos().y() = "<<p1.pos().y()<<" p2.pos().y() = "<<p2.pos().y()<<'\n';
        //exit(-1);
    }
    if (p1.pos().z() != p2.pos().z()){
        std::cout<<"id = "<<id<<" "<<"Particles z pos differ, p1.pos().z() = "<<p1.pos().z()<<" p2.pos().z() = "<<p2.pos().z()<<'\n';
        //exit(-1);
    }
    if (p1.hv().x() != p2.hv().x()){
        std::cout<<"id = "<<id<<" "<<"Particles hvx pos differ, p1.hv().x() = "<<p1.hv().x()<<" p2.hv().x() = "<<p2.hv().x()<<'\n';
        //exit(-1);
    }
    if (p1.hv().y() != p2.hv().y()){
        std::cout<<"id = "<<id<<" "<<"Particles hvy pos differ, p1.hv().y() = "<<p1.hv().y()<<" p2.hv().y() = "<<p2.hv().y()<<'\n';
        //exit(-1);
    }
    if (p1.hv().z() != p2.hv().z()){
        std::cout<<"id = "<<id<<" "<<"Particles hvz pos differ, p1.hv().z() = "<<p1.hv().z()<<" p2.hv().z() = "<<p2.hv().z()<<'\n';
        //exit(-1);
    }
    if (p1.vel().x() != p2.vel().x()){
        std::cout<<"id = "<<id<<" "<<"Particles vx pos differ, p1.vel().x() = "<<p1.vel().x()<<" p2.vel().x() = "<<p2.vel().x()<<'\n';
        //exit(-1);
    }
    if (p1.vel().y() != p2.vel().y()){
        std::cout<<"id = "<<id<<" "<<"Particles vy pos differ, p1.vel().y() = "<<p1.vel().y()<<" p2.vel().y() = "<<p2.vel().y()<<'\n';
        //exit(-1);
    }
    if (p1.vel().z() != p2.vel().z()){
        std::cout<<"id = "<<id<<" "<<"Particles vz pos differ, p1.vel().z() = "<<p1.vel().z()<<" p2.vel().z() = "<<p2.vel().z()<<'\n';
        //exit(-1);
    }

}

void find_elem(int e, std::vector<int> &v){
    ///Only works if both vectors have same size
    int found = 0;
    for (int i = 0; i<v.size();i++){
        if(e == v[i]){
            found = 1;
        }
    }
    if (found==0){
        std::cout<<"Id "<<e<<" not found in grid block\n";
        //exit(-1);
    }
}

void check_trace(std::string trz, Grid &grid, std::vector<Particle> &particles, std::vector<double> &densities, std::vector<Acceleration> &accelerations){
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
    Acceleration a;

    for (int i = 0; i < num_blocks;i++){
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));//NOLINT
        std::cout<<"Entering block: "<<i<<" particles in block = "<<particles_in_block<<'\n';

        if(grid.blocks[i].size()!=particles_in_block){
            std::cout<<"Number of particles for block "<<i<<" mismatch: "<<"grid.blocks["<<i<<"].size() = "<<grid.blocks[i].size()<<" particles in block = "<<particles_in_block<<'\n';
            exit(-1);
        }

        for (int p = 0; p<particles_in_block;p++){
            file.read(reinterpret_cast<char*>(&id), sizeof(long));//NOLINT
            //cout<<"Particle "<<id<<" in block["<<i<<"] : ";
            find_elem(id,grid.blocks[i]);
            file.read(reinterpret_cast<char*>(&part.pos().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.pos().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.pos().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.hv().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.hv().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.hv().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.vel().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.vel().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&part.vel().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&d), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&a.x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&a.ay), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char*>(&a.az), sizeof(double));//NOLINT


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

std::vector <Particle> receive_trace(std::string trz, std::vector<double> &densities, std::vector<Acceleration> &accelerations) {
    std::ifstream file(trz, std::ios::binary);
    if (!file.is_open()) { //Check error opening
        std::cout << "Error: Cannot open trace file: " << trz << " for reading";
        exit(-1);
    }

    int num_blocks;
    Grid grid;
    std::vector <Particle> particles;
    file.read(reinterpret_cast<char *>(&num_blocks), sizeof(int));//NOLINT
    long particles_in_block;
    long id;
    Particle part;
    double d;
    Acceleration a;
    std::cout<<"\n"<<num_blocks<<"\n";
    for (int i = 0; i < num_blocks; i++) {
        file.read(reinterpret_cast<char *>(&particles_in_block), sizeof(long));//NOLINT
        std::cout << "Entering block: " << i << " particles in block = " << particles_in_block << '\n';

        for (int p = 0; p < particles_in_block; p++) {
            file.read(reinterpret_cast<char *>(&id), sizeof(long));//NOLINT
            file.read(reinterpret_cast<char *>(&part.pos().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.pos().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.pos().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.hv().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.hv().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.hv().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.vel().x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.vel().y()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&part.vel().z()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&d), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&a.x()), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&a.ay), sizeof(double));//NOLINT
            file.read(reinterpret_cast<char *>(&a.az), sizeof(double));//NOLINT
            accelerations.push_back(a);
            densities.push_back(d);
            particles.push_back(part);
        }
    }
    return particles;
}

void write_to_file(const std::string& output_file_address,std::vector<Particle> particles, Initial_Values &initialValues){
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
    for (int loop_i = 0; loop_i < particles.size(); loop_i++){//NOLINT
        part_coord_x = static_cast<float>(particles[loop_i].pos().x());
        part_coord_y = static_cast<float>(particles[loop_i].pos().y());
        part_coord_z = static_cast<float>(particles[loop_i].pos().z());
        hvx = static_cast<float>(particles[loop_i].hv().x());
        hvy = static_cast<float>(particles[loop_i].hv().y());
        hvz = static_cast<float>(particles[loop_i].hv().z());
        part_vel_x = static_cast<float>(particles[loop_i].vel().x());
        part_vel_y = static_cast<float>(particles[loop_i].vel().y());
        part_vel_z = static_cast<float>(particles[loop_i].vel().z());
        output_file.write(reinterpret_cast<const char*>(&part_coord_x), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&part_coord_y), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&part_coord_z), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&hvx), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&hvy), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&hvz), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&part_vel_x), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&part_vel_y), sizeof(float));//NOLINT
        output_file.write(reinterpret_cast<const char*>(&part_vel_z), sizeof(float));//NOLINT
        //output_file.write("\n", sizeof(char ));
        std::cout<<"\n"<<particles[loop_i].pos().x()<<"\n";
    }
    output_file.close();
}
void particles_motion(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations){
    for (int loop_i = 0; loop_i < particles.size(); loop_i++){

        double move_x = particles[loop_i].hv().x()*time_step + accelerations[loop_i].x()*pow(time_step,2);
        double move_y = particles[loop_i].hv().y()*time_step + accelerations[loop_i].ay*pow(time_step,2);
        double move_z = particles[loop_i].hv().z()*time_step + accelerations[loop_i].az*pow(time_step,2);

        int old_block = find_block(particles[loop_i],grid.size);

        particles[loop_i].pos().x() += move_x;
        particles[loop_i].pos().y() += move_y;
        particles[loop_i].pos().z() += move_z;
        particles[loop_i].vel().x() = particles[loop_i].hv().x() + (accelerations[loop_i].x()*time_step)/2;
        particles[loop_i].vel().y() = particles[loop_i].hv().y() + (accelerations[loop_i].ay*time_step)/2;
        particles[loop_i].vel().z() = particles[loop_i].hv().z() + (accelerations[loop_i].az*time_step)/2;
        particles[loop_i].hv().x() = particles[loop_i].hv().x() + accelerations[loop_i].x()*time_step;
        particles[loop_i].hv().y() = particles[loop_i].hv().y() + accelerations[loop_i].ay*time_step;
        particles[loop_i].hv().z() = particles[loop_i].hv().z() + accelerations[loop_i].az*time_step;

        int new_block = find_block(particles[loop_i],grid.size);

        if (old_block!=new_block){
            grid.blocks[new_block].push_back(loop_i);

            auto part_x = grid.blocks[old_block].begin();
            while ( *part_x != loop_i){part_x++;}

            grid.blocks[old_block].erase(part_x);
        }

    }
}

void densities_increase(std::vector<Particle> &particles, Grid &grid, std::vector<double> &densities){ /// Cambiar p por part, porque ya hay una varibale gloabl p
    std::vector<int> contiguous_blocks;
    int particle_i_index;
    Particle pi;Particle pj;
    int c_block_index;
    int particle_j_index;
    for (int i = 0; i < grid.blocks.size();i++){ ///Go through all blocks
        contiguous_blocks = get_contiguous_blocks(i,grid.size); ///Get contiguous blocks to current block
        for (int p = 0; p < grid.blocks[i].size(); p++){ ///Go through each particle of the current block
            particle_i_index = grid.blocks[i][p];
            pi = particles[particle_i_index];
            for (int b = 0; b < contiguous_blocks.size(); b++){ ///Traverse the contiguous blocks
                c_block_index = contiguous_blocks[b]; /// Get the index of the contiguous block to traverse
                for (int j = 0; j < grid.blocks[c_block_index].size();j++){ /// Go through each particle in the contiguous block
                    particle_j_index = grid.blocks[c_block_index][j];
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

void densities_transform(std::vector<double> &densities){
    for (int i = 0; i < densities.size(); i++){
        densities[i] = (densities[i] + pow(h,6))* (315*m)/(64*std::numbers::pi* pow(h,9));
    }
}
void acceleration_transfer(std::vector<Particle> &particles, Grid &grid, std::vector<double> &densities, std::vector<Acceleration> &accelerations){
    for (int i = 0; i < grid.blocks.size();i++){ ///Go through all blocks
        std::vector<int> contiguous_blocks = get_contiguous_blocks(i,grid.size); ///Get contiguous blocks to current block
        for (int part = 0; part < grid.blocks[i].size(); part++){ ///Go through each particle of the current block
            int particle_i_index = grid.blocks[i][part];
            Particle pi = particles[particle_i_index];
            for (int b = 0; b < contiguous_blocks.size(); b++){ ///Traverse the contiguous blocks
                int c_block_index = contiguous_blocks[b]; /// Get the index of the contiguous block to traverse
                for (int j = 0; j < grid.blocks[c_block_index].size();j++){ /// Go through each particle in the contiguous block
                    int particle_j_index = grid.blocks[c_block_index][j];
                    Particle pj = particles[particle_j_index];
                    if (particle_i_index != particle_j_index) { /// Check pi != pj
                        double dist_squared = distance_squared(pi, pj);
                        if (dist_squared < (pow(h, 2))) {
                            double distij = sqrt(std::max(dist_squared, pow(10,-12))); /// In these 4 lines calculate distij as stated in project and update accelerations
                            accelerations[particle_i_index].x() += ((pi.pos().x() - pj.pos().x()) * (15 / (std::numbers::pi*pow(h,6))) * (3 * m * stiff_pressure/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*global_density) + (pj.vel().x() - pi.vel().x()) * (45/(std::numbers::pi*pow(h,6)) ) * viscosity * m) / (densities[particle_i_index] * densities[particle_j_index]);
                            accelerations[particle_i_index].ay += ((pi.pos().y() - pj.pos().y()) * (15 / (std::numbers::pi*pow(h,6))) * (3 * m * stiff_pressure/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*global_density) + (pj.vel().y() - pi.vel().y()) * (45/(std::numbers::pi*pow(h,6)) ) * viscosity * m) / (densities[particle_i_index] * densities[particle_j_index]);
                            accelerations[particle_i_index].az += ((pi.pos().z() - pj.pos().z()) * (15 / (std::numbers::pi*pow(h,6))) * (3 * m * stiff_pressure/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*global_density) + (pj.vel().z() - pi.vel().z()) * (45/(std::numbers::pi*pow(h,6)) ) * viscosity * m) / (densities[particle_i_index] * densities[particle_j_index]);
                        }
                    }
                }
            }
        }
    }
}

void accelerations_computation(std::vector<Particle> &particles, Grid &grid,std::vector <double> &densities, std::vector <Acceleration> &accelerations){
    //initialization of densities and accelerations
    densities.clear();
    accelerations.clear();
    struct Acceleration a;
    for (int i = 0; i < particles.size(); i++){
        densities.push_back(0);
        a.ax = 0;
        a.ay = -gravity;
        a.az = 0;
        accelerations.push_back(a);
    }

    //check_trace("../trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    densities_increase(particles,grid,densities);
    //check_trace("./trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    densities_transform(densities);
    acceleration_transfer(particles,grid,densities,accelerations);

}

void simulate(int nsteps, std::vector<Particle> &particles, Grid &grid){
    for(int i=0;i < nsteps; i++) {
        std::vector<double> densities;
        std::vector<Acceleration> accelerations;
        //Stages of Simulation
        //Stage 2: Accelerations computation
        accelerations_computation(particles, grid, densities, accelerations);
        //Stage 3: Particle Collisions
        particle_collision(particles,grid,accelerations);
        //stage 4: Particles motion
        particles_motion(particles, grid, accelerations);
        //Stage 5: Boundary collisions
        boundary_collision(particles,grid);
    }
}
