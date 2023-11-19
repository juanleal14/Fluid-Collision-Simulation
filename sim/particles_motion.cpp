#include "particles_motion.hpp"


void load_trace(std::string trz, Grid &grid_trz, Initial_Values &initialValues){
    std::ifstream file(trz, std::ios::binary);
    if (!file.is_open()) { //Check error opening
        std::cout<<"Error: Cannot open trace file: " << trz <<" for reading";
        exit (-1);
    }
    int num_blocks = 0;
    file.read(reinterpret_cast<char*>(&num_blocks), sizeof(int));//NOLINT
    std::vector<Block> const blocks(num_blocks);
    long particles_in_block = 0;
    long part_id = 0;
    std::cout<<"Grid.size = "<<grid_trz.blocks.size()<<" Num blocks = "<<num_blocks<<'\n';
    double value_double = 0;
    Particle particle_p;
    for (int loop_i = 0; loop_i<num_blocks;loop_i++){
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));//NOLINT
        for (int loop_p = 0; loop_p<particles_in_block;loop_p++){//NOLINT
            file.read(reinterpret_cast<char*>(&part_id), sizeof(long));//NOLINT
            int const part_id_int = static_cast<int>(part_id);
            particle_p.id = part_id_int;
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.pos.set_x(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.pos.set_y(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.pos.set_z(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.hv.set_x(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.hv.set_y(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.hv.set_z(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.v.set_x(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.v.set_y(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.v.set_z(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.density = value_double;
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.acceleration.set_x(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.acceleration.set_y(value_double);
            file.read(reinterpret_cast<char*>(&value_double), sizeof(double));//NOLINT
            particle_p.acceleration.set_z(value_double);
            grid_trz.blocks[loop_i].push_back(particle_p);
        }
    }
    std::cout<<"\nTrace loaded\n";
}


/*void compare_accelerations(Particle &p1, Particle &p2){
    long id = p1.id;
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
}*/

void compare_particle(Particle &p1, Particle &p2){
    long id = p1.id;
    bool if_found =0;
    if (p1.density != p2.density){
        std::cout<<"id = "<<id<<" "<<"Particles densities differ, p1.dens = "<<p1.density<<" p2.dens = "<<p2.density<<" Difference = "<<p1.density-p2.density<<'\n';
        if_found =1;
    }if (p1.pos.x() != p2.pos.x()){
        std::cout<<"id = "<<id<<" "<<"Particles x pos differ, p1.pos.x() = "<<p1.pos.x()<<" p2.pos.x() = "<<p2.pos.x()<<" Difference = "<<p1.pos.x()-p2.pos.x()<<'\n';
        if_found =1;
    }if (p1.pos.y() != p2.pos.y()) {
        std::cout << "id = " << id << " " << "Particles y pos differ, p1.pos.y() = " << p1.pos.y() << " p2.pos.y() = "<< p2.pos.y() <<" Difference = "<<p1.pos.y()-p2.pos.y()<< '\n';
        if_found =1;
    }if (p1.pos.z() != p2.pos.z()){
        std::cout<<"id = "<<id<<" "<<"Particles z pos differ, p1.pos.z() = "<<p1.pos.z()<<" p2.pos.z() = "<<p2.pos.z()<<" Difference = "<<p1.pos.z()-p2.pos.z()<<'\n';
        if_found =1;
    }if (p1.hv.x() != p2.hv.x()){
        std::cout<<"id = "<<id<<" "<<"Particles hvx differ, p1.hv.x() = "<<p1.hv.x()<<" p2.hv.x() = "<<p2.hv.x()<<" Difference = "<<p1.hv.x()-p2.hv.x()<<'\n';
        if_found =1;
    }if (p1.hv.y() != p2.hv.y()){
        std::cout<<"id = "<<id<<" "<<"Particles hvy differ, p1.hv.y() = "<<p1.hv.y()<<" p2.hv.y() = "<<p2.hv.y()<<" Difference = "<<p1.hv.y()-p2.hv.y()<<'\n';
        if_found =1;
    }if (p1.hv.z() != p2.hv.z()){
        std::cout<<"id = "<<id<<" "<<"Particles hvz differ, p1.hv.z() = "<<p1.hv.z()<<" p2.hv.z() = "<<p2.hv.z()<<" Difference = "<<p1.hv.z()-p2.hv.z()<<'\n';
        if_found =1;
    }if (p1.v.x() != p2.v.x()){
        std::cout<<"id = "<<id<<" "<<"Particles vx differ, p1.v.x() = "<<p1.v.x()<<" p2.v.x() = "<<p2.v.x()<<" Difference = "<<p1.v.x()-p2.v.x()<<'\n';
        if_found =1;
    }if (p1.v.y() != p2.v.y()) {
        std::cout << "id = " << id << " " << "Particles vy differ, p1.v.y() = " << p1.v.y() << " p2.v.y() = "<< p2.v.y() << " Difference = "<<p1.v.y()-p2.v.y()<<'\n';
        if_found =1;
    }if (p1.v.z() != p2.v.z()){
        std::cout << "id = " << id << " " << "Particles vz differ, p1.v.z() = " << p1.v.z() << " p2.v.z() = "<< p2.v.z() << " Difference = "<<p1.v.z()-p2.v.z()<<'\n';
        if_found =1;
    }if (p1.acceleration.x() != p2.acceleration.x()){
        std::cout<<"id = "<<id<<" "<<"Particles ax differ, p1.a.x() = "<<p1.acceleration.x()<<" p2.a.x() = "<<p2.acceleration.x()<<" Difference = "<<p1.acceleration.x()-p2.acceleration.x()<<'\n';
        if_found =1;
    }if (p1.acceleration.y() != p2.acceleration.y()) {
        std::cout << "id = " << id << " " << "Particles ay differ, p1.a.y() = " << p1.acceleration.y() << " p2.a.y() = "<< p2.acceleration.y() << " Difference = "<<p1.acceleration.y()-p2.acceleration.y()<<'\n';
        if_found =1;
    }if (p1.acceleration.z() != p2.acceleration.z()){
        std::cout << "id = " << id << " " << "Particles az differ, p1.a.z() = " << p1.acceleration.z() << " p2.a.z() = "<< p2.acceleration.z() << " Difference = "<<p1.acceleration.z()-p2.acceleration.z()<<'\n';
        if_found =1;
    }if (if_found==1){
        //exit(-1);
    }
}

Particle find_elem(long id, Block &block){
    ///Only works if both vectors have same size
    for (auto loop_i : block){
        if(id == loop_i.id){
            return loop_i;
        }
    }

    std::cout<<"Id "<<id<<" not found in grid block\n";
    exit(-1);
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
        //exit(-1);
    }

    long particles_in_block;
    long id;
    Particle part;
    double d;
    Vect3<double> a(0,0,0); //Por qué da error?
    double write_value;
    int counter = 0;
    for (auto current_block: grid.blocks) {
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));//NOLINT
        if (particles_in_block!=0)
          std::cout<<"Entering block: "<<counter<<" particles in block = "<<particles_in_block<<'\n';

        if(current_block.size()!=particles_in_block){
            std::cout<<"Number of particles for block "<<" mismatch: "<<"grid["<<"].size() = "<<" particles in block = "<<particles_in_block<<'\n';
            exit(-1);
        }
        for (int i = 0; i < current_block.size();i++){
            file.read(reinterpret_cast<char*>(&id), sizeof(long));//NOLINT
            Particle part_in_our_grid = find_elem(id,current_block);
            //std::cout<<"Particle "<<part_in_our_grid.id<<" in block["<<"] : ";


            part.id = id;

            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.pos.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.pos.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value),sizeof(double));//NOLINT
            part.pos.set_z(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.hv.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.hv.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.hv.set_z(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.v.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.v.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.v.set_z(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.density = write_value;
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.acceleration.set_x(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.acceleration.set_y(write_value);
            file.read(reinterpret_cast<char*>(&write_value), sizeof(double));//NOLINT
            part.acceleration.set_z(write_value);
            /*if((current_particle != part)){
                std::cout << "\nParticles " << current_particle.id << " and " << part.id << " are not the same." << "I am in block: " << counter<<"\n";
                exit(-1);
            }*/
            compare_particle(part_in_our_grid,part);
        }
        counter +=1;
    }
    std::cout<<"\nTrace is equal to current state of the simulation\n";
    file.close();
}

void write_to_file(const std::string& output_file_address,Grid &grid, Initial_Values &initialValues){
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
        for (int loop_i = 0; loop_i < current_block.size(); loop_i++) {//NOLINT
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

        }
    }
    output_file.close();
}
void particles_motion(Grid &grid) {
    for (const auto& current_block: grid.blocks) {
        for (auto particle : current_block) {
            double const move_x = particle.hv.x() * time_step +
                            particle.acceleration.x() * pow(time_step, 2);
            double const move_y = particle.hv.y() * time_step +
                            particle.acceleration.y() * pow(time_step, 2);
            double const move_z = particle.hv.z() * time_step +
                            particle.acceleration.z() * pow(time_step, 2);
            //optimization
            particle.pos.set(particle.pos.x()+move_x,particle.pos.y()+move_y,particle.pos.z()+move_z);
            particle.v.set(particle.hv.x() + (particle.acceleration.x() * time_step)*.5,particle.hv.y() + (particle.acceleration.y() * time_step)*.5,particle.hv.z() + (particle.acceleration.z() * time_step)*.5);
            particle.hv.set(particle.hv.x() + particle.acceleration.x() * time_step,particle.hv.y() + particle.acceleration.y() * time_step,particle.hv.z() + particle.acceleration.z() * time_step);

        }
    }
}



void density_transform(Particle & particle, Initial_Values& initialValues){
    particle.density = (particle.density + pow(initialValues.getH(),6))* (315*initialValues.getM())/(64*std::numbers::pi* pow(initialValues.getH(),9));
}

void increase_d (Particle &p1, Particle &p2, double h){
    double const dist = p1.distance_to(p2); /// ∥pi − ⃗pj∥2
    if (dist < pow(h,2)) {
        double const increment = pow(pow(h, 2) - dist, 3); /// ∆ρij
        p1.density += increment;  /// ρi = ρi + ∆ρij
        p2.density += increment;  /// ρj = ρj + ∆ρij
    }
}


  void densities_increase(Grid &grid, Initial_Values &initialValues){ ///Ordering way 2
    std::vector<int> contiguous_blocks;
    double const h_val = initialValues.getH();
    for (int i_current_b = 0; i_current_b < grid.blocks.size();i_current_b++){ ///Go through all blocks
        contiguous_blocks = get_contiguous_blocks(i_current_b,grid); ///Get contiguous blocks to current block
        for (int get_cont_block = 2;  get_cont_block < contiguous_blocks.size(); get_cont_block++){ ///Traverse the contiguous blocks
            int const cont_block = contiguous_blocks[get_cont_block];
            for (int particle_current = 0; particle_current < grid[i_current_b].size(); particle_current++) {
              for (int particle_cont = 0; particle_cont < grid[cont_block].size(); particle_cont++) {  /// Go through each particle in the contiguous block
                increase_d(grid[i_current_b][particle_current], grid[cont_block][particle_cont],h_val);
              }
            }
        }
        /// Traverse particle inside current block
        for (int particle_current = 0; particle_current < grid[i_current_b].size(); particle_current++){
            for (int second_particle = particle_current + 1; second_particle < grid[i_current_b].size(); second_particle++){
              increase_d(grid[i_current_b][particle_current], grid[i_current_b][second_particle], h_val);
            }

            //Aqui meter resto de funciones
            //density_transform(grid[i_current_b][particle_current],initialValues); /// We can apply directly density transform

        }
    }
}

void increase_a (Particle& p1, Particle& p2, double h, double m){
    double const dist_sqrd = p1.distance_to(p2); /// ∥pi − ⃗pj∥2
    if (dist_sqrd < pow(h,2)){
        //std::cout<<"Entered i = "<<p1.id<<" j = "<<p2.id<<'\n';
        double const distij = sqrt(std::max(dist_sqrd,pow(10,-12)));
        //Vect3<double> const increment = ((p1.pos - p2.pos) * (15/(std::numbers::pi * pow(h,6))) * ((3*m*stiff_pressure)/2) * (pow(h-updated_dist,2)/updated_dist) * (p1.density + p2.density - 2* global_density) + (p2.v - p1.v) * (45/(std::numbers::pi * pow(h,6))) * viscosity * m)/(p1.density * p2.density);
        //                            = ((pi.px - pj.px) * (15 / (M_PI*pow(h,6))) *                      (3 * m * ps/2) *             pow(h-distij,2)/distij *                  (densities[particle_i_index] + densities[particle_j_index] - 2*p) + (pj.vx - pi.vx) * (45/(M_PI*pow(h,6)) )* nu * m)/(densities[particle_i_index] * densities[particle_j_index]);
        //((pi.px - pj.px) * (15 / (M_PI*pow(h,6))) * (3 * m * ps/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*p) + (pj.vx - pi.vx) * (45/(M_PI*pow(h,6)) )* nu * m)/(densities[particle_i_index] * densities[particle_j_index]);
        Vect3<double> const increment = ((p1.pos - p2.pos) * (15 / (std::numbers::pi * pow(h,6))) * (3 * m * stiff_pressure/2) * pow(h-distij,2)/distij * (p1.density + p2.density - 2*global_density) + (p2.v - p1.v) *  (45/(std::numbers::pi * pow(h,6))) * viscosity * m)/(p1.density * p2.density);
        //std::cout<<"Increment i = "<<p1.id<<" j = "<<p2.id<<" increment = ["<<increment[0]<<" , "<<increment[1]<<" ' "<<increment[2]<<']';
        p1.acceleration += increment;
        p2.acceleration -= increment;
    }
}

void accelerations_transfer(Grid &grid, Initial_Values &initialValues){ ///Ordering 2 acc transf
    std::vector<int> contiguous_blocks;
    double const h_val = initialValues.getH();
    double const m_val = initialValues.getM();
    for (int i_current_b = 0; i_current_b < grid.blocks.size();i_current_b++){ ///Go through all blocks
        contiguous_blocks = get_contiguous_blocks(i_current_b,grid); ///Get contiguous blocks to current block
        for (int get_cont_block = 2;  get_cont_block < contiguous_blocks.size(); get_cont_block++){ ///Traverse the contiguous blocks
            int const cont_block = contiguous_blocks[get_cont_block];
            for (int particle_current = 0; particle_current < grid[i_current_b].size(); particle_current++) {
              for (int particle_cont = 0; particle_cont < grid[cont_block].size(); particle_cont++) {  /// Go through each particle in the contiguous block
                //std::cout<<"Contiguous: ";
                increase_a(grid[i_current_b][particle_current], grid[cont_block][particle_cont],h_val,m_val);
                //std::cout<<'\n';
              }
            }
        }
        /// Traverse particle inside current block
        for (int particle_current = 0; particle_current < grid[i_current_b].size(); particle_current++){
            for (int second_particle = particle_current + 1; second_particle < grid[i_current_b].size(); second_particle++){
              //std::cout<<"Same: ";
              increase_a(grid[i_current_b][particle_current], grid[i_current_b][second_particle], h_val, m_val);
              //std::cout<<'\n';
            }

            //Aqui meter resto de funciones
        }
    }
}



void accelerations_computation(Grid &grid, Initial_Values &initialValues){
    //initialization of densities and accelerations
    densities_increase(grid,initialValues);
    check_trace("./trz/small/densinc-base-1.trz",grid);
    //check_trace("./trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    //densities_transform(grid,initialValues);
 //   acceleration_transfer(grid,initialValues);

}
/*
void simulate(int nsteps, Grid &grid, Initial_Values initialValues){
    for(int i=0;i < nsteps; i++) {
        //Stages of Simulation
        //Stage 2: Accelerations computation
        accelerations_computation(grid,initialValues);
        //Stage 3: Particle Collisions
        particle_collision(grid);
        //stage 4: Particles motion
        particles_motion(grid);
        //Stage 5: Boundary collisions
        boundary_collision(grid);
    }
}*/
void new_particles_motion(Particle &particle){
        Vect3<double> const move = particle.hv*time_step + particle.acceleration*pow(time_step,2);
        Vect3<double> hvnew = particle.acceleration*time_step;
        Vect3<double> const vnew = hvnew*.5;
        particle.pos += move;
        particle.v = particle.hv + vnew;
        particle.hv += hvnew;
}
/*
void new_simulate(int nsteps, Grid &grid, Initial_Values initialValues){
    for (const auto& current_block: grid.blocks) {
        for (auto particle : current_block) {
            accelerations_computation(grid,initialValues);
        }
    }
    for (const auto& current_block: grid.blocks) {
        for (auto particle : current_block) {
            new_particle_collision(belongs_to_boundary(particle,grid.size),particle);
            new_particles_motion(particle);
            new_boundary_collision(particle);
        }
    }
}*/
