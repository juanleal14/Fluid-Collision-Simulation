#include "progargs.hpp"
#include <math.h>


class Initial_Values {
public:
    // Constructors
    Initial_Values() : ppm(0.0), np(0), m(0.0), h(0.0) {}

    Initial_Values(float ppmValue, int npValue, double mValue, double hValue)
            : ppm(ppmValue), np(npValue), m(mValue), h(hValue) {}

    // Getter methods
    [[nodiscard]] float getPpm() const { return ppm; }
    [[nodiscard]] int getNp() const { return np; }
    [[nodiscard]] double getM() const { return m; }
    [[nodiscard]] double getH() const { return h; }

    // Setter methods
    void setPpm(float ppmValue) { ppm = ppmValue; }
    void setNp(int npValue) { np = npValue; }
    void setM(double mValue) { m = mValue; }
    void setH(double hValue) { h = hValue; }
private:
    float ppm;
    int np;
    double m;
    double h;
};



void check_command_errors(int argc,std::vector<std::string> arguments) {
    if (argc != 4) {
        std::cerr << "Error: Invalid number of arguments: " << argc << ".";
        exit(-1);
    }
    //Is integer
    try {
        // Convert the input parameter to an integer
        int const value = std::stoi(arguments[0]);
        // Use the integer value as needed
        std::cout << "The entered integer is: " << value << "\n";
        if (value < 0) {
            std::cerr << "Error: The entered integer is negative.\n";
            exit(-2);
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: The entered parameter is not an integer.\n";
        exit(-1);
    }
}


Initial_Values read_general_info(std::ifstream &file){
    // Read ppm and np
    //cap 10-11
    Initial_Values initialValues;
    float read_value = NAN;
    int read_value2 = 0;
    file.read(reinterpret_cast<char*>(&read_value), sizeof(float));//NOLINT
    initialValues.setPpm(read_value);
    file.read(reinterpret_cast<char*>(&read_value2), sizeof(int));//NOLINT
    initialValues.setNp(read_value2);
    initialValues.setM(global_density/pow(initialValues.getPpm(),3));
    initialValues.setH(const_r/initialValues.getPpm());
    std::cout << "ppm: " << initialValues.getPpm() << ", np: " << initialValues.getNp() << std::endl;
    return initialValues;
}

Particle read_particle(std::ifstream & file){ ///Read only one particle from file
    float value = 0;
    Particle particle;
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.pos.set_x(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.pos.set_y(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.pos.set_z(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.hv.set_x(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.hv.set_y(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.hv.set_z(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.v.set_x(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.v.set_y(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.v.set_z(static_cast<double>(value));

    return particle;
}




