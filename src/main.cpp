#include <iostream>
#include <filesystem>
#include <fstream>

#include "particle-system.h"
#include "vec3.h"
#include "data-parser.h"


int main() 
{   
    constexpr int N = 500;                  // Number of particles    
    ParticleSystem particleSystem(N); 
    
    // -----------------------------------
    // Parameters 
    constexpr char FILE_NAME[] = "lj.dat";  // output file name
    constexpr int TOTAL_STEPS = 500;        // total steps
    constexpr double dt = 0.0005;           // step time
    constexpr int writeEvery = 10;          // write to file every n steps
    
    constexpr double Ttarget = 1.0;       // target temperature (lj)
    constexpr double tau = 0.1;           // relaxation strength 
    constexpr double e = 1.0;             // epsilon (lj units)
    constexpr double sigma = 1.0;         // sigma (lj units)
    const double r0 = std::pow(2.0, 1.0 / 6.0) * sigma;
    const double a = std::sqrt(2.0) * r0;  // lattice constant ( = sqrt(2) * 1.122)


    // ==== initialize parameteres ====
    particleSystem.initializeParameteres(a);
    
    // remove net momentum
    // it should : sum(vel[i] == 0)
    particleSystem.removeNetMomentum();

    // ==== rescale cube =======
    // (rescale cube BEFORE initialize neighborList and computeForces)
    int cells_per_dim = std::ceil(std::cbrt(N / 4.0));
    double L = cells_per_dim * a; // cubes length
    
    particleSystem.xMin = particleSystem.yMin = particleSystem.zMin = 0.0;
    particleSystem.xMax = particleSystem.yMax = particleSystem.zMax = L;

    // volume, density calculation 
    particleSystem.calcVolumeDensity();

    // ==== neighbor list =====
    double rc = 2.5 * sigma;    // r cutoff = 2.5 * 1.0 = 2.5
    double rSkin = 0.3;         // r skin
    double rlist = rc + rSkin;   // neighbor list radius

    // apply pbc
    for (auto& p : particleSystem.particles)
    {
        particleSystem.applyPBC(p.position);
    }

    // initalize neighbour list and forces
    particleSystem.buildNeighborList(rlist);
    particleSystem.computeForces_Epot_Virial(e, sigma);

    // initiallize acceleration
    particleSystem.initiallizeAcc();

    // ======= open file ===========     
    std::filesystem::create_directories("../output");  // creates folder in md-gp-potential/output 
    
    std::ofstream meltFile(std::string("../output/") + FILE_NAME); // creates or opens file

    if (!meltFile)
    {
        throw std::runtime_error("Failed to open file: " + std::string(FILE_NAME));
    }

    meltFile << "step,time,Ekin,Epot,Etot,Temperature,Pressure \n";


    // ------------------------------------
    // RUN
    double time = 0.0;
    int step = 0;

    while (step < TOTAL_STEPS)
    {
        // reset Ekin, Epot, Virial
        particleSystem.resetEkinEpotVirial();

        // update pos wth velocity Verlet method
        particleSystem.updatePos(dt);

        // check if rebuild neighboorList
        if (particleSystem.needRebuild(rSkin))
        {
            particleSystem.buildNeighborList(rlist);
        }

        // update forces and E potential and Virial
        particleSystem.computeForces_Epot_Virial(e, sigma);

        // acceleration update
        particleSystem.updateAcc();

        // vel and Ekin update
        particleSystem.updateVelEkin(dt);

        // Temperature
        particleSystem.computeTemperature();

        // Pressure
        particleSystem.computePressure();

        // thermostat
        particleSystem.useThermostat(Ttarget, tau, dt);

        // write to file
        if (step % writeEvery == 0)
        {
            double Etot = particleSystem.Ekin + particleSystem.Epot;

            meltFile << step << ","
                << time << ","
                << particleSystem.Ekin << ","
                << particleSystem.Epot << ","
                << Etot << ","
                << particleSystem.temperature << ","
                << particleSystem.pressure << "\n";
        }

        step++;
        time += dt;
    }

    meltFile.close();

    return 0;
}