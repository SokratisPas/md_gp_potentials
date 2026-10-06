#include <iostream>
#include <filesystem>
#include <fstream>

#include "particle-system.h"
#include "vec3.h"
#include "data-parser.h"

/*
TO DO:
- include in gitignore the executable file.
- check the input in main

*/

int main(int argc, char* argv[]) 
{    
    // -----------------------------------
    // Parameters 
    std::filesystem::path projectRoot;
    std::filesystem::path current = std::filesystem::current_path();

    if (argc > 0 && argv[0] != nullptr)
    {
        std::filesystem::path executable = std::filesystem::absolute(argv[0]);
        projectRoot = executable.parent_path();
    }

    if (projectRoot.empty() || !std::filesystem::exists(projectRoot / "input" / "dataExample.txt"))
    {
        std::filesystem::path probe = current;
        for (int i = 0; i < 10; ++i)
        {
            if (std::filesystem::exists(probe / "input" / "dataExample.txt"))
            {
                projectRoot = probe;
                break;
            }

            const std::filesystem::path parent = probe.parent_path();
            if (parent == probe)
                break;
            probe = parent;
        }
    }

    if (projectRoot.empty())
    {
        throw std::runtime_error("Could not locate project root.");
    }

    const std::filesystem::path inputFilePath = projectRoot / "input" / "dataExample.txt";
    constexpr char OUTPUT_FILE_NAME[] = "lj.dat";  // output file name
    constexpr int TOTAL_STEPS = 500;        // total steps
    constexpr double dt = 0.0005;           // step time
    constexpr int writeEvery = 10;          // write to file every n steps
    
    constexpr double Ttarget = 1.0;       // target temperature (lj)
    constexpr double tau = 0.1;           // relaxation strength 
    
    const double rc = 2.5;            // r cutoff 
    const double rSkin = 0.3;         // r skin
    const double rlist = rc + rSkin;  // neighbor list radius
    const int updateListEvery = 50;   // update neighbor list extra


    XYZParser parser(rlist);
    std::vector<Snapshot> snapshots = parser.parse(inputFilePath.string());
        
    ParticleSystem particleSystem(snapshots); 
    
    // ==== initialize parameteres ====
    particleSystem.initializeParameteres();
    
    // remove net momentum
    // it should : sum(vel[i] == 0)
    particleSystem.removeNetMomentum();

    // volume, density calculation 
    particleSystem.calcVolumeDensity();

    // ==== neighbor list =====

    // apply pbc
    for (auto& atom : particleSystem.snapshots[0].atoms)
    {
        particleSystem.applyPBC(atom.position);
    }

    // initalize neighbour list and forces
    particleSystem.buildNeighborList(rlist);    // use rlist   
    particleSystem.computeForces(rc);           // use rc not rlist

    // initiallize acceleration (we first need to compute forces)
    particleSystem.initiallizeAcc();

    // ======= open file ===========     
    const std::filesystem::path outputDir = projectRoot / "output";
    std::filesystem::create_directories(outputDir);

    std::ofstream meltFile(outputDir / OUTPUT_FILE_NAME);

    if (!meltFile)
    {
        throw std::runtime_error("Failed to open file: " + std::string(OUTPUT_FILE_NAME));
    }

    meltFile << "step,time,Ekin,Epot,Etot,Temperature,Pressure \n";


    // ------------------------------------
    // RUN
    double time = 0.0;
    int step = 0;

    while (step < TOTAL_STEPS)
    {
        // reset Ekin, Epot, Virial ( Sum(r_ij * F_ij) ) 
        particleSystem.resetEkinEpotVirial();

        // update pos wth velocity Verlet method
        particleSystem.updatePos(dt);
        
        // rebuild neighbor list every updateListEvery steps
        if (step % updateListEvery == 0)
        {
            particleSystem.buildNeighborList(rlist);
        }
        // check if rebuild neighboorList
        else if (particleSystem.needRebuild(rSkin))
        {
            particleSystem.buildNeighborList(rlist);
        }

        // update forces 
        particleSystem.computeForces(rc);

        // acceleration update
        particleSystem.updateAcc();

        // vel and Ekin update
        particleSystem.updateVelEkin(dt);

        // Temperature
        particleSystem.computeTemperature();

        // calc virial
        particleSystem.calcVirial();

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