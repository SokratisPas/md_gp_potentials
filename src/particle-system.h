#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <functional>

#include "vec3.h"
#include "data-parser.h"

// =======================================
// Particle System Class
class ParticleSystem {
public:
    int partNumber = 1;         // Number of particles

    double Ekin = 0.0;          // Kinetik energy
    double Epot = 0.0;          // potential energy
    double temperature = 0.0;   // temp
    double pressure = 0.0;      // Pressure
    double virial = 0.0;        // Virial Pot Energy
    double volume = 1.0;        // volume
    double density = 1.0;       // Density
    int dof = 1;                // degrees of Freedom
    std::function<double(double)> forceFunc;
    std::function<double(double)> potentialFunc;

    // box dimensions { a11, b22, c33 }
    Vec3 boxVec = Vec3(0.0);

    std::vector<Snapshot> snapshots; // vector of snapshots
    std::vector<Vec3> lastBuildPos; // last positions when neighborList was built

    // Constructor
    explicit ParticleSystem(
        std::vector<Snapshot> snapshots,
        std::function<double(double)> forceFunc,
        std::function<double(double)> potentialFunc
    )
        :forceFunc(forceFunc), potentialFunc(potentialFunc), snapshots(std::move(snapshots))
    {
        if (this->snapshots.empty())
        {
            throw std::runtime_error("Error: Snapshots vector is empty.");
        }

        partNumber = this->snapshots[0].numberOfAtoms;
        dof = 3 * partNumber - 3; // degrees of freedom (we remove 3 bcz we remove net momentum)
    }
    
    double randomDouble(double start, double finish);
    void initializeParameteres();
    void removeNetMomentum();
    void applyPBC(Vec3& pos);
    Vec3 computeLJForce(
        const Vec3& atom1Pos,
        const Vec3& atom2Pos,
        const double cutoff,
        std::function<double(double)> forceFunc);
    void computeForces(const double cutoff);
    void computeEpot(const double cutoff);
    void calcVolumeDensity();
    void initiallizeAcc();
    void resetEkinEpotVirial();
    void updatePos(const double dt);
    bool needRebuild(const double skin);
    void updateAcc();
    void updateVelEkin(const double dt);
    void computeTemperature();
    void computePressure();
    void useThermostat(const double Ttarget, const double tau, const double dt);
    void buildNeighborList(const double rlist);
    void calcVirial();
};



// ----------------------------------------
// calc random double in [start, finish]
double ParticleSystem::randomDouble(double start, double finish)
{
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(start, finish);
    return dist(gen);
}

// ---------------------------------------------------
// initialize vel, box 
void ParticleSystem::initializeParameteres() 
{
    for (auto& atom : snapshots[0].atoms)
    {
        atom.velocity = Vec3(randomDouble(-0.1, 0.1), 
                                    randomDouble(-0.1, 0.1),
                                    randomDouble(-0.1, 0.1));        
    }

    // initialize box, the lattice is orthogonal so only a11, b22, c33
    boxVec = std::move(snapshots[0].box);

    // initialize lastBuildPos
    for (int i = 0; i < partNumber; ++i)
    {
        lastBuildPos.push_back(snapshots[0].atoms[i].position);
    }
}
    
// ----------------------------------------
// remove net momentum
// it should : sum(vel[i] == 0)
void ParticleSystem::removeNetMomentum()
{
    Vec3 velMean(0.0);
    for (auto& p : snapshots[0].atoms)
        velMean = velMean + p.velocity;

    velMean = velMean / static_cast<double>(snapshots[0].atoms.size());

    for (auto& p : snapshots[0].atoms)
        p.velocity = p.velocity - velMean;
}
    
// ----------------------------------------
// aplly PBC 
void ParticleSystem::applyPBC(Vec3& pos)
{
    // keep x, y, x in cube boundaries 
    pos.x = std::fmod(std::fmod(pos.x, boxVec.x) + boxVec.x, boxVec.x);
    pos.y = std::fmod(std::fmod(pos.y, boxVec.y) + boxVec.y, boxVec.y);
    pos.z = std::fmod(std::fmod(pos.z, boxVec.z) + boxVec.z, boxVec.z);
}    

// ----------------------------------------
// particle-particle 
Vec3 ParticleSystem::computeLJForce(
    const Vec3& atom1Pos,
    const Vec3& atom2Pos,
    const double cutoff,
    std::function<double(double)> forceFunc
)
{
    // avoid division by zero 
    Vec3 rij = atom1Pos - atom2Pos;
    double distance = rij.length();

    if (distance == 0.0)
        return Vec3(0.0);
        
    // cut off radius
    if (distance > cutoff)
        return Vec3(0.0);


    double forceMagnitude = forceFunc(distance);

    return rij * (forceMagnitude / distance); // force vector
}

// ----------------------------------------
// calc all lj forces and Epot and Virial 
void ParticleSystem::computeForces(const double cutoff)
{
    for (int i = 0; i < partNumber; ++i)
    {
        snapshots[0].atoms[i].force = Vec3(0.0);

        for (auto& neighbor : snapshots[0].neighborLists[i])
        {
            Vec3 ijForce = computeLJForce(
                snapshots[0].atoms[i].position,
                neighbor.position,
                cutoff,
                forceFunc
            );

            snapshots[0].atoms[i].force = snapshots[0].atoms[i].force + ijForce;
        }
    }
}

// ----------------------------------------
// calculate Lennard-Jones potential energy
void ParticleSystem::computeEpot(const double cutoff)
{
    Epot = 0.0;

    const auto& atoms = snapshots[0].atoms;
    for (int i = 0; i < partNumber; ++i)
    {
        for (int j = i + 1; j < partNumber; ++j)
        {
            Vec3 rij = atoms[i].position - atoms[j].position;

            // Minimum image convention
            rij.x -= boxVec.x * std::round(rij.x / boxVec.x);
            rij.y -= boxVec.y * std::round(rij.y / boxVec.y);
            rij.z -= boxVec.z * std::round(rij.z / boxVec.z);

            const double distance = rij.length();
            if (distance > cutoff)
                continue;
            if (distance == 0.0)
                throw std::runtime_error("Cannot compute Lennard-Jones energy at zero separation.");

            Epot += potentialFunc(distance);
        }
    }
}

// ----------------------------------------
// calc volume, density
void ParticleSystem::calcVolumeDensity()
{
    volume = boxVec.x * boxVec.y * boxVec.z;
    density = partNumber / volume;
}

// ----------------------------------------
// initiallize acceleration
void ParticleSystem::initiallizeAcc()
{
    for (auto& atom : snapshots[0].atoms)
    {
        atom.acceleration = atom.force / atom.mass;   
        atom.accelerationPrev = atom.force / atom.mass;   
    }
}

// ----------------------------------------
// reset Ekin, Epot, Virial
void ParticleSystem::resetEkinEpotVirial()
{
    Ekin = 0.0, Epot = 0.0, virial = 0.0;
}

// ----------------------------------------
// update pos wth velocity Verlet method
void ParticleSystem::updatePos(const double dt)
{
    for (auto& atom : snapshots[0].atoms)
    {
        atom.position = atom.position + atom.velocity * dt +
            atom.accelerationPrev * 0.5 * (dt * dt);

        // periodic boundary conditions
        applyPBC(atom.position);

        // save prev acceleration
        atom.accelerationPrev = atom.acceleration;
    }
}

// ----------------------------------------
// check if neighboor list needs rebuild
// IDEA: If any particle has moved (distance > rs / 2), we rebuild the neighborList
bool ParticleSystem::needRebuild(const double skin)
{
    // max displacement^2 of any particle
    double maxDisp2 = 0.0;

    for (int i = 0; i < partNumber; i++)
    {
        Vec3 dr = snapshots[0].atoms[i].position - lastBuildPos[i];

        // minimum image
        dr.x -= boxVec.x * std::round(dr.x / boxVec.x);
        dr.y -= boxVec.y * std::round(dr.y / boxVec.y);
        dr.z -= boxVec.z * std::round(dr.z / boxVec.z);

        maxDisp2 = std::max(maxDisp2, dr.dot(dr));
    }

    return maxDisp2 > (skin * 0.5) * (skin * 0.5);
}

// ----------------------------------------
// acceleration update
void ParticleSystem::updateAcc()
{
    for (auto& atom : snapshots[0].atoms)
    {
        atom.acceleration = atom.force / atom.mass; 
    }
}

// ----------------------------------------
// vel and Ekin update
void ParticleSystem::updateVelEkin(const double dt)
{
    for (auto& atom : snapshots[0].atoms)
    {
        atom.velocity = 
            atom.velocity + (atom.accelerationPrev + atom.acceleration) * 0.5 * dt;

        Ekin += 0.5 * atom.mass * atom.velocity.dot(atom.velocity); 
    }
}

// ----------------------------------------
// calc T (kb = 1, lj units)
void ParticleSystem::computeTemperature()
{
    temperature = (2.0 * Ekin) / dof;
}

// ----------------------------------------
// calc pressure (= density * temp + sum(Rij*Fij) / (3V) )
void ParticleSystem::computePressure()
{
    pressure =  density * temperature + virial / (3.0 * volume);
}

// ----------------------------------------
/*
    Thermostat Berendsen
    labda = sqrt(1 + dt/tau * (Ttarget/T - 1))
    v_new = v_old * labda
*/ 
void ParticleSystem::useThermostat(const double Ttarget, const double tau, const double dt)
{
    // avoid division by zero or invalid values
    if (!std::isfinite(temperature) || temperature < 1e-8)
    {
        throw std::runtime_error("Error: Problem with temperature in Thermostat.");
    }

    // scaling factor
    double lambda = std::sqrt(1.0 + (dt / tau) * (Ttarget / temperature - 1.0));

    // change velocity of particles
    for (auto& atom : snapshots[0].atoms)
        atom.velocity = atom.velocity * lambda;
}

// ----------------------------------------
void ParticleSystem::buildNeighborList(const double rlist)
{
    snapshots[0].neighborLists.assign(partNumber, {});

    for (std::size_t i = 0; i < partNumber; ++i)
    {
        for (std::size_t j = 0; j < partNumber; ++j)
        {
            if (i == j)
                continue;

            Vec3 dr;

            dr.x =
                snapshots[0].atoms[j].position.x -
                snapshots[0].atoms[i].position.x;

            dr.y =
                snapshots[0].atoms[j].position.y -
                snapshots[0].atoms[i].position.y;

            dr.z =
                snapshots[0].atoms[j].position.z -
                snapshots[0].atoms[i].position.z;

            // Minimum image convention
            dr.x -= snapshots[0].box.x *
                    std::round(dr.x / snapshots[0].box.x);

            dr.y -= snapshots[0].box.y *
                    std::round(dr.y / snapshots[0].box.y);

            dr.z -= snapshots[0].box.z *
                    std::round(dr.z / snapshots[0].box.z);

            // Distance
            double distance =
                std::sqrt(
                    dr.x * dr.x +
                    dr.y * dr.y +
                    dr.z * dr.z
                );


            // Cutoff rlist = rc + rSkin
            if (distance <= rlist)
            {
                snapshots[0].neighborLists[i].push_back(snapshots[0].atoms[j]);
            }
        }
    }
}

// ----------------------------------------
void ParticleSystem::calcVirial()
{
    virial = 0.0;

    for (int i = 0; i < partNumber; ++i)
    {
        for (auto& neighbor : snapshots[0].neighborLists[i])
        {
            Vec3 dr = snapshots[0].atoms[i].position - neighbor.position;

            // Minimum image convention
            dr.x -= snapshots[0].box.x * std::round(dr.x / snapshots[0].box.x);
            dr.y -= snapshots[0].box.y * std::round(dr.y / snapshots[0].box.y);
            dr.z -= snapshots[0].box.z * std::round(dr.z / snapshots[0].box.z);

            // dr_ij * Fij
            virial += dr.dot(snapshots[0].atoms[i].force);
        }
    }
}