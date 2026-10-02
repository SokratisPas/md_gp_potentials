#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <random>

#include "vec3.h"



// =======================================
// Particle Class
struct Particle {
    Vec3 position = Vec3(0.0);       // possition
    Vec3 velocity = Vec3(0.0);       // velocity
    Vec3 acceleration = Vec3(0.0);   // acceleration

    Vec3 force = Vec3(0.0);          // force

    Vec3 accelerationPrev = Vec3(0.0);   // previous acceleration

    double mass = 1.0;                          // mass (lj)
};


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

    double dof = 1.0;           // degrees of Freedom

    // Cube dimensions 
    double xMin = 0.0, xMax = 1.0;
    double yMin = 0.0, yMax = 1.0;
    double zMin = 0.0, zMax = 1.0;

    // vector of particles
    std::vector<Particle> particles;

    // each particle[i] has a vector that includes all the neighbor particles
    std::vector<std::vector<int>> neighborList;

    // last pos of particle[i] 
    std::vector<Vec3> lastBuildPos;

    // Constructor
    ParticleSystem(int in_partNumber)
    {
        partNumber = in_partNumber;

        // degrees of freedom (we remove 3 bcz we remove net momentum)
        dof = 3.0 * static_cast<double>(partNumber) - 3.0;

        particles.resize(partNumber);
        neighborList.resize(partNumber);
        lastBuildPos.resize(partNumber);
    }

    // ----------------------------------------
    // calc random double in [start, finish]
    double randomDouble(double start, double finish)
    {
        static std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(start, finish);
        return dist(gen);
    }

    // ----------------------------------------
    // particles positions in FCC latice
    std::vector<std::array<double, 3>> generateFCC(const int n, const double a)
    {
        std::vector<std::array<double, 3>> positions;

        // FCC has 4 atoms per unit cell
        int cells_per_dim = std::ceil(std::cbrt(n / 4.0));

        // Basis positions inside one FCC cell
        std::vector<std::array<double, 3>> basis = {
            {0.0, 0.0, 0.0},
            {0.0, 0.5, 0.5},
            {0.5, 0.0, 0.5},
            {0.5, 0.5, 0.0}
        };

        for (int i = 0; i < cells_per_dim; i++) {
            for (int j = 0; j < cells_per_dim; j++) {
                for (int k = 0; k < cells_per_dim; k++) {

                    for (const auto& b : basis) {
                        if (positions.size() >= n) return positions;

                        positions.push_back({
                            a * (i + b[0]),
                            a * (j + b[1]),
                            a * (k + b[2])
                            });
                    }
                }
            }
        }
        return positions;
    }

    // ---------------------------------------------------
    // initialize Pos, Vel, lastBuildPos (NOT Acceleration)
    void initializeParameteres(const double laticeConst) 
    {
        // initial pos vector
        std::vector<std::array<double, 3>> initialPos = 
            generateFCC(partNumber, laticeConst); 

        for (int i = 0; i < partNumber; i++)
        {
            // initial pos
            particles[i].position = Vec3(initialPos[i][0], initialPos[i][1], initialPos[i][2]);

            // initial small random vel
            particles[i].velocity = Vec3(randomDouble(-0.1, 0.1), 
                                            randomDouble(-0.1, 0.1),
                                            randomDouble(-0.1, 0.1));

            // initial lastBuildPos
            lastBuildPos[i] = particles[i].position;
        }
    }

    // ----------------------------------------
    // remove net momentum
    // it should : sum(vel[i] == 0)
    void removeNetMomentum()
    {
        Vec3 velMean(0.0);
        for (auto& p : particles)
            velMean = velMean + p.velocity;

        velMean = velMean / static_cast<double>(particles.size());

        for (auto& p : particles)
            p.velocity = p.velocity - velMean;
    }
      
    // ----------------------------------------
    // aplly PBC 
    void applyPBC(Vec3& pos)
    {
        // cube lengths
        double Lx = xMax - xMin;
        double Ly = yMax - yMin;
        double Lz = zMax - zMin;

        // keep x, y, x in cube boundaries 
        pos.x = xMin + std::fmod(std::fmod(pos.x - xMin, Lx) + Lx, Lx);
        pos.y = yMin + std::fmod(std::fmod(pos.y - yMin, Ly) + Ly, Ly);
        pos.z = zMin + std::fmod(std::fmod(pos.z - zMin, Lz) + Lz, Lz);
    }

    // ----------------------------------------
    // Build Neighbour List
    void buildNeighborList(const double rlist)
    {
        double rlist2 = rlist * rlist;

        int N = partNumber;

        // clear neighborList
        for (auto& nl : neighborList)
            nl.clear();

        // build cells (using rlist as cell size)
        double cellSize = rlist;

        // number of cels
        int nx = static_cast<int>((xMax - xMin) / cellSize);
        int ny = nx;
        int nz = nx;

        // debug line for n <= 0
        if (nx <= 0 || ny <= 0 || nz <= 0)
        {
            throw std::runtime_error("Error: Number of cells some dim is non-positive.");
        }

        std::vector<std::vector<int>> cells(nx * ny * nz);

        // make 1D coords to 3D coords
        auto getCellIndex = [&](int ix, int iy, int iz)
            {
                return ix + nx * (iy + ny * iz);
            };

        // assign particles to cells
        for (int i = 0; i < N; i++)
        {
            Vec3 p = particles[i].position;

            // cells indices
            int ix = static_cast<int>((p.x - xMin) / cellSize) % nx;
            int iy = static_cast<int>((p.y - yMin) / cellSize) % ny;
            int iz = static_cast<int>((p.z - zMin) / cellSize) % nz;

            // catch negative indices
            if (ix < 0) ix += nx;
            if (iy < 0) iy += ny;
            if (iz < 0) iz += nz;

            cells[getCellIndex(ix, iy, iz)].push_back(i);
        }

        // build neighbor list
        for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
                for (int iz = 0; iz < nz; iz++)
                {
                    int cellA = getCellIndex(ix, iy, iz);

                    // check only neighbors of cellA
                    for (int dx = -1; dx <= 1; dx++)
                        for (int dy = -1; dy <= 1; dy++)
                            for (int dz = -1; dz <= 1; dz++)
                            {
                                int jx = (ix + dx + nx) % nx;
                                int jy = (iy + dy + ny) % ny;
                                int jz = (iz + dz + nz) % nz;

                                int cellB = getCellIndex(jx, jy, jz);

                                for (int i : cells[cellA])
                                {
                                    for (int j : cells[cellB])
                                    {
                                        // avoid double counting
                                        if (i >= j) continue;

                                        Vec3 r = particles[i].position - particles[j].position;

                                        // minimum image
                                        double Lx = xMax - xMin;
                                        double Ly = yMax - yMin;
                                        double Lz = zMax - zMin;

                                        r.x -= Lx * std::round(r.x / Lx);
                                        r.y -= Ly * std::round(r.y / Ly);
                                        r.z -= Lz * std::round(r.z / Lz);

                                        double dist2 = r.dot(r); // r * r
 
                                        if (dist2 < rlist2)
                                        {
                                            neighborList[i].push_back(j);
                                            neighborList[j].push_back(i);
                                        }
                                    }
                                }
                            }
                }

        // store positions
        for (int i = 0; i < N; i++)
            lastBuildPos[i] = particles[i].position;
    }

    // ----------------------------------------
    // particle-particle lj force
    Vec3 computeLJForce(
        const Vec3& pi, 
        const Vec3& pj, 
        const double epsilon, 
        const double sigma)
    {
        Vec3 r = pi - pj; // distance vec of pair

        // box dimensions
        double Lx = xMax - xMin;
        double Ly = yMax - yMin;
        double Lz = zMax - zMin;

        // minimum image convention
        // IDEA:
        // use the smallest distance between 2 particles 
        // r / L callc how many boxes is the distance 
        // and we remove this distance from r
        r.x -= Lx * std::round(r.x / Lx);
        r.y -= Ly * std::round(r.y / Ly);
        r.z -= Lz * std::round(r.z / Lz);

        double dist2 = r.dot(r);

        // avoid division by zero 
        if (dist2 < 1e-12)
            return Vec3(0.0);
            
        // cut off radius
        double rc = 2.5 * sigma;
        if (dist2 > (rc * rc))
            return Vec3(0.0);

        // distance magnitude
        double dist = std::sqrt(dist2);

        // inverse distance
        double inv_r = 1.0 / std::sqrt(dist2);

        double s_over_r = sigma * inv_r;
        double s2 = s_over_r * s_over_r;
        double s6 = s2 * s2 * s2;
        double s12 = s6 * s6;

        // E potential update
        Epot += 4.0 * epsilon * (s12 - s6);

        // LJ force magnitude
        double f = 24.0 * epsilon * (2.0 * s12 - s6) * inv_r;

        return r * (f * inv_r);
    }

    // ----------------------------------------
    // calc all lj forces and Epot and Virial 
    void computeForces_Epot_Virial(const double epsilon, const double sigma)
    {
        int N = partNumber;

        // set forces = 0
        for (auto& p : particles)
            p.force = Vec3(0.0);

        // calc new lj forces
        for (int i = 0; i < N; i++)
        {
            for (int j : neighborList[i])
            {
                if (i >= j) continue;

                Vec3 ijForce = computeLJForce(particles[i].position, particles[j].position, epsilon, sigma);

                particles[i].force = particles[i].force + ijForce;
                particles[j].force = particles[j].force - ijForce;


                Vec3 rij = particles[i].position - particles[j].position;
                virial += rij.dot(ijForce);
            }
        }
    }

    // ----------------------------------------
    // calc volume, density
    void calcVolumeDensity()
    {
        double volume = (xMax - xMin) * (yMax - yMin) * (zMax - zMin);
        double density = partNumber / volume;
    }

    // ----------------------------------------
    // initiallize acceleration
    void initiallizeAcc()
    {
        for (auto& p : particles)
        {
            p.acceleration = p.force;       // m = 1
            p.accelerationPrev = p.force;   // m = 1
        }
    }

    // ----------------------------------------
    // reset Ekin, Epot, Virial
    void resetEkinEpotVirial()
    {
        Ekin = 0.0, Epot = 0.0, virial = 0.0;
    }

    // ----------------------------------------
    // update pos wth velocity Verlet method
    void updatePos(const double dt)
    {
        for (auto& part : particles)
        {
            part.position = part.position + part.velocity * dt +
                part.accelerationPrev * 0.5 * (dt * dt);

            // periodic boundary conditions
            applyPBC(part.position);

            // save prev acceleration
            part.accelerationPrev = part.acceleration;
        }
    }

    // ----------------------------------------
    // check if neighboor list needs rebuild
    // IDEA: If any particle has moved (distance > rs / 2), we rebuild the neighborList
    bool needRebuild(const double skin)
    {
        int N = partNumber;

        // max displacement^2 of any particle
        double maxDisp2 = 0.0;

        for (int i = 0; i < N; i++)
        {
            Vec3 dr = particles[i].position - lastBuildPos[i];

            // minimum image
            double Lx = xMax - xMin;
            double Ly = yMax - yMin;
            double Lz = zMax - zMin;

            dr.x -= Lx * std::round(dr.x / Lx);
            dr.y -= Ly * std::round(dr.y / Ly);
            dr.z -= Lz * std::round(dr.z / Lz);

            maxDisp2 = std::max(maxDisp2, dr.dot(dr));
        }

        return maxDisp2 > (skin * 0.5) * (skin * 0.5);
    }

    // ----------------------------------------
    // acceleration update
    void updateAcc()
    {
        for (auto& part : particles)
        {
            part.acceleration = part.force; // m = 1
        }
    }

    // ----------------------------------------
    // vel and Ekin update
    void updateVelEkin(const double dt)
    {
        for (auto& part : particles)
        {
            part.velocity = 
                part.velocity + (part.accelerationPrev + part.acceleration) * 0.5 * dt;

            Ekin += 0.5 * part.velocity.dot(part.velocity); // m = 1
        }
    }

    // ----------------------------------------
    // calc T (kb = 1, lj units)
    void computeTemperature()
    {
        temperature = (2.0 * Ekin) / dof;
    }

    // ----------------------------------------
    // calc pressure (= density * temp + sum(Rij*Fij) / (3V) )
    void computePressure()
    {
        pressure =  density * temperature + virial / (3.0 * volume);
    }

    // ----------------------------------------
    /*
        Thermostat Berendsen
        labda = sqrt(1 + dt/tau * (Ttarget/T - 1))
        v_new = v_old * labda
    */ 
    void useThermostat(const double Ttarget, const double tau, const double dt)
    {
        // avoid division by zero or invalid values
        if (!std::isfinite(temperature) || temperature < 1e-8)
        {
            throw std::runtime_error("Error: Problem with temperature in Thermostat.");
        }

        // scaling factor
        double lambda = std::sqrt(1.0 + (dt / tau) * (Ttarget / temperature - 1.0));

        // change velocity of particles
        for (auto& p : particles)
            p.velocity = p.velocity * lambda;
    }
};