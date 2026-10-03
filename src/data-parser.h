#pragma once

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

#include "vec3.h"


// =============================================================
struct Neighbor
{
    std::string type;
    double distance;
};

// =============================================================
struct Atom
{
    std::string type;
    Vec3 position;
    std::vector<Neighbor> neighbors;
};

// =============================================================
struct Snapshot
{
    int numberOfAtoms;
    Vec3 box;
    std::vector<Atom> atoms;
    std::string configType;
    double energy;
};

// =============================================================
class XYZParser
{
public:
    
    explicit XYZParser(const double cutoff)
        : cutoff(cutoff)
    {
    }

    // vector of snapshots
    std::vector<Snapshot> parse(const std::string& filename);  

private:

    const double cutoff;

    Vec3 parseLattice(const std::string& line);
    std::string parseConfigType(const std::string& line);
    double parseEnergy(const std::string& line);
    void buildNeighborLists(Snapshot& snapshot);
};


// ------------------------------------------------
std::vector<Snapshot> XYZParser::parse(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error(
            "Could not open XYZ file: " + filename
        );
    }

    std::vector<Snapshot> snapshots;

    std::string line;

    while (std::getline(file, line))
    {
        // Skip empty lines
        if (line.empty())
            continue;

        // number of atoms
        int numberOfAtoms = std::stoi(line);

        
        // Metadata line (we read only energy, config_type, and lattice(box))
        if (!std::getline(file, line))
            throw std::runtime_error(
                "Unexpected end of file while reading metadata."
            );

        Snapshot snapshot;

        snapshot.numberOfAtoms = numberOfAtoms;
        snapshot.box = parseLattice(line);
        snapshot.configType = parseConfigType(line);
        snapshot.energy = parseEnergy(line);

        snapshot.atoms.reserve(numberOfAtoms);

        // Read atoms
        for (int i = 0; i < numberOfAtoms; ++i)
        {
            if (!std::getline(file, line))
                throw std::runtime_error(
                    "Unexpected end of file while reading atoms."
                );

            std::stringstream ss(line);

            Atom atom;

            ss >> atom.type;

            ss >> atom.position.x
                >> atom.position.y
                >> atom.position.z;

            Vec3 force;     // read force (not used)
            ss >> force.x
                >> force.y
                >> force.z;

            int atomicNumber;   // Read atomic number (not used)
            ss >> atomicNumber; 

            snapshot.atoms.push_back(atom);
        }

        // Calculate neighbors
        buildNeighborLists(snapshot);

        snapshots.push_back(snapshot);
    }

    return snapshots;
}

// ------------------------------------------------
Vec3 XYZParser::parseLattice(const std::string& line)
{
    std::size_t start = line.find("Lattice=\"");

    if (start == std::string::npos)
        throw std::runtime_error(
            "Lattice information not found."
        );

    start += 9;

    std::size_t end = line.find("\"", start);

    std::string lattice =
        line.substr(start, end - start);

    std::stringstream ss(lattice);

    double a1, a2, a3;
    double b1, b2, b3;
    double c1, c2, c3;

    ss >> a1 >> a2 >> a3
        >> b1 >> b2 >> b3
        >> c1 >> c2 >> c3;

    /*
    For these current files the lattice is orthogonal:

            a = (Lx, 0, 0)
            b = (0, Ly, 0)
            c = (0, 0, Lz)

    So we only need the diagonal values.
        */

    return {
        a1,
        b2,
        c3
    };
}

// ------------------------------------------------
std::string XYZParser::parseConfigType(const std::string& line)
{
    std::string key = "config_type=";

    std::size_t start = line.find(key);

    if (start == std::string::npos)
        return "";

    start += key.length();

    std::size_t end = line.find(" ", start);

    if (end == std::string::npos)
        end = line.length();

    return line.substr(start, end - start);
}

// ------------------------------------------------
double XYZParser::parseEnergy(const std::string& line)
{
    std::string key = "energy=";

    std::size_t start = line.find(key);

    if (start == std::string::npos)
        return 0.0;

    start += key.length();

    std::size_t end = line.find(" ", start);

    if (end == std::string::npos)
        end = line.length();

    return std::stod(
        line.substr(start, end - start)
    );
}

// ------------------------------------------------
void XYZParser::buildNeighborLists(Snapshot& snapshot)
{
    const int N = snapshot.numberOfAtoms;

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
                continue;

            Vec3 dr;

            dr.x =
                snapshot.atoms[j].position.x -
                snapshot.atoms[i].position.x;

            dr.y =
                snapshot.atoms[j].position.y -
                snapshot.atoms[i].position.y;

            dr.z =
                snapshot.atoms[j].position.z -
                snapshot.atoms[i].position.z;

            // Minimum image convention
            dr.x -= snapshot.box.x *
                    std::round(dr.x / snapshot.box.x);

            dr.y -= snapshot.box.y *
                    std::round(dr.y / snapshot.box.y);

            dr.z -= snapshot.box.z *
                    std::round(dr.z / snapshot.box.z);

            // Distance
            double distance =
                std::sqrt(
                    dr.x * dr.x +
                    dr.y * dr.y +
                    dr.z * dr.z
                );


            // Cutoff
            if (distance <= cutoff)
            {
                Neighbor neighbor;

                neighbor.type =
                    snapshot.atoms[j].type;

                neighbor.distance = distance;

                snapshot.atoms[i]
                    .neighbors
                    .push_back(neighbor);
            }
        }
    }
}