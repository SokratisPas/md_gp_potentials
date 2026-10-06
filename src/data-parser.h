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
struct Atom
{
    std::string type = "None";

    Vec3 position = Vec3(0.0);    
    Vec3 velocity = Vec3(0.0);     
    Vec3 acceleration = Vec3(0.0);   
    Vec3 force = Vec3(0.0);       

    Vec3 accelerationPrev = Vec3(0.0);   // previous acceleration

    double mass = 1.0;                   // mass (lj)
};

// =============================================================
struct Snapshot
{
    int numberOfAtoms;
    std::vector<Atom> atoms;
    std::vector<std::vector<Atom>> neighborLists; // vector of neighbors for each atom

    Vec3 box;
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

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(file, line))
    {
        if (!line.empty())
            lines.push_back(line);
    }

    std::vector<Snapshot> snapshots;

    std::size_t i = 0;
    while (i < lines.size())
    {
        const std::string& currentLine = lines[i];

        std::string metadataLine = currentLine;
        int numberOfAtoms = 0;

        std::stringstream firstTokenStream(currentLine);
        std::string firstToken;
        if (firstTokenStream >> firstToken &&
            std::all_of(firstToken.begin(), firstToken.end(), [](unsigned char ch)
            {
                return std::isdigit(ch);
            }))
        {
            numberOfAtoms = std::stoi(firstToken);
            ++i;

            if (i >= lines.size())
                throw std::runtime_error(
                    "Unexpected end of file while reading metadata."
                );

            metadataLine = lines[i];
        }

        if (metadataLine.find("Lattice=") == std::string::npos &&
            metadataLine.find("Properties=") == std::string::npos)
        {
            throw std::runtime_error(
                "Unexpected metadata line: " + metadataLine
            );
        }

        Snapshot snapshot;
        snapshot.numberOfAtoms = numberOfAtoms;
        snapshot.box = parseLattice(metadataLine);
        snapshot.configType = parseConfigType(metadataLine);
        snapshot.energy = parseEnergy(metadataLine);

        std::size_t atomStart = i + 1;
        std::size_t atomEnd = atomStart;

        if (numberOfAtoms > 0)
        {
            atomEnd = std::min(atomStart + static_cast<std::size_t>(numberOfAtoms), lines.size());
        }
        else
        {
            while (atomEnd < lines.size())
            {
                const std::string& candidate = lines[atomEnd];
                if (candidate.find("Lattice=") != std::string::npos ||
                    candidate.find("Properties=") != std::string::npos)
                {
                    break;
                }
                ++atomEnd;
            }
        }

        snapshot.atoms.reserve(static_cast<std::size_t>(std::max(numberOfAtoms, static_cast<int>(atomEnd - atomStart))));

        for (std::size_t atomIndex = atomStart; atomIndex < atomEnd; ++atomIndex)
        {
            std::stringstream ss(lines[atomIndex]);

            Atom atom;

            ss >> atom.type;

            ss >> atom.position.x
                >> atom.position.y
                >> atom.position.z;

            Vec3 force;
            ss >> force.x
                >> force.y
                >> force.z;

            int atomicNumber = 0;
            ss >> atomicNumber;

            snapshot.atoms.push_back(atom);
        }

        if (numberOfAtoms > 0 && static_cast<int>(snapshot.atoms.size()) != numberOfAtoms)
        {
            throw std::runtime_error(
                "Expected " + std::to_string(numberOfAtoms) +
                " atoms but found " + std::to_string(snapshot.atoms.size()) +
                " in file: " + filename
            );
        }

        snapshot.numberOfAtoms = static_cast<int>(snapshot.atoms.size());

        buildNeighborLists(snapshot);
        snapshots.push_back(snapshot);

        i = atomEnd;
    }

    if (snapshots.empty())
    {
        throw std::runtime_error("No valid snapshots found in file: " + filename);
    }

    return snapshots;
}

// ------------------------------------------------
Vec3 XYZParser::parseLattice(const std::string& line)
{
    std::size_t start = line.find("Lattice=\"");

    if (start == std::string::npos)
        throw std::runtime_error("Lattice information not found.");

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
    const std::size_t atomCount = snapshot.atoms.size();
    snapshot.neighborLists.assign(atomCount, {});

    for (std::size_t i = 0; i < atomCount; ++i)
    {
        for (std::size_t j = 0; j < atomCount; ++j)
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
                snapshot.neighborLists[i].push_back(snapshot.atoms[j]);
            }
        }
    }
}