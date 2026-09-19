#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <filesystem>
#include <sys/stat.h>
#include "autopas/AutoPas.h"
#include "SPHParticle.h"
// #include "src/TypeDefinitions.h" // Ensure this contains your ParticleType definition

namespace fs = std::filesystem;
using ParticleType = SPHParticle;

class SimpleVtkWriter {
public:
    SimpleVtkWriter(std::string sessionName, std::string outputFolder, int maxDigits, size_t num_probes,
                    const std::string &configPath = "")
        : _sessionName(std::move(sessionName)), 
          _outputFolder(std::move(outputFolder)), 
          _maxDigits(maxDigits) {

        _sessionFolderPath = _outputFolder + "/" + _sessionName + "/";
        _dataFolderPath = _sessionFolderPath + "data/";

        tryCreateFolder(_outputFolder, "./");
        tryCreateFolder(_sessionName, _outputFolder);
        tryCreateFolder("data", _sessionFolderPath);


        _csvPath = _sessionFolderPath + _sessionName + "_probes.csv";
        initializeProbeFile(num_probes);

        if (!configPath.empty() && fs::exists(configPath)) {
            std::string destConfigPath = _sessionFolderPath + _sessionName + ".yaml";
            try {
                fs::copy_file(configPath, destConfigPath, fs::copy_options::overwrite_existing);
            } catch (const fs::filesystem_error &e) {
                AutoPasLog(WARN, "Failed to copy config file to output folder");
            }
        }
    }

    void recordTimestep(size_t currentIteration, 
                        const autopas::AutoPas<ParticleType> &autoPasContainer,
                        const std::array<double, 3>& boxMin,
                        const std::array<double, 3>& boxMax,
                        const autopas::IteratorBehavior iterator = autopas::IteratorBehavior::owned) const {

        recordParticleStates(currentIteration, autoPasContainer, iterator);
        recordDomainSubdivision(currentIteration, boxMin, boxMax);
    }

    void recordProbes(size_t iteration, double simTime, size_t fileStep, double cutoff, double smoothingLength,
                      const autopas::AutoPas<ParticleType> &container,
                      const std::vector<std::array<double, 3>> &probes) {

        std::ofstream csv(_csvPath, std::ios::app);
        if (!csv.is_open()) return;

        size_t id = 1;

        csv << iteration << "," << simTime << "," << fileStep;

        for (const auto& probe : probes) {
            double probeDensity = 0.0;
            double probePressure = 0.0;
            std::array<double, 3> probeVel = {0.0, 0.0, 0.0};
            std::array<double, 3> probeForce = {0.0, 0.0, 0.0};

            probeMeasurement(container, probe, cutoff, smoothingLength, probeVel, probeForce, probeDensity, probePressure);

            csv << "," << probe[0] << "," << probe[1] << "," << probe[2] << ","
                << probeVel[0] << "," << probeVel[1] << "," << probeVel[2] << ","
                << probeForce[0] << "," << probeForce[1] << "," << probeForce[2] << ","
                << probeDensity << "," << probePressure ;
        }

        csv << "\n" ;
    }

private:
    std::string _sessionName;
    std::string _outputFolder;
    std::string _sessionFolderPath;
    std::string _dataFolderPath;
    std::string _csvPath;
    int _maxDigits;

    static void tryCreateFolder(const std::string &name, const std::string &location) {
        std::string path = location + "/" + name;
        mkdir(path.c_str(), 0777); 
    }

    void initializeProbeFile(size_t num_probes) {
        std::ofstream csv(_csvPath, std::ios::out);
        if (csv.is_open()) {
            csv << "iteration,time,file_step" ;
            for (size_t id = 0; id < num_probes; ++id) {
                csv << ",x_" << id << ",y_" << id << ",z_" << id << ",v_x_" << id << ",v_y_" << id << ",v_z_" << id 
                    << ",f_x_" << id << ",f_y_" << id << ",f_z_" << id << ",density_" << id << ",pressure_" << id ;
            }
            csv << "\n" ;
        } else {
            AutoPasLog(WARN, "Failed to initialize probe CSV file: {}", _csvPath);
        }
    }


    void generateFilename(const std::string &tag, size_t iteration, std::ostringstream &stream) const {
        stream << _dataFolderPath << _sessionName << "_" << tag << "_0_" 
               << std::setfill('0') << std::setw(_maxDigits) << iteration << ".vtu";
    }

    void createPvtuHeader(const std::string& tag, size_t iteration, bool isParticle) const {
        std::ostringstream filename;
        filename << _sessionFolderPath << _sessionName << "_" << tag << "_" 
                 << std::setfill('0') << std::setw(_maxDigits) << iteration << ".pvtu";

        std::ofstream file(filename.str(), std::ios::out | std::ios::binary);
        if (!file.is_open()) return;

        file << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n"
             << "<VTKFile byte_order=\"LittleEndian\" type=\"PUnstructuredGrid\" version=\"0.1\">\n"
             << "  <PUnstructuredGrid GhostLevel=\"0\">\n";

        if (isParticle) {
            file << "    <PPointData>\n"
                 << "      <PDataArray Name=\"velocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n"
                 << "      <PDataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/>\n"
                 << "      <PDataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\"/>\n"
                 << "      <PDataArray Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>\n"
                 << "      <PDataArray Name=\"mass\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>\n"
                 << "      <PDataArray Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\"/>\n"
                 << "    </PPointData>\n"
                 << "    <PPoints><PDataArray Name=\"positions\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"/></PPoints>\n";
        } else {
            file << "    <PCellData><PDataArray type=\"Int32\" Name=\"Rank\" /></PCellData>\n"
                 << "    <PPoints><DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\"></DataArray></PPoints>\n";
        }

        file << "    <Piece Source=\"./data/" << _sessionName << "_" << tag << "_0_" 
             << std::setfill('0') << std::setw(_maxDigits) << iteration << ".vtu\"/>\n"
             << "  </PUnstructuredGrid>\n</VTKFile>\n";
    }

    void recordParticleStates(size_t iteration, const autopas::AutoPas<ParticleType> &container, autopas::IteratorBehavior iterator) const {
        createPvtuHeader("Particles", iteration, true);

        std::ostringstream filename;
        generateFilename("Particles", iteration, filename);

        std::ofstream file(filename.str(), std::ios::out | std::ios::binary);
        if (!file.is_open()) throw std::runtime_error("Failed to open particle file");

        size_t numParticles = container.getNumberOfParticles(iterator);

        file << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n"
             << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">\n"
             << "  <UnstructuredGrid>\n"
             << "    <Piece NumberOfCells=\"0\" NumberOfPoints=\"" << numParticles << "\">\n"
             << "      <PointData>\n";

        // Velocities
        file << "        <DataArray Name=\"velocities\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
        for (auto p = container.begin(iterator); p.isValid(); ++p) {
            auto v = p->getV(); file << "        " << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        file << "        </DataArray>\n";

        // Forces
        file << "        <DataArray Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
        for (auto p = container.begin(iterator); p.isValid(); ++p) {
            auto f = p->getF(); file << "        " << f[0] << " " << f[1] << " " << f[2] << "\n";
        }
        file << "        </DataArray>\n";

        file << "        <DataArray Name=\"ids\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
        for (auto p = container.begin(iterator); p.isValid(); ++p) { file << "        " << p->getID() << "\n"; }
        file << "        </DataArray>\n";

        file << "        <DataArray Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">\n";
        for (auto p = container.begin(iterator); p.isValid(); ++p) { file << "        " << p->getDensity() << "\n"; }
        file << "        </DataArray>\n";

        file << "        <DataArray Name=\"mass\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">\n";
        for (auto p = container.begin(iterator); p.isValid(); ++p) { file << "        " << p->getMass() << "\n"; }
        file << "        </DataArray>\n";

        file << "        <DataArray Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">\n";
        for (auto p = container.begin(iterator); p.isValid(); ++p) { file << "        " << p->getPressure() + 100000.0 << "\n"; }
        file << "        </DataArray>\n";

        file << "      </PointData>\n<CellData/>\n<Points>\n"
             << "        <DataArray Name=\"positions\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";

        // Positions with boundary dynamic precision logic kept intact
        const auto boxMax = container.getBoxMax();
        for (auto p = container.begin(iterator); p.isValid(); ++p) {
            auto pos = p->getR();
            file << "        ";
            for(int d=0; d<3; ++d) {
                if (boxMax[d] - pos[d] < 0.1) {
                    // Quick fallback to high-precision formatting if near border
                    file << std::setprecision(16) << pos[d] << std::setprecision(6) << " ";
                } else {
                    file << pos[d] << " ";
                }
            }
            file << "\n";
        }

        file << "        </DataArray>\n      </Points>\n"
             << "      <Cells><DataArray Name=\"types\" NumberOfComponents=\"0\" format=\"ascii\" type=\"Float32\"/></Cells>\n"
             << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    }

    void recordDomainSubdivision(size_t iteration, const std::array<double, 3>& boxMin, const std::array<double, 3>& boxMax) const {
        createPvtuHeader("Ranks", iteration, false);

        std::ostringstream filename;
        generateFilename("Ranks", iteration, filename);

        std::ofstream file(filename.str(), std::ios::out | std::ios::binary);
        if (!file.is_open()) return;

        file << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n"
             << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">\n"
             << "  <UnstructuredGrid>\n"
             << "    <Piece NumberOfPoints=\"8\" NumberOfCells=\"1\">\n"
             << "      <CellData>\n"
             << "        <DataArray type=\"Int32\" Name=\"Rank\" format=\"ascii\">\n          0\n        </DataArray>\n"
             << "      </CellData>\n"
             << "      <Points>\n"
             << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"
             << "          " << boxMin[0] << " " << boxMin[1] << " " << boxMin[2] << "\n"
             << "          " << boxMin[0] << " " << boxMin[1] << " " << boxMax[2] << "\n"
             << "          " << boxMin[0] << " " << boxMax[1] << " " << boxMin[2] << "\n"
             << "          " << boxMin[0] << " " << boxMax[1] << " " << boxMax[2] << "\n"
             << "          " << boxMax[0] << " " << boxMin[1] << " " << boxMin[2] << "\n"
             << "          " << boxMax[0] << " " << boxMin[1] << " " << boxMax[2] << "\n"
             << "          " << boxMax[0] << " " << boxMax[1] << " " << boxMin[2] << "\n"
             << "          " << boxMax[0] << " " << boxMax[1] << " " << boxMax[2] << "\n"
             << "        </DataArray>\n"
             << "      </Points>\n"
             << "      <Cells>\n"
             << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n          0 1 2 3 4 5 6 7\n        </DataArray>\n"
             << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          8\n        </DataArray>\n"
             << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          11\n        </DataArray>\n"
             << "      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    }

    void probeMeasurement(const autopas::AutoPas<ParticleType> &container, std::array<double, 3> pos,
                          const double cutoff, const double smoothingLength, std::array<double, 3> &vel,
                          std::array<double, 3> &force, double &density, double &pressure) {
        using namespace autopas::utils::ArrayMath::literals;

        const std::array<double, 3> lowCorner = pos - cutoff;
        const std::array<double, 3> highCorner = pos + cutoff;
        for(auto part = container.getRegionIterator(lowCorner, highCorner); part.isValid(); ++part) {
            const std::array<double, 3> dr = pos - part->getR();
            const double W = SPHKernels::W(dr, smoothingLength);
            if (W > 0) {
                const double scale = W * part->getMass() / part->getDensity();
                vel += part->getV() * scale;
                force += part->getAcceleration() * scale;
                density += scale * part->getDensity();
                pressure += scale * part->getPressure();
            }
        }
    }
};
