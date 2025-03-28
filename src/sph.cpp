#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "SPHLibrary/autopassph.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"

using Particle = sphLib::SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

/**
 * The role of the main function is to instantiate an object of the Simulation
 * class which is actually responsible for the simulation and to run tests for
 * all classes.
 */
int main(int argc, char** argv) {}