// C++
#include <iostream>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
// SMASH
#include "smash/config.h"
#include "smash/experiment.h"
#include "smash/forwarddeclarations.h"
#include "smash/input_keys.h"
#include "smash/library.h"
#include "smash/setup_particles_decaymodes.h"
#include "smash/random.h"
// Project
#include "BlastWaveModus.h"

int main(int argc, char const *argv[])
{
    try
    {
        if (argc < 3)
            std::cout << "Usage: ./ThermalBlastMC </path/to/config/file> </path/to/output/directory>";
        std::string configFile = argv[1];
        std::string outputDir = argv[2];
        const std::filesystem::path configPath(configFile);
        if (!std::filesystem::exists(configPath))
            throw std::runtime_error("Invalid config file path!");

        const std::filesystem::path outputPath(outputDir);
        const std::string tabulationsPath("./tabulations");
        const std::string particlesFile(SMASH_INPUT_DIR "/particles.txt");
        const std::string decaymodesFile(SMASH_INPUT_DIR "/decaymodes.txt");

        // Ensure output path exists
        std::filesystem::create_directories(outputPath);

        // Setup SMASH
        smash::Configuration config(configPath.parent_path(), configPath.filename());
        int64_t seed = config.read(smash::InputKeys::gen_randomseed);
        if (seed < 0)
            config.set_value(smash::InputKeys::gen_randomseed, smash::random::generate_63bit_seed());
        const auto pd = smash::load_particles_and_decaymodes(particlesFile, decaymodesFile);
        config.set_value(smash::InputKeys::particles, pd.first);
        config.set_value(smash::InputKeys::decaymodes, pd.second);
        std::string smashVersion = SMASH_VERSION;
        smash::initialize_particles_decays_and_tabulations(config, smashVersion, tabulationsPath);
        smash::set_default_loglevel(config.take(smash::InputKeys::log_default));
        auto loggerConfig = config.extract_sub_configuration(smash::InputSections::logging, smash::Configuration::GetEmpty::Yes);
        if (!loggerConfig.is_empty())
            loggerConfig.enclose_into_section(smash::InputSections::logging);
        smash::create_all_loggers(std::move(loggerConfig));

        // Copy config to output directory
        std::ofstream(outputPath / "config.yaml")
            << "# " << SMASH_VERSION << '\n'
            << "# System   : " << CMAKE_SYSTEM << '\n'
            << "# Compiler : " << CMAKE_CXX_COMPILER_ID << ' '
            << CMAKE_CXX_COMPILER_VERSION << '\n'
            << "# Build    : " << CMAKE_BUILD_TYPE << '\n'
            << "# Date     : " << BUILD_DATE << '\n'
            << config.to_string() << '\n';

        // Create experiment
        auto experiment = smash::Experiment<BlastWaveModus>(config, outputPath);
        // Run!!!
        experiment.run();
    }
    catch (std::exception &e)
    {
        std::cout << "Caught " << e.what() << std::endl;
        std::cout << "Type " << typeid(e).name() << std::endl;
    }
    catch (...)
    {
        std::cout << "Caught exeption!" << std::endl;
    }
    return 0;
}
