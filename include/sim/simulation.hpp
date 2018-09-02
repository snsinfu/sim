// Simulation drivers

// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef INCLUDED_SIM_SIMULATION_HPP
#define INCLUDED_SIM_SIMULATION_HPP

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "system.hpp"


namespace sim
{
    struct newtonian_dynamics_config
    {
        sim::scalar timestep = 1;
        sim::step simulation_length = 1;
    };

    inline void simulate_newtonian_dynamics(
        sim::system& system,
        sim::newtonian_dynamics_config const& config
    )
    {
        auto const particle_count = system.particle_count();
        auto const positions = system.position_array();
        auto const velocities = system.velocity_array();
        auto const masses = system.mass_array();

        std::vector<sim::vector> forces(particle_count);

        system.compute_force(forces);

        for (sim::step stp = 0; stp < config.simulation_length; stp++) {
            for (sim::index i = 0; i < particle_count; i++) {
                velocities[i] += (config.timestep / 2) / masses[i] * forces[i];
                positions[i] += config.timestep * velocities[i];
            }

            system.compute_force(forces);

            for (sim::index i = 0; i < particle_count; i++) {
                velocities[i] += (config.timestep / 2) / masses[i] * forces[i];
            }
        }
    }

    namespace detail
    {
        inline sim::scalar compute_brownian_timestep(
            sim::scalar displacement,
            sim::scalar force,
            sim::scalar mobility,
            sim::scalar temperature
        )
        {
            if (force == 0) {
                return displacement * displacement * M_PI / (16 * mobility * temperature);
            }

            auto const alpha = 2.535;
            auto const fluctuation = alpha * temperature / force;
            auto const drift = std::hypot(fluctuation, displacement) - fluctuation;

            return drift / (mobility * force);
        }
    }

    struct brownian_dynamics_config
    {
        sim::scalar timestep = 1;
        sim::scalar spacestep = 0;
        sim::scalar temperature = 1;
        sim::step simulation_length = 1;
        std::uint32_t random_seed = 0;
    };

    inline void simulate_brownian_dynamics(
        sim::system& system,
        sim::brownian_dynamics_config const& config
    )
    {
        auto const particle_count = system.particle_count();
        auto const positions = system.position_array();
        auto const mobilities = system.mobility_array();

        std::normal_distribution<sim::scalar> normal;
        std::mt19937 random_engine{config.random_seed};

        random_engine.discard(1000000);

        std::vector<sim::vector> forces(particle_count);
        std::vector<sim::vector> previous_weiner(particle_count);

        auto const generate_weiner = [&](sim::scalar mu, sim::scalar dt) {
            return std::sqrt(2 * config.temperature * mu * dt) * sim::vector {
                normal(random_engine),
                normal(random_engine),
                normal(random_engine)
            };
        };

        system.compute_force(forces);

        for (sim::index i = 0; i < particle_count; i++) {
            previous_weiner[i] = generate_weiner(mobilities[i], config.timestep);
        }

        for (sim::step stp = 0; stp < config.simulation_length; stp++) {
            system.compute_force(forces);

            sim::scalar timestep = config.timestep;

            if (config.spacestep > 0) {
                for (sim::index i = 0; i < particle_count; i++) {
                    auto const dt = detail::compute_brownian_timestep(
                        config.spacestep,
                        forces[i].norm(),
                        mobilities[i],
                        config.temperature
                    );
                    if (dt < timestep) {
                        timestep = dt;
                    }
                }
            }

            for (sim::index i = 0; i < particle_count; i++) {
                auto const weiner = generate_weiner(mobilities[i], timestep);
                auto const mean_weiner = 0.5 * (weiner + previous_weiner[i]);
                positions[i] += mobilities[i] * forces[i] * timestep + mean_weiner;
                previous_weiner[i] = weiner;
            }
        }
    }
}

#endif
