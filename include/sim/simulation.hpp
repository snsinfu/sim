// Simulation drivers

// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef INCLUDED_SIM_SIMULATION_HPP
#define INCLUDED_SIM_SIMULATION_HPP

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
}

#endif
