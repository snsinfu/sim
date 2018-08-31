// Customizable forcefield implementations

// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef INCLUDED_SIM_FORCEFIELD_TEMPLATE_HPP
#define INCLUDED_SIM_FORCEFIELD_TEMPLATE_HPP

#include "forcefield.hpp"
#include "system.hpp"
#include "vendor/neighbor_searcher.hpp"


namespace sim
{
    template<typename Derived>
    class pair_forcefield : public sim::forcefield
    {
    public:
        sim::scalar compute_energy(sim::system const& system) override
        {
            sim::scalar energy = 0;

            auto const particle_count = system.particle_count();
            auto const positions = system.position_array();

            for (sim::index j = 0; j < particle_count; j++) {
                auto const position_j = positions[j];

                for (sim::index i = 0; i < j; i++) {
                    auto r_ij = positions[i] - position_j;
                    auto potential = derived().pair_potential(system, i, j);

                    energy += potential.evaluate_energy(r_ij);
                }
            }

            return energy;
        }

        void compute_force(sim::system const& system, sim::array_view<sim::vector> forces) override
        {
            auto const particle_count = system.particle_count();
            auto const positions = system.position_array();

            for (sim::index j = 0; j < particle_count; j++) {
                auto const position_j = positions[j];
                sim::vector reaction;

                for (sim::index i = 0; i < j; i++) {
                    auto r_ij = positions[i] - position_j;
                    auto potential = derived().pair_potential(system, i, j);
                    auto force = potential.evaluate_force(r_ij);

                    forces[i] += force;
                    reaction -= force;
                }
                forces[j] += reaction;
            }
        }

    private:
        inline Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
    };
}

#endif
