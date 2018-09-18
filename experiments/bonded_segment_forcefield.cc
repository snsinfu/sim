#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "../include/sim/forcefield.hpp"
#include "../include/sim/potential.hpp"
#include "../include/sim/system.hpp"
#include "../include/sim/typedef.hpp"


namespace
{
    template<typename Derived>
    class bonded_segment_forcefield : public sim::forcefield
    {
    public:
        // bonded_segments contain inclusive index ranges of sequentially-bonded particles.
        std::vector<std::pair<sim::index, sim::index>> bonded_segments;

        // compute_energy computes the total bond potential energy of the system.
        sim::scalar compute_energy(sim::system const& system) override
        {
            sim::scalar sum = 0;

            auto& derived = static_cast<Derived&>(*this);

            foreach_pair(system, [&](sim::vector r_ij, sim::index i, sim::index j) {
                auto potential = derived.bonded_segment_potential(system, i, i + 1);
                sum += potential.evaluate_energy(r_ij);
            });

            return sum;
        }

        // compute_force computes bond force acting on each particle and accumulates the results
        // into forces_array.
        void compute_force(sim::system const& system, sim::array_view<sim::vector> force_array) override
        {
            auto& derived = static_cast<Derived&>(*this);

            foreach_pair(system, [&](sim::vector r_ij, sim::index i, sim::index j) {
                auto potential = derived.bonded_segment_potential(system, i, i + 1);
                auto force = potential.evaluate_force(r_ij);
                force_array[i] += force;
                force_array[i + 1] -= force;
            });
        }

    private:
        // foreach_pair iterates over all bonded pairs and invokes func for each pair with the
        // displacement between particle positions and the indices of the particles.
        template<typename Func>
        inline void foreach_pair(sim::system const& system, Func func) const
        {
            auto const position_array = system.position_array();

            for (auto const& segment : bonded_segments) {
                for (sim::index i = segment.first; i < segment.second; i++) {
                    auto const r_ij = position_array[i] - position_array[i + 1];
                    func(r_ij, i, i + 1);
                }
            }
        }
    };

    template<typename PotentialFunc>
    void force_bonded_segment(
        sim::system& system,
        std::vector<std::pair<sim::index, sim::index>> const& segments,
        PotentialFunc potential
    )
    {
        using Potential = decltype(potential(sim::index(0), sim::index(0)));

        class impl : public bonded_segment_forcefield<impl>
        {
        public:
            explicit impl(PotentialFunc potential)
                : potential_{potential}
            {
            }

            Potential bonded_segment_potential(sim::system const&, sim::index i, sim::index j) const
            {
                return potential_(i, j);
            }

        private:
            PotentialFunc potential_;
        };

        auto forcefield = std::make_shared<impl>(potential);
        forcefield->bonded_segments = segments;
        system.add_forcefield(forcefield);
    }
}


int main()
{
    sim::system system;

    for (int i = 0; i < 100; i++) {
        system.add_particle({
            .position = {sim::scalar(i), 0, 0}
        });
    }

    force_bonded_segment(
        system,
        {{0, 49}, {50, 99}},
        [](sim::index i, sim::index j) {
            auto sigma_i = (i % 2 ? 0.1 : 0.2);
            auto sigma_j = (j % 2 ? 0.1 : 0.2);
            auto sigma = std::sqrt(sigma_i * sigma_i + sigma_j * sigma_j);
            return sim::harmonic_potential{.spring_constant = 10 / (sigma * sigma)};
        }
    );

    std::cout << system.compute_potential_energy() << '\n';
}
