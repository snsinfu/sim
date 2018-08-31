#include <iostream>
#include <memory>
#include <random>

#include <sim/all.hpp>


class my_forcefield : public sim::pair_forcefield<my_forcefield>
{
public:
    inline auto pair_potential(sim::system const&, sim::index, sim::index)
    {
        return sim::harmonic_potential {.spring_constant = 1};
    }
};

int main()
{
    sim::system system;

    for (int i = 0; i < 100; i++) {
        system.add_particle(sim::basic_properties {
            .position = {i / 100.0, 0, 0}
        });
    }

    system.add_forcefield(std::make_shared<my_forcefield>());

    auto const report_energy = [&]() {
        std::cout << "Energy: "
                  << system.compute_energy()
                  << " (K = "
                  << system.compute_kinetic_energy()
                  << " | V = "
                  << system.compute_potential_energy()
                  << ")\n";
    };

    report_energy();

    sim::simulate_newtonian_dynamics(system, {
        .timestep = 0.001,
        .simulation_length = 1000000
    });

    report_energy();
}
