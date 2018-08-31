#include <iostream>
#include <memory>
#include <random>

#include <sim/all.hpp>


class my_forcefield : public sim::forcefield
{
public:
    sim::scalar spring_constant = 1;

    sim::scalar compute_energy(sim::system const& system) override
    {
        auto const particle_count = system.particle_count();
        auto const positions = system.position_array();

        sim::scalar energy = 0;

        for (sim::index i = 0; i < particle_count; i++) {
            energy += 0.5 * spring_constant * positions[i].squared_distance({});
        }

        return energy;
    }

    void compute_force(sim::system const& system, sim::array_view<sim::vector> forces) override
    {
        auto const particle_count = system.particle_count();
        auto const positions = system.position_array();

        for (sim::index i = 0; i < particle_count; i++) {
            forces[i] += -spring_constant * positions[i].vector();
        }
    }
};

int main()
{
    sim::system system;

    for (int i = 0; i < 100; i++) {
        system.add_particle({
            .position = {i / 100.0, 0, 0}
        });
    }

    system.add_forcefield(std::make_shared<my_forcefield>());

    std::cout << "Energy: "
              << system.compute_energy()
              << " (K = "
              << system.compute_kinetic_energy()
              << " | V = "
              << system.compute_potential_energy()
              << ")\n";

    sim::simulate_newtonian_dynamics(system, {
        .timestep = 0.001,
        .simulation_length = 1000000
    });

    std::cout << "Energy: "
              << system.compute_energy()
              << " (K = "
              << system.compute_kinetic_energy()
              << " | V = "
              << system.compute_potential_energy()
              << ")\n";
}
