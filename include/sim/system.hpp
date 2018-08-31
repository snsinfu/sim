// system - collection of particles, their properties and forcefield

// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef INCLUDED_SIM_SYSTEM_HPP
#define INCLUDED_SIM_SYSTEM_HPP

#include <memory>
#include <typeindex>
#include <utility>
#include <unordered_map>
#include <vector>


#include "forcefield.hpp"
#include "typedef.hpp"


namespace sim
{
    template<typename T>
    struct basic_property_traits
    {
        using value_type = T;

        static constexpr value_type default_value()
        {
            return value_type{};
        }
    };

    template<typename T>
    struct property_traits : basic_property_traits<T>
    {
    };

    template<typename T>
    using property_value_t = typename property_traits<T>::value_type;

    namespace detail
    {
        struct mass_property
        {
        };

        struct position_property
        {
        };

        struct velocity_property
        {
        };
    }

    template<>
    struct property_traits<detail::mass_property> : sim::basic_property_traits<sim::scalar>
    {
        static constexpr value_type default_value()
        {
            return 1;
        }
    };

    template<>
    struct property_traits<detail::position_property> : sim::basic_property_traits<sim::point>
    {
    };

    template<>
    struct property_traits<detail::velocity_property> : sim::basic_property_traits<sim::vector>
    {
    };

    namespace detail
    {
        class property_base
        {
        public:
            virtual ~property_base() = default;
            virtual void resize(sim::index n) = 0;
        };

        template<typename T>
        class property : public property_base
        {
        public:
            using value_type = sim::property_value_t<T>;

            void resize(sim::index n) override
            {
                values_.resize(n, sim::property_traits<T>::default_value());
            }

            sim::array_view<value_type> view() noexcept
            {
                return values_;
            }

        private:
            std::vector<value_type> values_;
        };

        class particle_array
        {
        public:
            sim::index size() const noexcept
            {
                return size_;
            }

            void resize(sim::index n)
            {
                for (auto& prop : properties_) {
                    prop.second->resize(n);
                }
                size_ = n;
            }

            template<typename T>
            sim::array_view<sim::property_value_t<T>> require_property_array()
            {
                std::type_index const key = typeid(T);

                if (properties_.find(key) == properties_.end()) {
                    auto prop = std::make_unique<detail::property<T>>();
                    prop->resize(size_);
                    properties_.emplace(key, std::move(prop));
                }

                return property_array<T>();
            }

            template<typename T>
            sim::array_view<sim::property_value_t<T>> property_array()
            {
                return dynamic_cast<detail::property<T>&>(*properties_.at(typeid(T))).view();
            }

            template<typename T>
            sim::array_view<sim::property_value_t<T> const> property_array() const
            {
                return dynamic_cast<detail::property<T>&>(*properties_.at(typeid(T))).view();
            }

        private:
            sim::index size_ = 0;
            std::unordered_map<std::type_index, std::unique_ptr<detail::property_base>> properties_;
        };
    }

    struct basic_properties
    {
        sim::scalar mass = sim::property_traits<detail::mass_property>::default_value();
        sim::point position;
        sim::vector velocity;
    };

    class system
    {
    public:
        system()
        {
            require_property_array<detail::mass_property>();
            require_property_array<detail::position_property>();
            require_property_array<detail::velocity_property>();
        }

        void add_particle(basic_properties const& props = {})
        {
            auto const index = particles_.size();

            particles_.resize(particles_.size() + 1);

            mass_array()[index] = props.mass;
            position_array()[index] = props.position;
            velocity_array()[index] = props.velocity;
        }

        sim::index particle_count() const noexcept
        {
            return particles_.size();
        }

        template<typename T>
        sim::array_view<sim::property_value_t<T>> require_property_array()
        {
            return particles_.require_property_array<T>();
        }

        template<typename T>
        sim::array_view<sim::property_value_t<T>> property_array()
        {
            return particles_.property_array<T>();
        }

        template<typename T>
        sim::array_view<sim::property_value_t<T> const> property_array() const
        {
            return particles_.property_array<T>();
        }

        sim::array_view<sim::scalar> mass_array() noexcept
        {
            return property_array<detail::mass_property>();
        }

        sim::array_view<sim::scalar const> mass_array() const noexcept
        {
            return property_array<detail::mass_property>();
        }

        sim::array_view<sim::point> position_array() noexcept
        {
            return property_array<detail::position_property>();
        }

        sim::array_view<sim::point const> position_array() const noexcept
        {
            return property_array<detail::position_property>();
        }

        sim::array_view<sim::vector> velocity_array() noexcept
        {
            return property_array<detail::velocity_property>();
        }

        sim::array_view<sim::vector const> velocity_array() const noexcept
        {
            return property_array<detail::velocity_property>();
        }

        void add_forcefield(std::shared_ptr<sim::forcefield> forcefield)
        {
            forcefield_.add_component(forcefield);
        }

        sim::scalar compute_kinetic_energy()
        {
            sim::scalar kinetic = 0;

            auto const n = particle_count();
            auto const masses = mass_array();
            auto const velocities = velocity_array();

            for (sim::index i = 0; i < n; i++) {
                kinetic += 0.5 * masses[i] * velocities[i].squared_norm();
            }

            return kinetic;
        }

        sim::scalar compute_potential_energy()
        {
            return forcefield_.compute_energy(*this);
        }

        sim::scalar compute_energy()
        {
            return compute_kinetic_energy() + compute_potential_energy();
        }

        void compute_force(sim::array_view<sim::vector> forces)
        {
            for (auto& force : forces) {
                force = {};
            }
            return forcefield_.compute_force(*this, forces);
        }

    private:
        detail::particle_array particles_;
        sim::composite_forcefield forcefield_;
    };
}

#endif
