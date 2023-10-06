"""
Script to simulate a binary star with a stellar wind, using a hydro code.
"""

import argparse
import matplotlib.pyplot as plt

from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.ext.orbital_elements import (
    new_binary_from_orbital_elements,
    get_orbital_elements_from_binary,
)
from amuse.couple.bridge import Bridge
from amuse.ext.stellar_wind import new_stellar_wind
# stellar evolution codes
from amuse.community.seba import Seba
# gravity codes
try:
    from amuse.community.huayno import Huayno
except ImportError:
    Huayno = None
try:
    from amuse.community.hermite import Hermite
except ImportError:
    Hermite = None
try:
    from amuse.community.ph4 import Ph4
except ImportError:
    Ph4 = None
from amuse.community.smalln import Smalln
# hydro codes
from amuse.community.fi import Fi
try:
    from amuse.community.gadget2 import Gadget2
except:
    Gadget2 = None
try:
    from amuse.community.phantom import Phantom
except:
    Phantom = None
try:
    from amuse.community.arepo import Arepo
except ImportError:
    Arepo = None
# from amuse.ext.ekster.gas_class import GasCode


def new_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mp", dest="mass_primary",
        default=40 | units.MSun, type=units.MSun,
        help="Mass of the primary star",
    )
    parser.add_argument(
        "--ms", dest="mass_secondary",
        default=20 | units.MSun, type=units.MSun,
        help="Mass of the secondary star",
    )
    parser.add_argument(
        "-a", dest="semi_major_axis",
        default=100 | units.RSun, type=units.RSun,
        help="Semi-major axis",
    )
    parser.add_argument(
        "-e", dest="eccentricity",
        default=0, type=float,
        help="Eccentricity",
    )
    parser.add_argument(
        "--evo", dest="evo",
        default="seba", type=str,
        help="Stellar evolution code to use",
    )
    parser.add_argument(
        "--grav", dest="grav",
        default="smalln", type=str,
        help="Gravity code to use",
    )
    parser.add_argument(
        "--hydro", dest="hydro",
        default="phantom", type=str,
        help="Hydro code to use",
    )

    return parser.parse_args()


def new_stellar_evolution_instance(name):
    name = name.lower()
    if name == "seba":
        instance = Seba()
    else:
        raise Exception("Unknown stellar evolution code")
    return instance


def new_gravity_instance(name, converter):
    name = name.lower()
    if name == "huayno":
        instance = Huayno(converter, redirection="none")
    elif name == "hermite":
        instance = Hermite(converter)
    elif name == "ph4":
        instance = Ph4(converter)
    elif name == "smalln":
        instance = Smalln(converter)
    else:
        raise Exception("Unknown gravity code")
    return instance


def new_hydro_instance(name, *args, **kwargs):
    name = name.lower()
    if name == "phantom":
        code = Phantom
    elif name == "gadget2":
        code = Gadget2
    elif name == "arepo":
        code = Arepo
    elif name == "fi":
        code = Fi
    else:
        raise Exception("Unknown hydro code")
    return GasCode(sph_code=code, *args, **kwargs)


class BinaryWithWind:
    def __init__(
        self,
        mass1,
        mass2,
        semi_major_axis,
        eccentricity,
        evo_code_name="Seba",
        grav_code_name="Smalln",
        hydro_code_name="Fi",
    ):
        self.debug = 1
        sph_particle_mass = 0.0001 | units.MSun
        self.model_time = 0 | units.yr
        self.stars = new_binary_from_orbital_elements(
            mass1,
            mass2,
            semi_major_axis,
            eccentricity,
        )
        self.binary = Particles(
            semi_major_axis=semi_major_axis,
            eccentricity=eccentricity,
            child1=self.stars[0],
            child2=self.stars[1],
        )
        self.converter = nbody_system.nbody_to_si(
            self.stars.total_mass(),
            semi_major_axis,
        )
        self.stellar_evolution = new_stellar_evolution_instance(
            evo_code_name,
        )
        self.stars_in_stellar_evolution = \
            self.stellar_evolution.particles.add_particles(
                self.stars
            )
        self.binary_in_stellar_evolution = \
            self.stellar_evolution.binaries.add_particle(
                self.binary
            )
        self.gravity = new_gravity_instance(
            grav_code_name, self.converter
        )
        self.stars_in_gravity = self.gravity.particles.add_particles(
            self.stars
        )
        # self.hydro = new_hydro_instance(hydro_code_name, self.converter)

        self.system = Bridge()
        self.system.add_system(self.gravity, ())
        # self.system.add_system(self.gravity, (self.hydro,))
        # self.system.add_system(self.hydro, (self.gravity,))

        self.wind = new_stellar_wind(
            sph_particle_mass,
            # target_gas=self.hydro.gas_particles,
            # timestep=gas_timestep,
            mode="simple",
            # mode="accelerating",
            derive_from_evolution=True,
        )
        self.wind.particles.add_particles(self.stars)

        self.wind_gas = Particles()

        self.stars_in_stellar_evolution.new_channel_to(self.stars).copy()

    def evolve_model(self, time):
        "Evolve model forward to 'time'"
        time_step = time - self.model_time
        if time_step <= 0 | units.yr:
            print("error: time step too small or negative")
            return
        # if self.debug:
        #     print(f"EVOLVING SE TO {self.model_time + time_step/2}")
        self.stellar_evolution.evolve_model(
            self.model_time + time_step/2
        )
        self.stellar_evolution.particles.new_channel_to(
            self.stars
        ).copy()
        # if self.debug:
        #     print(f"EVOLVING WIND TO {self.model_time + time_step/2}")
        self.stars.new_channel_to(self.wind.particles).copy()
        print(
            f"{self.wind.particles[0].lost_mass} "
            f"{self.wind.particles[1].lost_mass} "
        )
        self.wind.evolve_model(
            self.model_time + time_step/2
        )
        # Only add new wind particles here - just before doing a hydro step
        if self.wind.has_new_wind_particles():
            # if self.debug:
            #     print("CREATING NEW WIND PARTICLES")
            wind_gas = self.wind.create_wind_particles()
            self.wind_gas.add_particles(wind_gas)
            # self.hydro.gas_particles.add_particles(wind_gas)
            # print(self.hydro.gas_particles.total_mass())
        self.stars.new_channel_to(
            self.gravity.particles
        ).copy()
        # if self.debug:
        #     print(f"EVOLVING SYSTEM TO {time}")
        self.system.evolve_model(time)
        self.gravity.particles.new_channel_to(
            self.stars
        ).copy()
        (
            mass1, mass2, semi_major_axis, eccentricity,
            true_anomaly, inclination, long_asc_node, arg_per
        ) = get_orbital_elements_from_binary(self.stars)
        self.binary.semi_major_axis = semi_major_axis
        self.binary.eccentricity = eccentricity
        self.binary.true_anomaly = true_anomaly
        self.binary.inclination = inclination
        self.binary.long_asc_node = long_asc_node
        self.binary.arg_per = arg_per
        self.binary.new_channel_to(
            self.stellar_evolution.binaries
        ).copy()
        self.stellar_evolution.evolve_model(time)
        self.stellar_evolution.particles.new_channel_to(
            self.stars
        ).copy()
        self.model_time = time
        return


def main():
    args = new_argument_parser()

    system = BinaryWithWind(
        args.mass_primary,
        args.mass_secondary,
        args.semi_major_axis,
        args.eccentricity,
        evo_code_name=args.evo,
        grav_code_name=args.grav,
        hydro_code_name=args.hydro,
    )

    t = [] | units.yr
    m1 = [] | units.MSun
    m2 = [] | units.MSun
    a = [] | units.RSun
    e = []
    transients = [0, 13, 14] | units.stellar_type
    while (
        system.stars[0].stellar_type not in transients
        or system.stars[1].stellar_type not in transients
    ):
        time_step = system.stars.time_step.min()
        system.evolve_model(system.model_time + time_step)
        t.append(system.model_time)
        m1.append(system.stars[0].mass)
        m2.append(system.stars[1].mass)
        a.append(system.binary[0].semi_major_axis)
        e.append(system.binary[0].eccentricity)
        print(
            f"{t[-1]} "
            f"{m1[-1].in_(units.MSun)} "
            f"{m2[-1].in_(units.MSun)} "
            f"{a[-1].in_(units.RSun)} "
            f"{e[-1]} "
            f"{system.stars[0].stellar_type.value_in(units.stellar_type)} "
            # f"{system.stars[0].stellar_type} "
            f"{system.stars[1].stellar_type.value_in(units.stellar_type)} "
            # f"{system.stars[1].stellar_type} "
            f"{len(system.wind_gas)} "
            f"{system.wind.particles[0].lost_mass} "
            f"{system.wind.particles[1].lost_mass} "
        )
    fig = plt.figure(figsize=(8, 8))
    ax = []
    time_label = f'time ({t.unit})'
    for i in range(4):
        ax.append(fig.add_subplot(2, 2, 1+i))
    # ax 0: semi-major axes
    ax[0].plot(
        t.value_in(system.model_time.unit),
        a.value_in(units.RSun),
    )
    ax[0].set_xlabel(time_label)
    ax[0].set_ylabel('semi-major axis (RSun)')
    # ax[0].set_xscale('log')
    # ax[0].set_yscale('log')

    # ax 1: eccentricity
    ax[1].plot(
        t.value_in(system.model_time.unit),
        e,
    )
    ax[1].set_xlabel(time_label)
    ax[1].set_ylabel('eccentricity')
    # ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    # ax 2: masses
    ax[2].plot(
        t.value_in(system.model_time.unit),
        m1.value_in(units.MSun),
    )
    ax[2].plot(
        t.value_in(system.model_time.unit),
        m2.value_in(units.MSun),
    )
    ax[2].set_xlabel(time_label)
    ax[2].set_ylabel(f'Mass ({units.MSun})')
    # ax[2].set_xscale('log')
    ax[2].set_yscale('log')

    # ax 3: ?

    plt.savefig('plot.pdf')
    return


if __name__ == "__main__":
    main()
