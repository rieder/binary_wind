"""
Script to simulate a binary star with a stellar wind, using a hydro code.
"""

import argparse
import matplotlib.pyplot as plt

from amuse.datamodel import Particle
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
        default="huayno", type=str,
        help="Gravity code to use",
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


def new_hydro_instance(name, converter):
    name = name.lower()
    if name == "phantom":
        instance = Phantom(converter)
    elif name == "gadget2":
        instance = Gadget2(converter)
    elif name == "arepo":
        instance = Arepo(converter)
    elif name == "fi":
        instance = Fi(converter, mode="openmp")
    else:
        raise Exception("Unknown hydro code")
    return instance


def main():
    args = new_argument_parser()
    stars = new_binary_from_orbital_elements(
        args.mass_primary,
        args.mass_secondary,
        args.semi_major_axis,
        args.eccentricity,
    )

    converter = nbody_system.nbody_to_si(
        stars.total_mass(),
        args.semi_major_axis,
    )
    system = Bridge()
    grav = new_gravity_instance(args.grav, converter)
    system.add_system(grav, ())
    stars_in_grav = grav.particles.add_particles(stars)
    evo = new_stellar_evolution_instance(args.evo)

    binary = Particle()
    binary.semi_major_axis = args.semi_major_axis
    binary.eccentricity = args.eccentricity
    binary.child1 = stars[0]
    binary.child2 = stars[1]
    stars_in_evo = evo.particles.add_particles(stars)
    binary_in_evo = evo.binaries.add_particle(binary)

    channel_from_evo = stars_in_evo.new_channel_to(stars)
    channel_from_grav = stars_in_grav.new_channel_to(stars)
    channel_to_grav = stars.new_channel_to(stars_in_grav)
    channel_to_evo = stars.new_channel_to(stars_in_evo)
    channel_from_evo.copy()

    time = 0 | units.yr
    time_step = stars.time_step.min()
    t = [] | units.yr
    m1 = [] | units.MSun
    m2 = [] | units.MSun
    a = [] | units.RSun
    e = []
    transients = [0, 13, 14] | units.stellar_type
    while (
        stars[0].stellar_type not in transients
        or stars[1].stellar_type not in transients
    ):
        time += time_step/2
        evo.evolve_model(time)
        channel_from_evo.copy()
        channel_to_grav.copy()
        time += time_step/2
        system.evolve_model(time)
        channel_to_evo.copy()
        evo.evolve_model(time)
        channel_from_grav.copy()
        channel_from_evo.copy()
        channel_to_grav.copy()
        time_step = stars.time_step.min()
        (
            mass1, mass2, semimajor_axis, eccentricity,
            true_anomaly, inclination, long_asc_node, arg_per
        ) = get_orbital_elements_from_binary(stars)
        t.append(time)
        m1.append(mass1)
        m2.append(mass2)
        a.append(semimajor_axis)
        e.append(eccentricity)
        print(
            f"{time} "
            f"{mass1.in_(units.MSun)} "
            f"{mass2.in_(units.MSun)} "
            f"{semimajor_axis.in_(units.RSun)} "
            f"{eccentricity} "
            f"{stars[0].stellar_type.value_in(units.stellar_type)} "
            f"{stars[0].stellar_type} "
            f"{stars[1].stellar_type.value_in(units.stellar_type)} "
            f"{stars[1].stellar_type} "
        )
    fig = plt.figure(figsize=(8, 8))
    ax = []
    time_label = f'time ({t.unit})'
    for i in range(4):
        ax.append(fig.add_subplot(2, 2, 1+i))
    # ax 0: semi-major axes
    ax[0].plot(
        t.value_in(time.unit),
        a.value_in(units.RSun),
    )
    ax[0].set_xlabel(time_label)
    ax[0].set_ylabel('semi-major axis (RSun)')
    ax[0].set_xscale('log')
    # ax[0].set_yscale('log')

    # ax 1: eccentricity
    ax[1].plot(
        t.value_in(time.unit),
        e,
    )
    ax[1].set_xlabel(time_label)
    ax[1].set_ylabel('eccentricity')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    # ax 2: masses
    ax[2].plot(
        t.value_in(time.unit),
        m1.value_in(units.MSun),
    )
    ax[2].plot(
        t.value_in(time.unit),
        m2.value_in(units.MSun),
    )
    ax[2].set_xlabel(time_label)
    ax[2].set_ylabel(f'Mass ({units.MSun})')
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')

    # ax 3: ?

    plt.savefig('plot.pdf')
    return


if __name__ == "__main__":
    main()
