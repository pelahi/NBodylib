/*! \file SwiftParticle.h
 *  \brief header file for the SWIFT particle type.

*/

#ifndef SWIFT_PARTICLE_H
#define SWIFT_PARTICLE_H

#ifdef SWIFTINTERFACE

namespace Swift
{
    #define SWIFT_STRUCT_ALIGNMENT 32
    #define SWIFT_STRUCT_ALIGN __attribute__((aligned(SWIFT_STRUCT_ALIGNMENT)))

    /* SWIFT enum of part types. Should match VELOCIraptor type values. */
    enum part_type {
        swift_type_gas = 0,
        swift_type_dark_matter = 1,
        swift_type_star = 4,
        swift_type_black_hole = 5,
        swift_type_count
    } __attribute__((packed));

    enum swift_vel_chemistry_element {
        chemistry_element_H = 0,
        chemistry_element_He,
        chemistry_element_C,
        chemistry_element_N,
        chemistry_element_O,
        chemistry_element_Ne,
        chemistry_element_Mg,
        chemistry_element_Si,
        chemistry_element_Fe,
        chemistry_element_count
    }__attribute__((packed));

    /* SWIFT/VELOCIraptor particle. */
    struct swift_vel_part {

      /*! Particle ID. */
      long long id;

      /*! Particle position. */
      double x[3];

      /*! Particle velocity. */
      float v[3];

      /*! Particle mass. */
      #ifndef NOMASS
      float mass;
      #endif

      /*! Gravitational potential */
      float potential;

      /*! Internal energy of gas particle */
      float u;

      /*! Temperature of a gas particle */
      float T;

      /*! Type of the #gpart (DM, gas, star, ...) */
      enum part_type type;

      /*! swift task */
      int task;

      /*! swift index */
      int index;

    };

    /* SWIFT/VELOCIraptor chemistry data. */
    struct swift_vel_chemistry_data {

        /*! Fraction of the particle mass in a given element */
        float metal_mass_fraction[chemistry_element_count];

        /*! Fraction of the particle mass in *all* metals */
        float metal_mass_fraction_total;

        /*! Mass coming from SNIa */
        float mass_from_SNIa;

        /*! Fraction of total gas mass in metals coming from SNIa */
        float metal_mass_fraction_from_SNIa;

        /*! Mass coming from AGB */
        float mass_from_AGB;

        /*! Fraction of total gas mass in metals coming from AGB */
        float metal_mass_fraction_from_AGB;

        /*! Mass coming from SNII */
        float mass_from_SNII;

        /*! Fraction of total gas mass in metals coming from SNII */
        float metal_mass_fraction_from_SNII;

        /*! Fraction of total gas mass in Iron coming from SNIa */
        float iron_mass_fraction_from_SNIa;

    };

    /* SWIFT/VELOCIraptor chemistry of black holes. */
    struct swift_vel_chemistry_bh_data {

        /*! Mass in a given element */
        float metal_mass[chemistry_element_count];

        /*! Mass in *all* metals */
        float metal_mass_total;

        /*! Mass coming from SNIa */
        float mass_from_SNIa;

        /*! Mass coming from AGB */
        float mass_from_AGB;

        /*! Mass coming from SNII */
        float mass_from_SNII;

        /*! Metal mass coming from SNIa */
        float metal_mass_from_SNIa;

        /*! Metal mass coming from AGB */
        float metal_mass_from_AGB;

        /*! Metal mass coming from SNII */
        float metal_mass_from_SNII;

        /*! Iron mass coming from SNIa */
        float iron_mass_from_SNIa;
    };


    /* SWIFT/VELOCIraptor gas particle. */
    struct swift_vel_gas_part {
        /*! Particle smoothing length. */
        float h;

        /*! Particle internal energy. */
        float u;

        /*! Time derivative of the internal energy. */
        float u_dt;

        /*! Particle density. */
        float rho;

        /*! Particle pressure (weighted) */
        float pressure_bar;

        /* Store viscosity information in a separate struct. */
        struct {

        /*! Particle velocity divergence */
        float div_v;

        /*! Particle velocity divergence from previous step */
        float div_v_previous_step;

        /*! Artificial viscosity parameter */
        float alpha;

        /*! Signal velocity */
        float v_sig;

        } viscosity;

        /* Store thermal diffusion information in a separate struct. */
        struct {

        /*! del^2 u, a smoothed quantity */
        float laplace_u;

        /*! Thermal diffusion coefficient */
        float alpha;

        } diffusion;

        /* Chemistry information */
        struct swift_vel_chemistry_data chemistry_data;

        /*! swift index */
        int index;
    };


    /* SWIFT/VELOCIraptor star particle. */
    struct swift_vel_star_part {

        /*! Birth time (or scalefactor)*/
        float birth_time;

        /*! Birth density */
        float birth_density;

        /*! Birth temperature */
        float birth_temperature;

        /*! Feedback energy fraction */
        float f_E;

        /*! Chemistry structure */
        struct swift_vel_chemistry_data chemistry_data;

        /*! swift index */
        int index;
    };

    /* SWIFT/VELOCIraptor black hole particle. */
    struct swift_vel_bh_part {

        /*! Formation time (or scale factor)*/
        float formation_time;

        /*! Subgrid mass of the black hole */
        float subgrid_mass;

        /*! Total accreted mass of the black hole (including accreted mass onto BHs
        * that were merged) */
        float total_accreted_mass;

        /*! Energy reservoir for feedback */
        float energy_reservoir;

        /*! Instantaneous accretion rate */
        float accretion_rate;

        /*! Density of the gas surrounding the black hole. */
        float rho_gas;

        /*! Smoothed sound speed of the gas surrounding the black hole. */
        float sound_speed_gas;

        /*! Smoothed velocity (peculiar) of the gas surrounding the black hole */
        float velocity_gas[3];

        /*! Curl of the velocity field around the black hole */
        float circular_velocity_gas[3];

        /*! Number of seeds in this BH (i.e. itself + the merged ones) */
        int cumulative_number_seeds;

        /*! Total number of BH merger events (i.e. not including all progenies) */
        int number_of_mergers;

        /*! Chemistry information (e.g. metal content at birth, swallowed metal
        * content, etc.) */
        struct swift_vel_chemistry_bh_data chemistry_data;

        /*! swift index */
        int index;
    };
}

#endif

#endif
