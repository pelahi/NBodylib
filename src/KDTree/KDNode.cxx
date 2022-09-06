/*! \file KDNode.cxx
 *  \brief This file contains implementation of Node functions

*/

#include <iostream>
#include <KDNode.h>

namespace NBody
{
    /// \defgroup GetPartDataForTree 
    /// Get particle data relevant to tree 
    //@{
    #pragma omp declare simd
    DoublePos_t Node::get_particle_pos_jth(Particle &p, int j) {
        return p.GetPosition(j);
    }
    #pragma omp declare simd
    DoublePos_t Node::get_particle_vel_jth(Particle &p, int j) {
        return p.GetVelocity(j);
    }
    #pragma omp declare simd
    DoublePos_t Node::get_particle_phs_jth(Particle &p, int j) {
        return p.GetPhase(j);
    }
    void Node::set_particle_pos_jth(Particle &p, int j, Double_t x) {
        p.SetPosition(j, x);
    }
    void Node::set_particle_vel_jth(Particle &p, int j, Double_t x) {
        p.SetVelocity(j, x);
    }
    void Node::set_particle_phs_jth(Particle &p, int j, Double_t x) {
        p.SetPhase(j, x);
    }
    Double_t Node::D2BetweenParticles(Particle &p1, Particle &p2, int dim) 
    {
        Double_t d2=0;
        for (auto i=0;i<dim;i++) {
            auto diff = (this->*get_part_data_jth)(p1,i)-(this->*get_part_data_jth)(p2,i);
            d2 += diff * diff;
        }
        return d2;
    }
    Double_t Node::D2ToParticle(Double_t *x, Particle &p2, int dim) 
    {
        Double_t d2=0;
        for (auto i=0;i<dim;i++) {
            auto diff = x[i]-(this->*get_part_data_jth)(p2,i);
            d2 += diff * diff;
        }
        return d2;
    }
    //@}

}
