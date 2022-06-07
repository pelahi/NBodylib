#include <iostream>
#include <numeric>
#include <vector>
#include <random> 
#include <chrono>

#include <NBodyMath.h>
#include <NBody.h>
#include <KDTree.h>


#define start_time() auto linestart = __LINE__; auto start = std::chrono::high_resolution_clock::now();
#define end_time() auto lineend = __LINE__; auto end = std::chrono::high_resolution_clock::now();
#define get_time_taken() std::chrono::duration<double> elapsed = end-start; std::cout<<__func__<<"@L"<<linestart<<"-L"<<lineend<<" time taken "<<elapsed.count()<<std::endl;

std::vector<NBody::Particle> generate_vector(std::size_t size)
{
    std::cout<<"Particle class requires "<<sizeof(NBody::Particle)<<" bytes "<<std::endl;
    std::cout<<"Total number of particles "<<size<<std::endl;
    std::cout<<"Total memory required is "<<size*sizeof(NBody::Particle)/1024./1024./1024.<<" GB "<<std::endl;
    std::vector<NBody::Particle> parts(size);
    start_time();
#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    std::default_random_engine generator;
    std::uniform_real_distribution<double> pos(-1,1);
    std::uniform_real_distribution<double> vel(-1,1);
    std::vector<NBody::DoublePos_t> x(3), v(3);
    
#if defined(USEOPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto &p:parts) {
        for (auto j = 0; j < 3; j++) {
            x[j] = pos(generator);
            v[j] = vel(generator);
        }
        p.SetPosition(x.data());
        p.SetVelocity(v.data());
    }
#if defined(USEOPENMP)
}
#endif
    end_time();
    get_time_taken();
    return parts;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " num_particles \n";
        return 1;
    }
    // generate vector 
    std::size_t size = atoi(argv[1]);
    auto parts = generate_vector(size);
}