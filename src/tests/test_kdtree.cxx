#include <iostream>
#include <numeric>
#include <vector>
#include <random> 
#include <chrono>
#include <algorithm>
#include <string>
#include <map>
#include <tuple>

#include <NBodyMath.h>
#include <NBody.h>
#include <KDTree.h>

#include <omp.h>

#define PERIOD 1.0

#define start_time() auto linestart = __LINE__; auto start = std::chrono::high_resolution_clock::now();
#define end_time() auto lineend = __LINE__; auto end = std::chrono::high_resolution_clock::now();
#define get_time_taken() std::chrono::duration<double> elapsed = end-start; std::cout<<__func__<<"@L"<<linestart<<"-L"<<lineend<<" time taken "<<elapsed.count()<<std::endl;
#define log_time_taken(linestart, lineend, elapsed_time) std::cout<<__func__<<","<<__FILE__<<"@L"<<linestart<<"-L"<<lineend<<" time taken "<<elapsed_time<<std::endl;

int tpos=NBody::KDTree::TPHYS;
int tvel=NBody::KDTree::TVEL;
int tphs=NBody::KDTree::TPHS;

std::vector<Int_t> get_quants(Int_t N) {
    std::vector<Int_t> quants = {static_cast<Int_t>(0), static_cast<Int_t>(N*0.25), static_cast<Int_t>(N*0.5), static_cast<Int_t>(N*0.75), static_cast<Int_t>(N-1)};
    return quants;
}

std::map<std::string, std::tuple<int, int, double, double, double>> TreeTypes() 
{

    std::map<std::string, std::tuple<int, int, double, double, double>> treetypes;
    treetypes.insert(std::pair<std::string, std::tuple<int, int, double, double, double>>("Physical", {tpos, 16, -1, 0, 0}));
    treetypes.insert(std::pair<std::string, std::tuple<int, int, double, double, double>>("Physical Rdist", {tpos, 0, 0.01, 0, 0.01}));
    treetypes.insert(std::pair<std::string, std::tuple<int, int, double, double, double>>("Physical Rdist Adaptfac", {tpos, 0, 0.01, 0.1, 0.01}));
    treetypes.insert(std::pair<std::string, std::tuple<int, int, double, double, double>>("Velocity", {tvel, 16, -1, 0, 0}));
    treetypes.insert(std::pair<std::string, std::tuple<int, int, double, double, double>>("Phase", {tphs, 16, -1, 0, 0}));
    return treetypes;
}

std::vector<NBody::Particle> generate_vector(std::size_t size, double fac = 0.5, int nsub = 10)
{
    double period = PERIOD;
    Int_t nend = size*fac;
    auto disp = 0.1 * period / pow(static_cast<Double_t>(size), 1.0/3.0);
    std::vector<NBody::Particle> parts(size);
    if (nsub == 0) fac = 1.0;

    std::cout<<"Particle class requires "<<sizeof(NBody::Particle)<<" bytes "<<std::endl;
    std::cout<<"Total number of particles "<<size<<std::endl;
    std::cout<<"Total memory required is "<<size*sizeof(NBody::Particle)/1024./1024./1024.<<" GB "<<std::endl;
    std::cout<<"Uniform background is "<<fac<<" of total number of particles, "<<nend<<" particles"<<std::endl;
    std::cout<<nsub<<" substructures evenly containing "<<1.0 - fac<<" of total number of particles, "<<size-nend<<" particles, spanning a factor of  "<<disp/period<<" of the period"<<std::endl;

    start_time();
#if defined(_OPENMP)
#pragma omp parallel \
default(none) shared(parts, period, nend)
#endif
    {
    unsigned int seed = 4322;
#ifdef _OPENMP
    seed *= (omp_get_thread_num()+1);
#endif
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> pos(0, period);
    std::uniform_real_distribution<double> vel(0, period);
    std::vector<NBody::DoublePos_t> x(3), v(3);
    
#if defined(_OPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto i=0;i<nend; i++) {
        for (auto j = 0; j < 3; j++) {
            x[j] = pos(generator);
            v[j] = vel(generator);
        }
        parts[i].SetPosition(x.data());
        parts[i].SetVelocity(v.data());
    }
    }

    if (nsub == 0) return parts;
    // generate random positions for the subs
    unsigned int sub_seed = 4322;
    std::default_random_engine sub_generator(sub_seed);
    std::uniform_real_distribution<double> sub_pos(0, period);
    std::uniform_real_distribution<double> sub_vel(0, period);
    std::vector<NBody::DoublePos_t> sub_x(3*nsub), sub_v(3*nsub);
    for (auto i=0;i<nsub; i++) {
        for (auto j = 0; j < 3; j++) {
            sub_x[j+i*3] = sub_pos(sub_generator);
            sub_v[j+i+3] = sub_vel(sub_generator);
        }
    }
    
    Int_t ninsubs = (size - nend)/nsub;
#if defined(_OPENMP)
#pragma omp parallel \
default(none) shared(parts, period, nsub, sub_x, sub_v, disp, nend, ninsubs, size)
#endif
{
    #pragma omp for schedule(static,1)
    for (auto isub=0; isub<nsub; isub++)
    {
        unsigned int seed = 4322;
        std::default_random_engine generator(seed);
        std::normal_distribution<double> gauspos(0, disp);
        std::normal_distribution<double> gausvel(0, disp);
        std::vector<NBody::DoublePos_t> x(3), v(3);
        auto istart = nend + isub * ninsubs;
        auto iend = istart + ninsubs;
        if (isub == nsub -1) iend = size;
        for (auto i=istart;i<iend; i++) {
            for (auto j = 0; j < 3; j++) {
                x[j] = sub_x[j+isub*3]+gauspos(generator);
                v[j] = sub_v[j+isub*3]+gausvel(generator);
                if (x[j] > period) x[j] -= period;
                if (x[j]< 0 ) x[j] += period;
                if (v[j] > period) v[j] -= period;
                if (v[j]< 0 ) v[j] += period;
            }
            parts[i].SetPosition(x.data());
            parts[i].SetVelocity(v.data());
        }
    }
    }
    for (auto i=0;i<size;i++) parts[i].SetID(i);
    end_time();
    get_time_taken();
    return parts;
}

NBody::KDTree * build_kdtree(std::vector<NBody::Particle> &parts,
    Int_t b = 16, 
    int treetype = NBody::KDTree::TPHYS,
    double rdist2fac = -1, 
    double adaptivefactor = 0.0,
    double bfac = 0.0 
)
{
    auto N = static_cast<Int_t>(parts.size());
    if (b < 1) b = bfac*N;
    NBody::KDTree *tree = nullptr;
    // set some default stuff 
    std::vector<double> period(6);
    for (auto &p:period) p=1;
    Double_t *P = NULL; 
    // P=period.data();
    Double_t **m = nullptr;
    int split = 0; 
    int aniso = 0;
    int scalespace = 0;
    bool iBuildInParallel = true;
    bool iKeepInputOrder = false;

    std::cout<<"Building tree with [leaf node size, radius2, adaptivefactor]=["<<b<<", "<<rdist2fac<<", "<<adaptivefactor<<"]"<<std::endl;
    start_time();
    tree = new NBody::KDTree(parts.data(), N, 
            b, treetype, 
            tree->KEPAN, 100,
            split, aniso, scalespace,
            P, m,
            iBuildInParallel,
            iKeepInputOrder, 
            rdist2fac * PERIOD, 
            adaptivefactor
        );
    std::cout<<"Tree has "<<tree->GetNumNodes()<<" nodes, and "<<tree->GetNumLeafNodes()<<" leaf nodes "<<std::endl;
    end_time();
    get_time_taken();
    return tree;
}

void kdtree_test_NN(NBody::KDTree *&tree, 
    std::vector<NBody::Particle> &parts, 
    int num_nn = 16
)
{
    std::cout<<__func__<<" :: Search tree for NN = "<<num_nn<<std::endl;
    auto N = parts.size();
    std::vector<double> rdistave(N);
    std::vector<double> times(N);
    auto quants = get_quants(N);
    double sum = 0;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
    {
    std::vector<Int_t> nn(num_nn);
    std::vector<Double_t> nnd2(num_nn);
#if defined(_OPENMP)
    #pragma omp for schedule(guided) reduction(+:sum)
#endif
        for (Int_t i = 0; i < N; i++)
        {
            auto s = std::chrono::high_resolution_clock::now();
            tree->FindNearest(i, nn.data(), nnd2.data(), num_nn);
            auto e = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> el = e - s;
            times[i] = el.count();
            sum += el.count();
            rdistave[i] = 0;
            for (auto &d2:nnd2) {
                rdistave[i] += d2;
            }
            rdistave[i] /= static_cast<double>(num_nn);
        }
    }
    std::sort(times.begin(), times.end());
    std::cout<<"Processing times per particle statistics:";
    for (auto &x:quants) std::cout<<times[x]<<" ";
    std::cout<<" total "<<sum<<std::endl;
    std::sort(rdistave.begin(), rdistave.end());
    std::cout<<"NN dist^2 stats: ";
    for (auto &x:quants) std::cout<<rdistave[x]<<" ";
    std:cout<<std::endl;
}

void kdtree_test_ballsearch(NBody::KDTree *&tree, 
    std::vector<NBody::Particle> &parts, 
    double rdist = 0.1
)
{
    std::cout<<__func__<<" :: Search tree with distance "<<rdist<<std::endl;
    auto N = parts.size();
    std::vector<double> ninr2(N);
    auto r2 = rdist*rdist;
    std::vector<double> times(N);
    auto quants = get_quants(N);
    double sum = 0;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
    {
#if defined(_OPENMP)
    #pragma omp for schedule(guided) reduction(+:sum)
#endif
        for (Int_t i=0;i<N;i++) 
        {
            auto s = std::chrono::high_resolution_clock::now();
            auto neighbours = tree->SearchBallPosTagged(i, r2);
            auto e = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> el = e - s;
            times[i] = el.count();
            sum += el.count();
            ninr2[i] = neighbours.size();
        }
    }
    std::sort(times.begin(), times.end());
    std::cout<<"Processing times per particle statistics: ";
    for (auto &x:quants) std::cout<<times[x]<<" ";
    std::cout<<" total "<<sum<<std::endl;

    std::sort(ninr2.begin(), ninr2.end());
    std::cout<<"Number fo neighbours within dist^2 stats: ";
    for (auto &x:quants) std::cout<<ninr2[x]<<" ";
    std:cout<<std::endl;
}

void kdtree_test_FOF(NBody::KDTree *&tree, 
    std::vector<NBody::Particle> &parts, 
    Int_t minnum = 20,
    Double_t rdist = 0.1
)
{
    auto N = parts.size();
    Int_t numgroup = 0;
    // rdist = 1.0*PERIOD/pow(static_cast<Double_t>(N), 1.0/3.0);
    rdist *= PERIOD;
    std::cout<<__func__<<" :: FOF groups with distance "<<rdist<<" and min number of members "<<minnum<<std::endl;
    start_time();
    auto fofarr = tree->FOF(rdist, numgroup, minnum);
    end_time();
    std::cout<<"Number of groups "<<numgroup<<std::endl;
    if (numgroup>0) {
        std::vector<Int_t> fof(fofarr, fofarr + N);
        std::vector<Int_t> numingroup(numgroup+1);
        for (auto x:fof) numingroup[x]++;
        std::sort(numingroup.begin()+1, numingroup.end());
        std::vector<Int_t> quants = get_quants(numgroup);
        std::cout<<"Number of members in groups stats: ";
        for (auto &x:quants) std::cout<<numingroup[x+1]<<" ";
        std:cout<<std::endl;
        std::cout<<"Number of particles not in groups "<<numingroup[0]<<std::endl;
    }
    get_time_taken();
    delete[] fofarr;
}

int main(int argc, char *argv[])
{
    std::cerr << "Usage: " << argv[0] << " num_particles \n";
    std::size_t size = 1000000;
    if (argc == 2)
    {
        size = atoi(argv[1]);
    }
#ifdef _OPENMP 
    std::cout<<"Running with OpenMP "<<_OPENMP<<std::endl;
    std::cout<<"Using "<<omp_get_max_threads()<<" threads and "<<omp_get_max_active_levels()<<" level "<<std::endl;
#endif
    // generate vector 
    auto parts = generate_vector(size);

    std::cout<<std::endl;
    auto treetypes = TreeTypes();
    for (auto &m:treetypes) {
        std::cout<<"Building and testing a "<<m.first<<std::endl;
        std::cout<<"=========================="<<std::endl;
        auto [ttype, bsize, rdist2, adaptfac, bfac] = m.second;
        auto tree = build_kdtree(parts, bsize, ttype, rdist2, adaptfac, bfac);
        kdtree_test_NN(tree, parts);
        kdtree_test_ballsearch(tree, parts);
        // NBody::nbody_total_time = 0;
        // for (auto &x:NBody::nbody_counter) x = 0;
        kdtree_test_FOF(tree, parts);
        // std::cout<<NBody::nbody_total_time<<std::endl;
        // for (auto &x:NBody::nbody_counter) std::cout<<x<<" ";std::cout<<std::endl;
        delete tree;
        std::cout<<std::endl;
    }

}