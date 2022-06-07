#include <iostream>
#include <numeric>
#include <vector>
#include <random> 
#include <chrono>
#include <algorithm>
#include <string>
#include <map>

#include <NBodyMath.h>
#include <NBody.h>
#include <KDTree.h>

#include <omp.h>



#define start_time() auto linestart = __LINE__; auto start = std::chrono::high_resolution_clock::now();
#define end_time() auto lineend = __LINE__; auto end = std::chrono::high_resolution_clock::now();
#define get_time_taken() std::chrono::duration<double> elapsed = end-start; std::cout<<__func__<<"@L"<<linestart<<"-L"<<lineend<<" time taken "<<elapsed.count()<<std::endl;

int tpos=NBody::KDTree::TPHYS;
int tvel=NBody::KDTree::TVEL;
int tphs=NBody::KDTree::TPHS;

std::map<std::string, int> TreeTypes() 
{
    std::map<std::string, int> treetypes;
    treetypes.insert(std::pair<std::string, int>("Physical", tpos));
    treetypes.insert(std::pair<std::string, int>("Velocity", tvel));
    treetypes.insert(std::pair<std::string, int>("Phase", tphs));
    return treetypes;
}

std::vector<NBody::Particle> generate_vector(std::size_t size)
{
    std::cout<<"Particle class requires "<<sizeof(NBody::Particle)<<" bytes "<<std::endl;
    std::cout<<"Total number of particles "<<size<<std::endl;
    std::cout<<"Total memory required is "<<size*sizeof(NBody::Particle)/1024./1024./1024.<<" GB "<<std::endl;
    std::vector<NBody::Particle> parts(size);
    start_time();
#if defined(_OPENMP)
#pragma omp parallel default(shared)
{
#endif
    unsigned int seed = 4322;
#ifdef _OPENMP
    seed *= (omp_get_thread_num()+1);
#endif
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> pos(-1,1);
    std::uniform_real_distribution<double> vel(-1,1);
    std::vector<NBody::DoublePos_t> x(3), v(3);
    
#if defined(_OPENMP)
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
#if defined(_OPENMP)
}
#endif
    end_time();
    get_time_taken();
    return parts;
}

NBody::KDTree * build_kdtree(std::vector<NBody::Particle> &parts,
    Int_t b = 16, 
    int treetype = NBody::KDTree::TPHYS,
    double rdist2 = -1
)
{
    NBody::KDTree *tree = nullptr;
    // set some default stuff 
    std::vector<double> period(6);
    for (auto &p:period) p=1;
    Double_t *P = NULL; // P=period.data();
    Double_t **m = nullptr;
    int split = 0; 
    int aniso = 0;
    int scalespace = 0;
    bool iBuildInParallel = true;
    bool iKeepInputOrder = false;

    std::cout<<"Building tree with "<<b<<" leaf node size"<<std::endl;
    start_time();
    tree = new NBody::KDTree(parts.data(), static_cast<Int_t>(parts.size()), b, treetype, 
            tree->KEPAN, 100,
            split, aniso, scalespace,
            P, m,
            iBuildInParallel,
            iKeepInputOrder
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
    std::cout<<"Search tree for NN = "<<num_nn<<std::endl;
    auto N = parts.size();
    std::vector<double> rdistave(N);
    start_time();
#if defined(_OPENMP)
#pragma omp parallel default(shared)
{
#endif
    std::vector<Int_t> nn(num_nn);
    std::vector<Double_t> nnd2(num_nn);
#if defined(_OPENMP)
    #pragma omp for schedule(guided)
#endif
    for (Int_t i=0;i<N;i++) 
    {
        tree->FindNearest(i, nn.data(), nnd2.data(), num_nn);
        rdistave[i] = 0;
        for (auto &d2:nnd2) {
            rdistave[i] += d2;
        }
        rdistave[i] /= static_cast<double>(num_nn);
    }
#if defined(_OPENMP)
}
#endif
    std::sort(rdistave.begin(), rdistave.end());
    std::cout<<"Median distance of NN search is "<<rdistave[N/2]<<std::endl;
    end_time();
    get_time_taken();
}

void kdtree_test_ballsearch(NBody::KDTree *&tree, 
    std::vector<NBody::Particle> &parts, 
    double rdist = 0.1
)
{
    std::cout<<"Search tree with distance "<<rdist<<std::endl;
    auto N = parts.size();
    std::vector<double> ninr2(N);
    start_time();
    auto r2 = rdist*rdist;
#if defined(_OPENMP)
#pragma omp parallel default(shared)
{
#endif
#if defined(_OPENMP)
    #pragma omp for schedule(guided)
#endif
    for (Int_t i=0;i<N;i++) 
    {
        auto neighbours = tree->SearchBallPosTagged(i, r2);
        ninr2[i] = neighbours.size();
    }
#if defined(_OPENMP)
}
#endif
    std::sort(ninr2.begin(), ninr2.end());
    std::cout<<"Median number of neighbours within search radius "<<ninr2[N/2]<<std::endl;
    end_time();
    get_time_taken();
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " num_particles \n";
        return 1;
    }
#ifdef _OPENMP 
    std::cout<<"Running with OpenMP "<<_OPENMP<<std::endl;
    std::cout<<"Using "<<omp_get_max_threads()<<" threads and "<<omp_get_max_active_levels()<<" level "<<std::endl;
#endif
    // generate vector 
    std::size_t size = atoi(argv[1]);
    auto parts = generate_vector(size);

    auto treetypes = TreeTypes();

    for (auto &m:treetypes) {
        std::cout<<"Building and testing a "<<m.first<<std::endl;
        std::cout<<"=========================="<<std::endl;
        auto tree = build_kdtree(parts, 16, m.second);
        kdtree_test_NN(tree,parts);
        kdtree_test_ballsearch(tree,parts);
        delete tree;
        std::cout<<std::endl;
    }

}