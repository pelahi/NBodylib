/*! \file KDTree.cxx
 *  \brief This file contains subroutines involving tree construction

    NOTE: when possible, openmp parallelism implemented otherwise have MPI or simple serial code
    NOTE: openmp reductions are not implemented for min, max searches as not part of OpenMP API yet.
    and realistically, because one would have to implement a critical stop to check global minimum
    not really worth it. With MPI, obviously could sent bucket list to appropriate threads and min/max
    for that thread and do global broadcast to find global min/max but again overhead means it might not be
    worth it.
    I could also generate a array of size numthreads and then find min/max in that thread
    Must implement check that size of array worthwhile running in parallel
    \todo I can also reduce tree build time related to select the median point. A simple heuristic to avoid coding a complex linear-time median-finding algorithm, or using an O(n log n) sort of all n points, is to use sort to find the median of a fixed number of randomly selected points to serve as the splitting plane. In practice, this technique often results in nicely balanced trees.

*/

#include <KDTree.h>

namespace NBody
{

    #ifdef USEOPENMP
    // define openmp reduction operations
    // Currently issue is that it appears that most gcc compilers have not implemented omp delcare omp
    // clang 9 works but not gcc-9 with the typical gomp lib
    /*
    #pragma omp declare reduction( \
        vec_plus : \
        std::vector<Double_t> : \
        std::transform(omp_in.cbegin(), omp_in.cend(), omp_out.cbegin(), omp_out.begin(), \
        [](Double_t x, Double_t y)-> Double_t {return x+y;})) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
    #pragma omp declare reduction( \
        vec_max : \
        std::vector<Double_t> : \
        std::transform(omp_in.cbegin(), omp_in.cend(), omp_out.cbegin(), omp_out.begin(), \
        [](Double_t x, Double_t y)-> Double_t {return std::max<Double_t>(x,y);})) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
    #pragma omp declare reduction( \
        vec_min : \
        std::vector<Double_t> : \
        std::transform(omp_in.cbegin(), omp_in.cend(), omp_out.cbegin(), omp_out.begin(), \
        [](Double_t x, Double_t y)-> Double_t {return std::min<Double_t>(x,y);})) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
    */
    #endif

    /// \name wrappers for gettting the desired info from particles
    //@{
    void KDTree::GetParticlePos(const Particle &p, vector<Double_t> &x) {
        for (auto j = 0; j < ND; j++) x[j] = p.GetPosition(j);
    }
    void KDTree::GetParticleVel(const Particle &p, vector<Double_t> &x) {
        for (auto j = 0; j < ND; j++) x[j] = p.GetVelocity(j);
    }
    void KDTree::GetParticlePhs(const Particle &p, vector<Double_t> &x) {
        for (auto j = 0; j < ND; j++) x[j] = p.GetPhase(j);
    }
    void KDTree::GetParticleProj(const Particle &p, vector<Double_t> &x) {
        for (auto j = 0; j < ND; j++) x[j] = p.GetPosition(j);
    }

    Double_t KDTree::GetParticleithPos(const Particle &p, int i) {
        return p.GetPosition(i);
    }
    Double_t KDTree::GetParticleithVel(const Particle &p, int i) {
        return p.GetVelocity(i);
    }
    Double_t KDTree::GetParticleithPhs(const Particle &p, int i) {
        return p.GetPhase(i);
    }
    Double_t KDTree::GetParticleithProj(const Particle &p, int i) {
        return p.GetPosition(i);
    }

    /// \name Find most spread dimension
    //@{

    void KDTree::Spreadest(Int_t start, Int_t end, Double_t bnd[6][2], vector<Double_t> &spread,
        KDTreeOMPThreadPool &otp)
    {
        vector<Double_t> minval(ND), maxval(ND), x(ND);
        (this->*getparticlepos)(bucket[start], x);
        for (auto j=0;j<ND;j++) minval[j] = maxval[j] = x[j];
        // currently reduction of vectors using min, max seems to have issues with precision
        // look like automatic conversion to floats, which is odd. Leaving code here
        // for completeness but using explicit critical section.
        /*
#ifdef USEOPENMP
#pragma omp parallel for \
default(shared) schedule(static) \
firstprivate(x) \
reduction(vec_min:minval) reduction(vec_max:maxval) \
num_threads(nthreads) if (nthreads>1)
#endif
        for (auto i = start + 1; i < end; i++)
        {
            (this->*getparticlepos)(bucket[i], x);
            for (auto j = 0; j < ND; j++)
            {
                minval[j] = min(minval[j],x[j]);
                maxval[j] = max(maxval[j],x[j]);
            }
        }
        */
        // using critical sections and local min/max
#ifdef USEOPENMP
        unsigned int nthreads;
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
        Int_t delta = ceil((end-start)/(double)nthreads);
        if (nthreads>1) {
#pragma omp parallel \
default(shared) \
firstprivate(x)
{
        vector<Double_t> localminval(ND);
        vector<Double_t> localmaxval(ND);
        // since this is nested thread id doesn't simply map to how
        // the local for loop is split so construct a tid to index map
        unordered_map<int, int> tidtoindex;
        unsigned int tid, count=0;
        #pragma omp critical
        {
            tid = omp_get_thread_num();
            tidtoindex[tid] = count++;
        }
        // determine region of for loop to process
        tid = tidtoindex[omp_get_thread_num()];
        auto localstart = start + delta * static_cast<Int_t>(tid);
        auto localend = localstart + delta;
        if (tid == nthreads-1) localend = end;
        (this->*getparticlepos)(bucket[localstart], x);
        for (auto j=0;j<ND;j++) localminval[j] = localmaxval[j] = x[j];

        #pragma omp for
        for (auto i = localstart + 1; i < localend; i++)
        {
            (this->*getparticlepos)(bucket[i], x);
            for (auto j = 0; j < ND; j++)
            {
                localminval[j] = min(localminval[j],x[j]);
                localmaxval[j] = max(localmaxval[j],x[j]);
            }
        }

        #pragma omp critical
        {
//             tid = tidtoindex[omp_get_thread_num()];
            for (auto j = 0; j < ND; j++)
            {
                minval[j] = min(localminval[j],minval[j]);
                maxval[j] = max(localmaxval[j],maxval[j]);
            }
        }
}
        }
        else
#endif
        {
            for (auto i = start + 1; i < end; i++)
            {
                (this->*getparticlepos)(bucket[i], x);
                for (auto j = 0; j < ND; j++)
                {
                    minval[j] = min(minval[j],x[j]);
                    maxval[j] = max(maxval[j],x[j]);
                }
            }
        }

        // use min max to set boundaries and spread
        for (auto j = 0; j < ND; j++) {
            bnd[j][0] = minval[j];
            bnd[j][1] = maxval[j];
            spread[j] = maxval[j] - minval[j];
        }
    }
    //@}

    /// \name Find the boundary of the data and return mean
    //@{
    void KDTree::BoundaryandMean(Int_t start, Int_t end, Double_t bnd[6][2], vector<Double_t> &mean,
        KDTreeOMPThreadPool &otp)
    {
        vector<Double_t> minval(ND), maxval(ND), x(ND);
        (this->*getparticlepos)(bucket[start], x);
        for (auto j=0;j<ND;j++) {minval[j] = maxval[j] = x[j]; mean[j] = 0;}

        // currently reduction of vectors using min, max seems to have issues with precision
        // look like automatic conversion to floats, which is odd. Mean which is vector sum reduction seems fine.
        // Leaving code here for completeness but using explicit critical section
        /*
#ifdef USEOPENMP
        unsigned int nthreads;
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
#pragma omp parallel for \
default(shared) schedule(static) \
firstprivate(x) \
reduction(vec_min:minval) reduction(vec_max:maxval) reduction(vec_plus:mean) \
num_threads(nthreads) if (nthreads>1)
#endif
        for (auto i = start; i < end; i++)
        {
            (this->*getparticlepos)(bucket[i], x);
            for (auto j = 0; j < ND; j++)
            {
                minval[j] = min(minval[j],x[j]);
                maxval[j] = max(maxval[j],x[j]);
                mean[j] += x[j];
            }
        }
        */

        // using critical sections and local min/max and mean
#ifdef USEOPENMP
        unsigned int nthreads;
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
        Int_t delta = ceil((end-start)/(double)nthreads);
        if (nthreads>1) {
#pragma omp parallel \
default(shared) \
firstprivate(x)
{
        vector<Double_t> localminval(ND), localmaxval(ND), localmean(ND,0);
        // since this is nested thread id doesn't simply map to how
        // the local for loop is split so construct a tid to index map
        unordered_map<int, int> tidtoindex;
        unsigned int tid, count=0;
        #pragma omp critical
        {
            tid = omp_get_thread_num();
            tidtoindex[tid] = count++;
        }
        // determine region of for loop to process
        tid = tidtoindex[omp_get_thread_num()];
        auto localstart = start + delta * static_cast<Int_t>(tid);
        auto localend = localstart + delta;
        if (tid == nthreads-1) localend = end;
        (this->*getparticlepos)(bucket[localstart], x);
        for (auto j=0;j<ND;j++) localminval[j] = localmaxval[j] = localmean[j] = x[j];
        #pragma omp for
        for (auto i = localstart + 1; i < localend; i++)
        {
            (this->*getparticlepos)(bucket[i], x);
            for (auto j = 0; j < ND; j++)
            {
                localminval[j] = min(localminval[j],x[j]);
                localmaxval[j] = max(localmaxval[j],x[j]);
                localmean[j] += x[j];
            }
        }

        #pragma omp critical
        {
            for (auto j = 0; j < ND; j++)
            {
                minval[j] = min(localminval[j],minval[j]);
                maxval[j] = max(localmaxval[j],maxval[j]);
                mean[j] += localmean[j];
            }
        }
}
        }
        else
#endif
        {
            for (auto i = start; i < end; i++)
            {
                (this->*getparticlepos)(bucket[i], x);
                for (auto j = 0; j < ND; j++)
                {
                    minval[j] = min(minval[j],x[j]);
                    maxval[j] = max(maxval[j],x[j]);
                    mean[j] += x[j];
                }
            }
        }

        auto norm = 1.0/(Double_t)(end-start);
        for (auto j = 0; j < ND; j++) {
            bnd[j][0] = minval[j];
            bnd[j][1] = maxval[j];
            mean[j] *= norm;
        }

    }
    //@}

    /// \name Find the dispersion in a dimension (biased variance using 1/N as opposed to 1/(N-1) so that if N=2, doesn't crash)
    //@{
    void KDTree::Dispersion(Int_t start, Int_t end, vector<Double_t> &mean, vector<Double_t> &disp,
        KDTreeOMPThreadPool &otp)
    {
        vector<Double_t> x(ND);
        for (auto &d:disp) d=0;
// #ifdef USEOPENMP
//         unsigned int nthreads;
//         nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
//         if (nthreads <1) nthreads=1;
// #pragma omp parallel for \
// default(shared) schedule(static) \
// firstprivate(x) \
// reduction(vec_plus:disp) \
// num_threads(nthreads) if (nthreads>1)
// #endif
//         for (auto i = start; i < end; i++) {
//             (this->*getparticlepos)(bucket[i], x);
//             for (auto j = 0; j < ND; j++) disp[j]+=pow(x[j]-mean[j], static_cast<Double_t>(2.0));
//         }
#ifdef USEOPENMP
        unsigned int nthreads;
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
        Int_t delta = ceil((end-start)/(double)nthreads);
        if (nthreads>1) {
#pragma omp parallel \
default(shared) \
firstprivate(x)
{
        vector<Double_t> localdisp(ND,0);
        // since this is nested thread id doesn't simply map to how
        // the local for loop is split so construct a tid to index map
        unordered_map<int, int> tidtoindex;
        unsigned int tid, count=0;
        #pragma omp critical
        {
            tid = omp_get_thread_num();
            tidtoindex[tid] = count++;
        }
        // determine region of for loop to process
        tid = tidtoindex[omp_get_thread_num()];
        auto localstart = start + delta * static_cast<Int_t>(tid);
        auto localend = localstart + delta;
        if (tid == nthreads-1) localend = end;

        #pragma omp for
        for (auto i = localstart; i < localend; i++)
        {
            (this->*getparticlepos)(bucket[i], x);
            for (auto j = 0; j < ND; j++) localdisp[j]+=pow(x[j]-mean[j], static_cast<Double_t>(2.0));
        }

        #pragma omp critical
        {
            for (auto j = 0; j < ND; j++)
            {
                disp[j] += localdisp[j];
            }
        }
}
        }
        else
#endif
        {
            for (auto i = start; i < end; i++)
            {
                (this->*getparticlepos)(bucket[i], x);
                for (auto j = 0; j < ND; j++) disp[j]+=pow(x[j]-mean[j], static_cast<Double_t>(2.0));
            }
        }

        auto norm = 1.0/(Double_t)(end-start);
        for (auto j = 0; j < ND; j++) {
            disp[j] *= norm;
        }
    }
    //@}

    /// \name Calculate the entropy in a given dimension. This can be used as a node splitting criterion
    /// This calculates Shannon Entropy, where the region is split into nbins=pow(N,1/3) (N is number of particles) where minimum nbins=1,
    /// and can be used instead of most spread dimension
    //@{
    void KDTree::Entropy(Int_t start, Int_t end, Int_t nbins,
       Double_t bnd[6][2], vector<Double_t> &spread, vector<Double_t> &entropy,
       KDTreeOMPThreadPool &otp)
   {
       vector<Double_t> low(ND), up(ND), x(ND), dx(ND);
       auto norm = 2.0/(Double_t)(end-start);
       vector<Double_t> nientropy(nbins*ND);
       Double_t mtot = 0.;
       for (auto j=0;j<ND;j++)
       {
            low[j] = bnd[j][0] - spread[j]*norm;
            up[j] = bnd[j][1] + spread[j]*norm;
            dx[j] = (up[j] - low[j])/(Double_t)nbins;
       }
       for (auto &ni:nientropy) ni = 0.;
// #ifdef USEOPENMP
//        unsigned int nthreads;
//        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
//        if (nthreads <1) nthreads=1;
// #pragma omp parallel for \
// default(shared) schedule(static) \
// firstprivate(x) \
// reduction(+:mtot) reduction(vec_plus:nientropy) \
// num_threads(nthreads) if (nthreads>1)
// #endif
//        for (auto i=start;i<end;i++)
//        {
//            auto mass=bucket[i].GetMass();
//            (this->*getparticlepos)(bucket[i], x);
//            for (auto j=0;j<ND;j++) {
//                auto ibin=static_cast<Int_t>((x[j]-low[j])/dx[j]);
//                nientropy[ibin+j*nbins] += mass;
//            }
//            mtot += mass;
//        }

#ifdef USEOPENMP
        unsigned int nthreads;
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
        Int_t delta = ceil((end-start)/(double)nthreads);
        if (nthreads>1) {
#pragma omp parallel \
default(shared) \
firstprivate(x)
{
        vector<Double_t> localnientropy = nientropy;
        Double_t localmtot=0;
        // since this is nested thread id doesn't simply map to how
        // the local for loop is split so construct a tid to index map
        unordered_map<int, int> tidtoindex;
        unsigned int tid, count=0;
        #pragma omp critical
        {
            tid = omp_get_thread_num();
            tidtoindex[tid] = count++;
        }
        // determine region of for loop to process
        tid = tidtoindex[omp_get_thread_num()];
        auto localstart = start + delta * static_cast<Int_t>(tid);
        auto localend = localstart + delta;
        if (tid == nthreads-1) localend = end;
        (this->*getparticlepos)(bucket[localstart], x);
        #pragma omp for
        for (auto i = localstart; i < localend; i++)
        {
           auto mass=bucket[i].GetMass();
           (this->*getparticlepos)(bucket[i], x);
           for (auto j=0;j<ND;j++) {
               auto ibin=static_cast<Int_t>((x[j]-low[j])/dx[j]);
               localnientropy[ibin+j*nbins] += mass;
           }
           localmtot += mass;
        }

        #pragma omp critical
        {
            mtot += localmtot;
            for (auto j=0u;j<nientropy.size();j++) nientropy[j] += localnientropy[j];
        }
}
        }
        else
#endif
        {
            for (auto i=start;i<end;i++)
            {
                auto mass=bucket[i].GetMass();
                (this->*getparticlepos)(bucket[i], x);
                for (auto j=0;j<ND;j++) {
                    auto ibin=static_cast<Int_t>((x[j]-low[j])/dx[j]);
                    nientropy[ibin+j*nbins] += mass;
                }
                mtot += mass;
            }
        }

        mtot=1.0/mtot;
        norm = 1.0/log10((Double_t)nbins);
        for (auto j=0;j<ND;j++)
        {
            auto offset = j*nbins;
            for (auto i=0;i<nbins;i++)
            {
                if (nientropy[i] == 0) continue;
                auto temp=nientropy[i+offset]*mtot;
                entropy[j] -= temp*log10(temp);
            }
            entropy[j] *= norm;
        }
   }

    /// \name Determine the median coordinates in some space
    //@{
    Double_t KDTree::Median(int d, Int_t k, Int_t start, Int_t end,
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Int_t left = start;
        Int_t right = end-1;
        Int_t i, j;
        Double_t x;
        Particle w;
        Particle *pval = nullptr;
        //produced a balanced tree
        if (balanced){
            while (left < right)
            {
                x = (this->*getparticleithpos)(bucket[k],d);
                swap(bucket[right],bucket[k]);
                pval = &bucket[k];
                i = left-1;
                j = right;
                while (1) {
                    while (i < j) if ((this->*getparticleithpos)(bucket[++i],d) >= x) break;
                    while (i < j) if ((this->*getparticleithpos)(bucket[--j],d) <= x) break;
                    swap(bucket[i],bucket[j]);
                    pval = &bucket[j];
                    if (j <= i) break;
                }
                w = *pval;
                bucket[j] = move(bucket[i]);
                bucket[i] = move(bucket[right]);
                bucket[right] = w;
                pval = nullptr;
                if (i >= k) right = i - 1;
                if (i <= k) left = i + 1;
            }
            return (this->*getparticleithpos)(bucket[k], d);
        }
        //requires that particle order is already balanced. Use with caution
        else
        {
            return (this->*getparticleithpos)(bucket[k], d);
            //printf("Note yet implemented\n");
            //exit(9);
        }
    }
    //@}
    int KDTree::DetermineSplitDim(Int_t start, Int_t end, Double_t bnd[6][2],
        KDTreeOMPThreadPool &otp)
    {
        int splitdim=0;
        Double_t cursplitvalue;
        Double_t nbins;
        vector<Double_t> spread(ND), mean(ND), splitvalue(ND);
        // vector<Double_t> entropybins;

        //if using shannon entropy criterion
        if(splittingcriterion==1) {
            if(end-start>8) nbins=ceil(pow((end-start),1./3.));
            else nbins=2;
            (this->*spreadfunc)(start, end, bnd, spread, otp);
            (this->*entropyfunc)(start, end, nbins, bnd, spread, splitvalue, otp);
        }
        //if using dispersion
        else if (splittingcriterion==2) {
            (this->*bmfunc)(start, end, bnd, mean, otp);
            (this->*dispfunc)(start, end, mean, splitvalue, otp);
        }
        // if just using spread
        else {
            (this->*spreadfunc)(start, end, bnd, splitvalue, otp);
        }

        splitdim=0; cursplitvalue = splitvalue[0];
        //splitdim=0; maxspread=0.0; minentropy=1.0;enflag=0;
        //for since entropy can only be used in cases where the subspace is not sparse or does not have lattice structure must check
        //the question is how? At the moment, I do not check for this, though the idea would be only check dimensions that meet the criteria
        //and if non of them meet it, then enflag still zero and perhaps, use most spread dimension
        for (auto j = 1; j < ND; j++)
        {
            if (splitvalue[j]>cursplitvalue) {
                splitdim = j;
                cursplitvalue = splitvalue[j];
            }
        }
        return splitdim;
    }

    //-- End of inline functions

    //-- Private functions used to build tree

    /// Recursively build the nodes of the tree.  This works by first finding the dimension under
    /// which the data has the most spread, and then splitting the data about the median
    /// in that dimension.  BuildNodes() is then called on each half. Once the size of the data is
    /// small enough, a leaf node is formed.
    Node *KDTree::BuildNodes(Int_t start, Int_t end, KDTreeOMPThreadPool &otp)
    {
        Double_t bnd[6][2];
        Int_t size = end - start;
        Int_tree_t id = 0;
        //if not building in parallel can set ids here and update number of nodes
        //otherwise, must set after construction
        if (ibuildinparallel == false) {
            id = numnodes;
            numnodes++;
        }
        if (size <= b)
        {
            if (ibuildinparallel == false) numleafnodes++;
            vector<Double_t> mean(ND);
            (this->*bmfunc)(start, end, bnd, mean, otp);
            return new LeafNode(id, start, end, bnd, ND);
        }
        else
        {
            bool irearrangeandbalance=true;
            if (ikeepinputorder) irearrangeandbalance=false;
            Int_t k = start + (size - 1) / 2;
            int splitdim = DetermineSplitDim(start, end, bnd, otp);
            Double_t splitvalue = (this->*medianfunc)(splitdim, k, start, end, otp, irearrangeandbalance);
             //run the node construction in parallel
            if (ibuildinparallel && otp.nactivethreads > 1) {
                //note that if OpenMP not defined then ibuildinparallel is false
#ifdef USEOPENMP
                vector<KDTreeOMPThreadPool> newotp = OMPSplitThreadPool(otp);
                Node *left, *right;
                #pragma omp parallel default(shared) num_threads(2)
                #pragma omp single
                {
                    #pragma omp task
                    left = BuildNodes(start, k+1, newotp[0]);
                    #pragma omp task
                    right = BuildNodes(k+1, end, newotp[1]);
                    #pragma omp taskwait
                }
                return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    left, right);
#endif
            }
            else {
                return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    BuildNodes(start, k+1, otp), BuildNodes(k+1, end, otp));
            }
        }
    }

    ///scales the space and calculates the corrected volumes
    ///note here this is not mass weighted which may lead to issues later on.
    void KDTree::ScaleSpace(){

        if (treetype==TPHYS||treetype==TPROJ)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xmean[j]+=bucket[i].GetPosition(j);
        else if (treetype==TVEL)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xmean[j]+=bucket[i].GetVelocity(j);
        else if (treetype==TPHS||treetype==TMETRIC)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xmean[j]+=bucket[i].GetPhase(j);

        for(int j=0;j<ND;j++) xmean[j]/=(Double_t)numparts;

        if (treetype==TPHYS||treetype==TPROJ)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xvar[j]+=(bucket[i].GetPosition(j)-xmean[j])*(bucket[i].GetPosition(j)-xmean[j]);
        else if (treetype==TVEL)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xvar[j]+=(bucket[i].GetVelocity(j)-xmean[j])*(bucket[i].GetVelocity(j)-xmean[j]);
        else if (treetype==TPHS||treetype==TMETRIC)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xvar[j]+=(bucket[i].GetPhase(j)-xmean[j])*(bucket[i].GetPhase(j)-xmean[j]);

        for(int j=0;j<ND;j++){xvar[j]=sqrt(xvar[j]/(Double_t)numparts);ixvar[j]=1./xvar[j];}
        if (treetype==TPHYS||treetype==TPROJ)
            for (Int_t i=0;i<numparts;i++)
    	        for (int j=0;j<ND;j++) bucket[i].SetPosition(j,bucket[i].GetPosition(j)*ixvar[j]);
        else if (treetype==TVEL)
            for (Int_t i=0;i<numparts;i++)
    	        for (int j=0;j<ND;j++) bucket[i].SetVelocity(j,bucket[i].GetVelocity(j)*ixvar[j]);
        else if (treetype==TPHS||treetype==TMETRIC)
            for (Int_t i=0;i<numparts;i++)
    	        for (int j=0;j<ND;j++) bucket[i].SetPhase(j,bucket[i].GetPhase(j)*ixvar[j]);
    }

    int KDTree::TreeTypeCheck(){
        if (treetype<TPHYS||treetype>TMETRIC) {
            printf("Error in type of tree specified\n");
            printf("Returning null pointer as root node and leaving particle list unchanged.\n");
            root=nullptr; return 0;
        }
        else if (treetype==TMETRIC&&metric==nullptr) {
            printf("For treetype=%d must pass metric\n",TMETRIC);
            printf("Returning null pointer as root node and leaving particle list unchanged.\n");
            root=nullptr; return 0;
        }
        else {
            bmfunc=&NBody::KDTree::BoundaryandMean;
            dispfunc=&NBody::KDTree::Dispersion;
            spreadfunc=&NBody::KDTree::Spreadest;
            entropyfunc=&NBody::KDTree::Entropy;
            medianfunc=&NBody::KDTree::Median;
            if (treetype==TPHYS)
            {
                // bmfunc=&NBody::KDTree::BoundaryandMeanPos;
                // dispfunc=&NBody::KDTree::DispersionPos;
                // spreadfunc=&NBody::KDTree::SpreadestPos;
                // entropyfunc=&NBody::KDTree::EntropyPos;
                // medianfunc=&NBody::KDTree::MedianPos;
                getparticlepos=&NBody::KDTree::GetParticlePos;
                getparticleithpos=&NBody::KDTree::GetParticleithPos;
            }
            else if (treetype==TVEL) {
                getparticlepos=&NBody::KDTree::GetParticleVel;
                getparticleithpos=&NBody::KDTree::GetParticleithVel;
                // bmfunc=&NBody::KDTree::BoundaryandMeanVel;
                // dispfunc=&NBody::KDTree::DispersionVel;
                // spreadfunc=&NBody::KDTree::SpreadestVel;
                // entropyfunc=&NBody::KDTree::EntropyVel;
                // medianfunc=&NBody::KDTree::MedianVel;
            }
            else if (treetype==TPHS) {
                getparticlepos=&NBody::KDTree::GetParticlePhs;
                getparticleithpos=&NBody::KDTree::GetParticleithPhs;
                // bmfunc=&NBody::KDTree::BoundaryandMeanPhs;
                // dispfunc=&NBody::KDTree::DispersionPhs;
                // spreadfunc=&NBody::KDTree::SpreadestPhs;
                // entropyfunc=&NBody::KDTree::EntropyPhs;
                // medianfunc=&NBody::KDTree::MedianPhs;
            }
            else if (treetype==TPROJ) {
                getparticlepos=&NBody::KDTree::GetParticleProj;
                getparticleithpos=&NBody::KDTree::GetParticleithProj;
                // bmfunc=&NBody::KDTree::BoundaryandMeanPos;
                // dispfunc=&NBody::KDTree::DispersionPos;
                // spreadfunc=&NBody::KDTree::SpreadestPos;
                // entropyfunc=&NBody::KDTree::EntropyPos;
                // medianfunc=&NBody::KDTree::MedianPos;
            }
            return 1;
        }
    }
    ///Calculate kernel quantities
    void KDTree::KernelConstruction(){
        if (kernres<100) {
            printf("Error in kernel resolution, <100, setting to 100\n");
            kernres=100;
        }
        Kernel=new Double_t[kernres];
        derKernel=new Double_t[kernres];

        //determine dimension for kernel normalization
        if (treetype==TPHYS) ND=3;
        else if (treetype==TVEL) ND=3;
    	else if (treetype==TPHS) ND=6;
        else if (treetype==TPROJ) ND=2;
        else if (treetype==TMETRIC) ND=6;
        kernnorm=1.0/ND*pow(M_1_PI,ND/2.)*gsl_sf_gamma(ND/2.+1.);

        //smoothing kernel type
        if (kernfunctype==KSPH) {
            kernfunc.function=WSPH;
            kernnorm*=ND*(ND+1.)*(ND+2.)*(ND+3.)/(6*(pow(2.,ND+1)-1.));
        }
        else if (kernfunctype==KGAUSS) {
            kernfunc.function=WGauss;
            kernnorm*=pow(0.5*M_1_PI,ND/2.);
        }
        else if (kernfunctype==KEPAN) {
            kernfunc.function=WEpan;
            kernnorm*=ND*(ND+2.)*pow(0.5,ND+1.);
        }
        else if (kernfunctype==KTH) kernfunc.function=WTH;
        else {
            printf("Error in kernel choice, using default SPH B-spline kernel\n");
            kernfunc.function=WSPH;;
            kernnorm*=ND*(ND+1.)*(ND+2.)*(ND+3.)/(6*(pow(2.,ND+1)-1.));
        }
        Double_t delta=2.0/(Double_t)(kernres-1);
        for (int i=0;i<kernres;i++) {
            Kernel[i]=kernnorm*kernfunc.function(i*delta,1.0);
        }
    }

    void KDTree::BuildNodeIDs()
    {
        numnodes = 0;
        numleafnodes = 0;
        UpdateNodeID(root);
    }

    void KDTree::UpdateNodeID(Node *node)
    {
        node->SetID(numnodes++);
        //walk tree increasing
        if (node->GetCount() <= b) {
            numleafnodes++;
            return;
        }
        else {
            UpdateNodeID(((SplitNode*)node)->GetLeft());
            UpdateNodeID(((SplitNode*)node)->GetRight());
        }
    }

    void KDTree::WalkNodesFromRoot()
    {
        cout<<"Walking tree"<<endl;
        WalkNode(root);
    }

    void KDTree::WalkNode(Node *node)
    {
        Int_t start, end;
        Int_tree_t id;
        start = node->GetStart();
        end = node->GetEnd();
        id = node->GetID();
        cout<<"At node "<<" "<<id<<" "<<start<<" "<<end<<" ";
        for (auto j=0;j<ND;j++)  cout<<"("<<node->GetBoundary(j,0)<<", "<<node->GetBoundary(j,1)<<")";
        cout<<endl;
        if (node->GetCount() > b) {
            WalkNode(((SplitNode*)node)->GetLeft());
            WalkNode(((SplitNode*)node)->GetRight());
        }
    }

    //-- End of private functions used to build the tree

    //-- Public constructors

    KDTree::KDTree(Particle *p, Int_t nparts, Int_t bucket_size,
      int ttype, int smfunctype, int smres,
      int criterion, int aniso, int scale,
      Double_t *Period, Double_t **m,
      bool iBuildInParallel,
      bool iKeepInputOrder)
    {
        iresetorder=true;
        ikeepinputorder = iKeepInputOrder;
        ibuildinparallel = false;
#ifdef USEOPENMP
        ibuildinparallel = iBuildInParallel;
        bool inested = omp_get_nested();
        int nthreads;
        #pragma omp parallel
        #pragma omp single
        {
            nthreads = omp_get_num_threads();
        }
        if (nthreads == 1) ibuildinparallel = false;
        if (inested == false) omp_set_nested(int(ibuildinparallel));
#endif
        numparts = nparts;
        numleafnodes=numnodes=0;
        bucket = p;
        b = bucket_size;
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (Period!=nullptr)
        {
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=Period[k];
        }
        else period=nullptr;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            // if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            KDTreeOMPThreadPool otp = OMPInitThreadPool();
            root=BuildNodes(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
            // if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_nested(inested);
#endif
    }

    KDTree::KDTree(System &s, Int_t bucket_size,
      int ttype, int smfunctype, int smres, int criterion, int aniso, int scale,
      Double_t **m,
      bool iBuildInParallel,
      bool iKeepInputOrder
    )
    {

        iresetorder=true;
        ikeepinputorder = iKeepInputOrder;
        ibuildinparallel = false;
#ifdef USEOPENMP
        ibuildinparallel = iBuildInParallel;
        bool inested = omp_get_nested();
        int nthreads;
        #pragma omp parallel
        #pragma omp single
        {
            nthreads = omp_get_num_threads();
        }
        if (nthreads == 1) ibuildinparallel = false;
        if (inested == false) omp_set_nested(int(ibuildinparallel));
#endif
        numparts = s.GetNumParts();
        numleafnodes=numnodes=0;
        bucket = s.Parts();
        b = bucket_size;
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (s.GetPeriod()[0]>0&&s.GetPeriod()[1]>0&&s.GetPeriod()[2]>0){
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=s.GetPeriod()[k];
        }
        else period=nullptr;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            // if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            KDTreeOMPThreadPool otp = OMPInitThreadPool();
            root=BuildNodes(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            // if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_nested(inested);
#endif
    }
    KDTree::~KDTree()
    {
	    if (root!=nullptr) {
            delete root;
            delete[] Kernel;
            delete[] derKernel;
            if (period!=NULL) delete[] period;
            if (iresetorder) std::sort(bucket, bucket + numparts, IDCompareVec);
            if (scalespace) {
                for (auto i=0;i<numparts;i++)
                    for (auto j=0;j<6;j++) bucket[i].SetPhase(j, bucket[i].GetPhase(j)*xvar[j]);
            }
        }
    }

    void KDTree::OverWriteInputOrder() {
        iresetorder=false;
        for (Int_t i=0;i<numparts;i++) bucket[i].SetID(i);
    }
    void KDTree::SetResetOrder(bool a) {iresetorder=a;}

    KDTreeOMPThreadPool KDTree::OMPInitThreadPool()
    {
        KDTreeOMPThreadPool ompthreadpool;
#ifdef USEOPENMP
        if (ibuildinparallel) {
            int nthreads;
            #pragma omp parallel
            #pragma omp single
            {
                nthreads = omp_get_num_threads();
            }
            ompthreadpool.nthreads = nthreads;
        }
        else {
            ompthreadpool.nthreads = 1;
        }
        ompthreadpool.nactivethreads = ompthreadpool.nthreads;
        ompthreadpool.activethreadids.resize(ompthreadpool.nactivethreads);
        for (auto i=0u;i<ompthreadpool.nactivethreads;i++) ompthreadpool.activethreadids[i]=i;
#endif
        return ompthreadpool;
    }
    vector<KDTreeOMPThreadPool> KDTree::OMPSplitThreadPool(KDTreeOMPThreadPool &ompthreadpool)
    {
        vector<KDTreeOMPThreadPool> newthreadpool(2);
#ifdef USEOPENMP
        if (ompthreadpool.nactivethreads >= 2)
        {
            newthreadpool[0].nthreads = ompthreadpool.nthreads;
            newthreadpool[1].nthreads = ompthreadpool.nthreads;
            newthreadpool[0].nactivethreads = ompthreadpool.nactivethreads/2;
            newthreadpool[1].nactivethreads = ompthreadpool.nactivethreads - newthreadpool[0].nactivethreads;
            newthreadpool[0].activethreadids.resize(newthreadpool[0].nactivethreads);
            newthreadpool[1].activethreadids.resize(newthreadpool[1].nactivethreads);
            for (auto i=0u;i<newthreadpool[0].nactivethreads;i++)
                newthreadpool[0].activethreadids[i] = ompthreadpool.activethreadids[i];
            for (auto i=0u;i<newthreadpool[1].nactivethreads;i++)
                newthreadpool[1].activethreadids[i] = ompthreadpool.activethreadids[i+newthreadpool[0].nactivethreads];
        }
#endif
        return newthreadpool;
    }
}
