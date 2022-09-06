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
#include <random>
#define XSORT

namespace NBody
{

    /// \defgroup GetPartDataForTree 
    /// Get particle data relevant to tree 
    //@{
    #pragma omp declare simd
    DoublePos_t KDTree::get_particle_pos_jth(Int_t i, int j) 
    {
        return bucket[i].GetPosition(j);
    }
    #pragma omp declare simd
    DoublePos_t KDTree::get_particle_vel_jth(Int_t i, int j) {
        return bucket[i].GetVelocity(j);
    }
    #pragma omp declare simd
    DoublePos_t KDTree::get_particle_phs_jth(Int_t i, int j) {
        return bucket[i].GetPhase(j);
    }
    std::vector<DoublePos_t> KDTree::get_particle_pos(Int_t i) 
    {
        std::vector<DoublePos_t> v(3);
        v[0]=bucket[i].X();
        v[1]=bucket[i].Y();
        v[2]=bucket[i].Z();
        return v;
    }
    std::vector<DoublePos_t> KDTree::get_particle_vel(Int_t i) 
    {
        std::vector<DoublePos_t> v(3);
        v[0]=bucket[i].Vx();
        v[1]=bucket[i].Vy();
        v[2]=bucket[i].Vz();
        return v;
    }
    std::vector<DoublePos_t> KDTree::get_particle_phs(Int_t i) 
    {
        std::vector<DoublePos_t> v(6);
        v[0]=bucket[i].X();
        v[1]=bucket[i].Y();
        v[2]=bucket[i].Z();
        v[3]=bucket[i].Vx();
        v[4]=bucket[i].Vy();
        v[5]=bucket[i].Vz();
        return v;
    }
    //@}

    /// \name Find most spread dimension
    //@{
    Double_t KDTree::Spreadest(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t minval = (this->*get_part_data_jth)(start, j);
        Double_t maxval = minval;
        Int_t i;
        unsigned int nthreads;
#ifdef USEOPENMP
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
#pragma omp parallel for \
default(shared) private(i) schedule(static) \
reduction(min:minval) reduction(max:maxval) num_threads(nthreads) if (nthreads>1)
#endif
        for (i = start + 1; i < end; i++)
        {
            auto x = (this->*get_part_data_jth)(i, j);
            if (x < minval) minval = x;
            if (x > maxval) maxval = x;
        }
        bnd[0]=minval;bnd[1]=maxval;
        return maxval - minval;
    }
    //@}
    /// \name Find the boundary of the data and return mean
    //@{
    inline Double_t KDTree::BoundaryandMean(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t mean=(this->*get_part_data_jth)(start, j);
        Double_t minval=mean, maxval=mean;
        Int_t i;
        unsigned int nthreads;
#ifdef USEOPENMP
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
#pragma omp parallel for \
default(shared) private(i) schedule(static) \
reduction(+:mean) reduction(min:minval) reduction(max:maxval) num_threads(nthreads) if (nthreads>1)
#endif
        for (i = start + 1; i < end; i++)
        {
            auto x = (this->*get_part_data_jth)(i, j);
            if (x < minval) minval = x;
            if (x > maxval) maxval = x;
            mean+=x;
        }
        bnd[0]=minval;bnd[1]=maxval;
        mean/=(Double_t)(end-start);
        return mean;
    }
    //@}
    /// \name Find the dispersion in a dimension (biased variance using 1/N as opposed to 1/(N-1) so that if N=2, doesn't crash)
    //@{
    inline Double_t KDTree::Dispersion(int j, Int_t start, Int_t end, Double_t mean,
        KDTreeOMPThreadPool &otp)
    {
        Double_t disp=0;
        Int_t i;
        unsigned int nthreads;
#ifdef USEOPENMP
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
#pragma omp parallel for \
default(shared) private(i) schedule(static) \
reduction(+:disp) num_threads(nthreads) if (nthreads>1)
#endif
        for (i = start; i < end; i++) 
        {
            auto x = (this->*get_part_data_jth)(i, j);
            disp+=(x-mean)*(x-mean);
        }
        disp/=(Double_t)(end-start);
        return disp;
    }
    //@}

/*
    //code that applies a correction to the boundary of a node, ibnd is initial boundary range estimate,
    //xbnd is the current parent nodes estimate
    //at the moment the code is not setup to correct for underestimation in outer and inner parts
    //(I'm not exactly sure how this corrects this, must check Enbid paper)
    inline void BoundaryCor(int j, Int_t count, Int_t dim, Int_t numparts, Double_t *ibnd, Double_t *xbnd){
        //factors in the correction (taken from Enbind where don't assume cubic cells)
        Double_t fac1=10.0, fac2=0.2, fac3=2.0;
        Double_t temp1,temp2,temp3;
        temp1=fac1*pow((Double_t)numparts,1.0/(Double_t)dim);
        temp2=fac2/temp1;
        if((ibnd[1]-ibnd[0])/(xbnd[1]-xbnd[0])> 1.0/(Double_t)numparts)
        {
            temp3=(ibnd[1]-ibnd[0])/(count-1);
            if(((xbnd[1]-ibnd[1])>(temp2*temp3))&&((ibnd[0]-xbnd[0])>(temp2*temp3)))
            {
                xbnd[1]=ibnd[1]+fac3*temp3;
                xbnd[0]=ibnd[0]-fac3*temp3;
            }
            else
            {
                if((xbnd[1]-ibnd[1])>(temp1*temp3)) xbnd[1]=ibnd[1]+fac3*temp3;
                if((ibnd[0]-xbnd[0])>(temp1*temp3)) xbnd[0]=ibnd[0]-fac3*temp3;
            }
        }
    }
*/
     /// \name Calculate the entropy in a given dimension. This can be used as a node splitting criterion
     /// This calculates Shannon Entropy, where the region is split into nbins=pow(N,1/3) (N is number of particles) where minimum nbins=1,
     /// and can be used instead of most spread dimension
     //@{
    Double_t KDTree::Entropy(int j, Int_t start, Int_t end,
        Double_t low, Double_t up, Double_t nbins, Double_t *nientropy,
        KDTreeOMPThreadPool &otp)
    {
        Int_t ibin, i;
        Double_t mtot=0.,entropy=0.;
        Double_t idx=static_cast<Double_t>(nbins)/(up-low);
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            auto x = (this->*get_part_data_jth)(i, j);
            ibin=static_cast<Int_t>(floor((x-low)*idx));
            nientropy[ibin]+=bucket[i].GetMass();
        }
        mtot=1.0/mtot;
        for (i=0;i<nbins;i++) {
            if (nientropy[i]>0) {
                Double_t temp=nientropy[i]*mtot;
                entropy-=temp*log10(temp);
            }
        }
        return entropy/log10((Double_t)nbins);
    }
    //@}

    /// \name Determine the median coordinates in some space
    //@{
    Double_t KDTree::Median(int d, Int_t &k, Int_t start, Int_t end, Double_t farthest, 
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Int_t left = start;
        Int_t right = end-1;
        Int_t i, j;
        Double_t x;
        Particle w;
        Particle *pval = NULL;
        //produced a balanced tree
        if (balanced){
            while (left < right)
            {
                x = (this->*get_part_data_jth)(k, d);
                swap(bucket[right],bucket[k]);
                pval = &bucket[k];
                i = left-1;
                j = right;
                while (1) {
                    while (i < j) 
                    {
                        auto xx = (this->*get_part_data_jth)(++i, d);
                        if (xx >= x) break;
                    }
                    while (i < j) 
                    {
                        auto xx = (this->*get_part_data_jth)(--j, d);
                        if (xx <= x) break;
                    }
                    swap(bucket[i],bucket[j]);
                    pval = &bucket[j];
                    if (j <= i) break;
                }
                w = *pval;
                bucket[j] = move(bucket[i]);
                bucket[i] = move(bucket[right]);
                bucket[right] = w;
                pval = NULL;
                if (i >= k) right = i - 1;
                if (i <= k) left = i + 1;
            }
            return (this->*get_part_data_jth)(k, d);
        }
        //requires that particle order is already balanced. Use with caution
        else
        {
            return (this->*get_part_data_jth)(k, d);
            //printf("Note yet implemented\n");
            //exit(9);
        }
    }
    //@}

    int KDTree::DetermineSplitDim(Int_t start, Int_t end, Double_t bnd[6][2],
        KDTreeOMPThreadPool &otp) {
        int splitdim=0;
        Double_t cursplitvalue;
        Double_t nbins;
        vector<Double_t> splitvalue(ND);
        vector<Double_t> entropybins;

        //if using shannon entropy criterion
        if(splittingcriterion==1) {
            if(end-start>8) nbins=ceil(pow((end-start),1./3.));
            else nbins=2;
            entropybins.resize(nbins+1);
        }
        for (auto j = 0; j < ND; j++)
        {
            if(splittingcriterion == KDTREE_SPLIT_ENTROPY) {
                splitvalue[j] = (this->*spreadfunc)(j, start, end, bnd[j], otp)+1e-32;//addition incase lattice and no spread
                Double_t low, up;
                low=bnd[j][0] - 2.0*(splitvalue[j])/(Double_t)(end-start);
                up=bnd[j][1] + 2.0*(splitvalue[j])/(Double_t)(end-start);
                splitvalue[j] = (this->*entropyfunc)(j, start, end, low, up, nbins, entropybins.data(), otp);
            }
            else if (splittingcriterion == KDTREE_SPLIT_DISPERSION) {
                splitvalue[j] = (this->*bmfunc)(j, start, end, bnd[j], otp);
                splitvalue[j] = (this->*dispfunc)(j, start, end, splitvalue[j], otp);
            }
            else {
                splitvalue[j] = (this->*spreadfunc)(j, start, end, bnd[j], otp);
            }
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


    //@{
    ///Determine whether to calculate inter-particle spacing and find maximum over just
    ///using the median as split index. 
    bool KDTree::UseMedianOverMaxInterparticleSpacing(Double_t nodefarthest, Int_t nodesize, Int_t bufferwidth) 
    {
        // if node is smaller than a * search distance, no need to adjust median splitting to maximum inter-particle
        // spacing splitting. 
        if (nodefarthest < rdist2daptwithfac) return true;
        // if node contains too few particles in which to search for max interparticle spacing 
        // no need to adjust median splitting to maximum inter-particle spacing 
        else if (bufferwidth<minadaptivemedianregionsize) return true;
        // check if nodesize is too small
        else if (nodesize < 2*b) return true;
        else return false;
    }

    ///adjust the sorted particle array so that the split is not at the median
    ///but approximative and splits at the particle with the largest distance
    ///between particles. This search is limited to a buffer region around
    ///the median index
    Double_t KDTree::AdjustMedianToMaximalDistance(int d,
        Int_t &splitindex, Int_t trueleft, Int_t trueright, Double_t farthest, 
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        UInt_tree_t truesize = trueright - trueleft;
        UInt_tree_t bufferwidth = truesize * adaptivemedianfac;
        Double_t splitvalue;
        // determine whether need to use Median over max interparcile spacing 
        if (UseMedianOverMaxInterparticleSpacing(farthest, truesize, bufferwidth)) {
            splitvalue = Median(d, splitindex, trueleft, trueright, farthest, otp, balanced);
            return splitvalue;
        }

        //now begin search for more optimal split point with large interparticle distance splitting. 
        UInt_tree_t left = splitindex - bufferwidth/2;
        UInt_tree_t right = splitindex + bufferwidth/2;
        UInt_tree_t size = right - left;
        splitvalue = Median(d, splitindex, trueleft, trueright, farthest, otp, balanced);
        //wonder if I should just sort particles (and whether to sort all particles)
        vector<KDTreeForSorting> x(size);
        for (auto i=0;i<size;i++) {
            x[i].val =  (this->*get_part_data_jth)(left+i, d);
            x[i].orgindex = left+i;
        }
        std::sort(x.begin(), x.end() , [](const KDTreeForSorting &a, const KDTreeForSorting &b) {
            return a.val < b.val;
        });

        // look at region around median and find that with maximum interparticle distance
        Double_t maxdist=0;
        UInt_tree_t maxi=0;
        Double_t startdist = x[size/2+1].val - x[size/2].val;
#ifdef USEOPENMP
        unsigned int nthreads;
        nthreads = min((unsigned int)(floor((size)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads < 1) nthreads=1;
        if (nthreads > 1) {
#pragma omp parallel \
default(shared) 
{
            Double_t localmaxdist;
            UInt_tree_t localmaxi;
            auto dist = (x[1].val - x[0].val);
            localmaxdist = dist;
            localmaxi = 0;
            #pragma omp for schedule(static) nowait
            for (UInt_tree_t i=1; i<size-1; i++)
            {
                dist = (x[i+1].val - x[i].val);
                if (dist > localmaxdist)
                {
                    localmaxdist = dist;
                    localmaxi = i;
                }
            }
            #pragma omp critical
            {
                if (localmaxdist > maxdist) {
                    maxdist = localmaxdist;
                    maxi = localmaxi;
                }
            }
}
        }
        else 
#endif 
        {
            auto dist = (x[1].val - x[0].val);
            for (UInt_tree_t i=1; i<size-1; i++)
            {
                dist = (x[i+1].val - x[i].val);
                if (dist > maxdist)
                {
                    maxdist = dist;
                    maxi = i;
                }
            }
        }

        splitindex = x[maxi].orgindex;
        splitvalue = x[maxi].val;
        // once split index found move particles. Only necessary if sort of particles not done 
        splitvalue = Median(d, splitindex, trueleft, trueright, farthest, otp, balanced);
        return splitvalue;
    }


    ///Calculate center and largest sqaured distance for node
    vector<Double_t> KDTree::DetermineCentreAndSmallestSphere(
        UInt_tree_t localstart, UInt_tree_t localend,
        Double_t &farthest,
        KDTreeOMPThreadPool &otp
        )
    {
        // get center
        vector<Double_t> center(ND,0), pos(ND);
        Double_t norm = 1.0/(static_cast<Double_t>(localend-localstart));
        Double_t maxr2 = 0.0;
        unsigned int nthreads = 1;
#ifdef USEOPENMP
        nthreads = min((unsigned int)(floor((localend - localstart)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
        UInt_tree_t delta = ceil((localend - localstart)/(double)nthreads);
        unordered_map<int, int> tidtoindex;
        vector<UInt_tree_t> threadlocalstart, threadlocalend;
#endif
        if (nthreads>1) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
            // since this is nested thread id doesn't simply map to how
            // the local for loop is split so construct a tid to index map
            int tid, count=0;
            #pragma omp critical
            {
                tid = omp_get_thread_num();
                tidtoindex[tid] = count++;
            }
            // determine region of for loop to process
            tid = tidtoindex[omp_get_thread_num()];
            threadlocalstart[tid] = localstart + delta * static_cast<Int_t>(tid);
            threadlocalend[tid] = threadlocalstart[tid] + delta;
            if (tid == nthreads-1) threadlocalend[tid] = localend;
            vector<Double_t> localcenter(ND,0);
            #pragma omp for nowait
            for (auto i = threadlocalstart[tid]; i < threadlocalend[tid]; i++)
            {
                for(auto j=0;j<ND;j++) localcenter[j] +=  (this->*get_part_data_jth)(i, j);
            }

            #pragma omp critical
            {
                for (auto j = 0; j < ND; j++) center[j] += localcenter[j];
            }
}
#endif
        }
        else
        {
            for(auto i=localstart; i<localend;i++)
            {
                for(auto j=0;j<ND;j++) center[j] += (this->*get_part_data_jth)(i, j);
            }
        }
        for (auto &c:center) c*= norm;

        //now find most distant particle
        if (nthreads>1) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
            int tid, count=0;
            tid = tidtoindex[omp_get_thread_num()];
            Double_t localmaxr2 = 0 ;
            #pragma omp for nowait
            for (auto i = threadlocalstart[tid]; i < threadlocalend[tid]; i++)
            {
                for(auto j=0;j<ND;j++) pos[j] = (this->*get_part_data_jth)(i, j);
                Double_t r2=0;
                for(auto j=0;j<ND;j++) r2+=(pos[j] - center[j])*(pos[j] - center[j]);
                localmaxr2 = std::max(localmaxr2, r2);
            }
            #pragma omp critical
            {
                maxr2 = std::max(maxr2, localmaxr2);
            }
}
#endif
        }
        else
        {
            //get largest distance
            for(auto i=localstart; i<localend;i++)
            {
                for(auto j=0;j<ND;j++) pos[j] = (this->*get_part_data_jth)(i, j);
                Double_t r2=0;
                for(auto j=0;j<ND;j++) r2+=(pos[j] - center[j])*(pos[j] - center[j]);
                maxr2 = std::max(maxr2, r2);
            }
        }
        farthest = maxr2;

        return center;
    }

    void KDTree::DetermineCentreAndSmallestSphere(
        UInt_tree_t localstart, UInt_tree_t localend,
        Node *&node,
        KDTreeOMPThreadPool &otp
    )
    {
        Double_t maxr2;
        vector<Double_t> center = DetermineCentreAndSmallestSphere(localstart, localend, maxr2, otp);
        for(auto j=0;j<ND;j++) node->SetCenter(j, center[j]);
        node->SetFarthest(maxr2);
    }


    ///Calculate maximum squared interparticle distance
    Double_t KDTree::DetermineMaxInterParticleSpacing(UInt_tree_t localstart, UInt_tree_t localend,
        int splitdim,
        KDTreeOMPThreadPool &otp
        )
    {
        // store position in desired dimension to find where to split particles
        Double_t maxinterdist = 0.0;
        UInt_tree_t size = (localend - localstart);
        vector<KDTreeForSorting> x(size);
        for (auto i=0; i<size; i++) {
            x[i].val = (this->*get_part_data_jth)(i + localstart, splitdim); 
            x[i].orgindex = i + localstart;
        }
        std::sort(x.begin(), x.end() , [](const KDTreeForSorting &a, const KDTreeForSorting &b) {
            return a.val < b.val;
        });
        unsigned int nthreads = 1;
#ifdef USEOPENMP
        nthreads = min(static_cast<unsigned int>(floor((size)/static_cast<float>(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
        UInt_tree_t delta = ceil((size)/(double)nthreads);
        unordered_map<int, int> tidtoindex;
        vector<UInt_tree_t> threadlocalstart, threadlocalend;
#endif
        if (nthreads>1) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
            // since this is nested thread id doesn't simply map to how
            // the local for loop is split so construct a tid to index map
            int tid, count=0;
            #pragma omp critical
            {
                tid = omp_get_thread_num();
                tidtoindex[tid] = count++;
            }
            // determine region of for loop to process
            tid = tidtoindex[omp_get_thread_num()];
            threadlocalstart[tid] = delta * static_cast<Int_t>(tid);
            threadlocalend[tid] = threadlocalstart[tid] + delta;
            if (tid == nthreads-1) threadlocalend[tid] = size;
            Double_t localmax = 0;
            #pragma omp for nowait
            for (auto i = threadlocalstart[tid]; i < threadlocalend[tid]-1; i++)
            {
                auto diff = (x[i+1].val - x[i].val);
                localmax = std::max(diff*diff,localmax);  
            }
            #pragma omp critical
            {
                maxinterdist = std::max(maxinterdist, localmax);
            }
}
#endif
        }
        else
        {
            for(auto i=0; i<size-1;i++)
            {
                auto diff = (x[i+1].val - x[i].val);
                maxinterdist = std::max(diff*diff,maxinterdist);
            }
        }
        return maxinterdist;
    }

    //@}

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
        int splitdim = -1;

        //if not building in parallel can set ids here and update number of nodes
        //otherwise, must set after construction
        if (ibuildinparallel == false) {
            id = numnodes;
            numnodes++;
        }
        bool isleafflag;
        vector<Double_t> center;
        Double_t localfarthest = 0, maxinterdist = 0;
        // if constructing adaptive tree where leaf nodes must be smaller than some size
        // calculate the farthest distance to the centre of the node
        if (rdist2adapt > 0) {
            center = DetermineCentreAndSmallestSphere(start, end, localfarthest, otp);
            isleafflag = ((size <= b && localfarthest < rdist2adapt) || (size <= bmin));
        }
        else {
            isleafflag = (size <= b);
        }
        if (isleafflag)
        {
            if (ibuildinparallel == false) numleafnodes++;
            for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j], otp);
            Node * leaf = new LeafNode(id, start, end,  bnd, ND, treetype);
            if (rdist2adapt > 0)
            {
                leaf->SetFarthest(localfarthest);
                for (int j=0;j<ND;j++) leaf->SetCenter(j,center[j]);
            }
            return leaf;
        }
        else
        {
            bool irearrangeandbalance=true;
            if (ikeepinputorder) irearrangeandbalance=false;
            if (splitdim == -1) splitdim = DetermineSplitDim(start, end, bnd, otp);
            Int_t splitindex = start + (size - 1) / 2;
            Double_t splitvalue = (this->*medianfunc)(splitdim, splitindex, start, end, localfarthest, otp, irearrangeandbalance);
             //run the node construction in parallel
            if (ibuildinparallel && otp.nactivethreads > 1) {
                //note that if OpenMP not defined then ibuildinparallel is false
                Node *left, *right;
#ifdef USEOPENMP
                vector<KDTreeOMPThreadPool> newotp = OMPSplitThreadPool(otp);
                #pragma omp parallel default(shared) num_threads(2)
                #pragma omp single
                {
                    #pragma omp task 
                    {
                        left = BuildNodes(start, splitindex+1, newotp[0]);
                        if (rdist2adapt>0) DetermineCentreAndSmallestSphere(start, splitindex+1, left, newotp[0]);
                    }
                    #pragma omp task
                    {
                        right = BuildNodes(splitindex+1, end, newotp[1]);
                        if (rdist2adapt>0) DetermineCentreAndSmallestSphere(splitindex+1, end, right, newotp[1]);
                    }
                    #pragma omp taskwait
                }
#endif
                return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    left, right, treetype);
            }
            else {
                if (rdist2adapt > 0) {
                    Node *left, *right;
                    left = BuildNodes(start, splitindex+1, otp);
                    DetermineCentreAndSmallestSphere(start, splitindex+1, left, otp);
                    right = BuildNodes(splitindex+1, end, otp);
                    DetermineCentreAndSmallestSphere(splitindex+1, end, right, otp);
                    return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    left, right, treetype);
                }
                else {
                    return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                        BuildNodes(start, splitindex+1, otp), BuildNodes(splitindex+1, end, otp));
                }
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
            root=NULL; return 0;
        }
        else if (treetype==TMETRIC&&metric==NULL) {
            printf("For treetype=%d must pass metric\n",TMETRIC);
            printf("Returning null pointer as root node and leaving particle list unchanged.\n");
            root=NULL; return 0;
        }
        else {
            bmfunc=&NBody::KDTree::BoundaryandMean;
            dispfunc=&NBody::KDTree::Dispersion;
            spreadfunc=&NBody::KDTree::Spreadest;
            entropyfunc=&NBody::KDTree::Entropy;
            if (adaptivemedianfac > 0) medianfunc=&NBody::KDTree::AdjustMedianToMaximalDistance;
            else medianfunc=&NBody::KDTree::Median;
            if (treetype==TPHYS)
            {
                get_part_data_jth=&NBody::KDTree::get_particle_pos_jth;
                get_part_data=&NBody::KDTree::get_particle_pos;
            }
            else if (treetype==TVEL) {
                get_part_data_jth=&NBody::KDTree::get_particle_vel_jth;
                get_part_data=&NBody::KDTree::get_particle_vel;
            }
            else if (treetype==TPHS) {
                get_part_data_jth=&NBody::KDTree::get_particle_phs_jth;
                get_part_data=&NBody::KDTree::get_particle_phs;
            }
            else if (treetype==TPROJ) {
                get_part_data_jth=&NBody::KDTree::get_particle_pos_jth;
                get_part_data=&NBody::KDTree::get_particle_pos;
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

    void KDTree::PeriodReflectionConstruction(){
        periodicreflectionindices.resize(ND);
        std::vector<int> indices(ND);
        for (auto i=0;i<ND;i++) indices[i] = i;
        for (auto i=1;i<=ND;i++) 
        {
            do {
                std::vector<int> mysub(indices.begin(), indices.begin()+i);
                std::sort(mysub.begin(), mysub.end());
                if (mysub.front() <= mysub.back()) periodicreflectionindices[i-1].insert(mysub);
            } while ( std::next_permutation(indices.begin(), indices.end()));
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
	    if(node->GetLeaf()){
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
        auto size = end - start;
        id = node->GetID();
        cout<<"At node "<<" "<<id<<" "<<node->GetLeaf()<<" "<<start<<" "<<end<<" "<<size<<" : ";
        if (!node->GetLeaf()) cout<<((SplitNode*)node)->GetCutDim()<<" "<<((SplitNode*)node)->GetCutValue()<<" : ";
        for (auto j=0;j<ND;j++)  cout<<"("<<node->GetBoundary(j,0)<<", "<<node->GetBoundary(j,1)<<"), ";
        auto dist=0.0; 
        for (auto j=0;j<ND;j++) dist+=pow(node->GetBoundary(j,0) - node->GetBoundary(j,1), 2.0);
        cout<<" : "<<sqrt(dist);
        cout<<" : "<<node->GetFarthest();
        cout<<endl;
	    if(!node->GetLeaf()){
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
      bool iKeepInputOrder,
      Double_t Rdistadapt,
      Double_t AdaptiveMedianFac
    )
    {
        iresetorder=true;
        ikeepinputorder = iKeepInputOrder;
        ibuildinparallel = false;
#ifdef USEOPENMP
        ibuildinparallel = iBuildInParallel;
        int inested = omp_get_max_active_levels();
        int nthreads;
        #pragma omp parallel
        #pragma omp single
        {
            nthreads = omp_get_num_threads();
        }
        if (nthreads == 1) ibuildinparallel = false;
        if (inested == 0 && ibuildinparallel) omp_set_max_active_levels(nthreads/2);
#endif
        numparts = nparts;
        numleafnodes=numnodes=0;
        bucket = p;
        b = bucket_size;
        bmin = std::max(static_cast<Int_t>(1),b/4);
        maxadaptivemedianregionsize = 8*b;
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (Rdistadapt > 0) rdist2adapt = Rdistadapt*Rdistadapt;
        else rdist2adapt = -1;
        // store rdist2apapt with and extra >1 factor which is useful for determining 
        // when to use non-median based splitting index 
        rdist2daptwithfac = 4.0*rdist2adapt;
        adaptivemedianfac = AdaptiveMedianFac;
        if (Period!=NULL)
        {
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=Period[k];
            PeriodReflectionConstruction();
        }
        else period=NULL;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            //if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            KDTreeOMPThreadPool otp = OMPInitThreadPool();
            root=BuildNodes(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
            //if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_max_active_levels(inested);
#endif
    }

    KDTree::KDTree(System &s, Int_t bucket_size,
      int ttype, int smfunctype, int smres, int criterion, int aniso, int scale,
      Double_t **m,
      bool iBuildInParallel,
      bool iKeepInputOrder,
      double Rdistadapt,
      Double_t AdaptiveMedianFac
    )
    {

        iresetorder=true;
        ikeepinputorder = iKeepInputOrder;
        ibuildinparallel = false;
#ifdef USEOPENMP
        ibuildinparallel = iBuildInParallel;
        int inested = omp_get_max_active_levels();
        int nthreads;
        #pragma omp parallel
        #pragma omp single
        {
            nthreads = omp_get_num_threads();
        }
        if (nthreads == 1) ibuildinparallel = false;
        if (inested == 0 && ibuildinparallel) omp_set_max_active_levels(nthreads/2);
#endif
        numparts = s.GetNumParts();
        numleafnodes=numnodes=0;
        bucket = s.Parts();
        b = bucket_size;
        bmin = std::max(static_cast<Int_t>(1),b/4);
        maxadaptivemedianregionsize = 8*b;
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (Rdistadapt > 0) rdist2adapt = Rdistadapt*Rdistadapt;
        else rdist2adapt = -1;
        // store rdist2apapt with and extra >1 factor which is useful for determining 
        // when to use non-median based splitting index 
        rdist2daptwithfac = 4.0*rdist2adapt;
        adaptivemedianfac = AdaptiveMedianFac;
        if (s.GetPeriod()[0]>0&&s.GetPeriod()[1]>0&&s.GetPeriod()[2]>0){
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=s.GetPeriod()[k];
            PeriodReflectionConstruction();
        }
        else period=NULL;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            //if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            KDTreeOMPThreadPool otp = OMPInitThreadPool();
            root=BuildNodes(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            //if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_max_active_levels(inested);
#endif
    }

    KDTree::~KDTree()
    {
	    if (root!=NULL) {
            delete root;
            delete[] Kernel;
            delete[] derKernel;
            if (period!=NULL) delete[] period;
            if (iresetorder) std::sort(bucket, bucket + numparts, IDCompareVec);
            if (scalespace) {
            for (Int_t i=0;i<numparts;i++)
                for (int j=0;j<3;j++) {
                    bucket[i].SetPosition(j,bucket[i].GetPosition(j)*xvar[j]);
                    bucket[i].SetVelocity(j,bucket[i].GetVelocity(j)*xvar[j+3]);
                }
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
