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

    // -- Inline functions that get called often when building the tree.

    /// \name Find most spread dimension
    //@{
    inline Double_t KDTree::SpreadestPos(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t minval = bucket[start].GetPosition(j);
        Double_t maxval = minval;
        Int_t i;
        unsigned int nthreads;;
#ifdef USEOPENMP
        nthreads = min((unsigned int)(floor((end-start)/float(KDTREEOMPCRITPARALLELSIZE))), otp.nactivethreads);
        if (nthreads <1) nthreads=1;
#pragma omp parallel for \
default(shared) private(i) schedule(static) \
reduction(min:minval) reduction(max:maxval) num_threads(nthreads) if (nthreads>1)
#endif
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < minval) minval = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > maxval) maxval = bucket[i].GetPosition(j);
        }
        bnd[0]=minval;bnd[1]=maxval;
        return maxval - minval;
    }
    inline Double_t KDTree::SpreadestVel(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t minval = bucket[start].GetVelocity(j);
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
            if (bucket[i].GetVelocity(j) < minval) minval = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > maxval) maxval = bucket[i].GetVelocity(j);
        }
        bnd[0]=minval;bnd[1]=maxval;
        return maxval - minval;
    }
    inline Double_t KDTree::SpreadestPhs(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t minval = bucket[start].GetPhase(j);
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
            if (bucket[i].GetPhase(j) < minval) minval = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > maxval) maxval = bucket[i].GetPhase(j);
        }
        bnd[0]=minval;bnd[1]=maxval;
        return maxval - minval;
    }
    //@}
    /// \name Find the boundary of the data and return mean
    //@{
    inline Double_t KDTree::BoundaryandMeanPos(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t mean=bucket[start].GetPosition(j);
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
            if (bucket[i].GetPosition(j) < minval) minval = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > maxval) maxval = bucket[i].GetPosition(j);
            mean+=bucket[i].GetPosition(j);
        }
        bnd[0]=minval;bnd[1]=maxval;
        mean/=(Double_t)(end-start);
        return mean;
    }
    inline Double_t KDTree::BoundaryandMeanVel(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t mean=bucket[start].GetVelocity(j);
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
            if (bucket[i].GetVelocity(j) < minval) minval = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > maxval) maxval = bucket[i].GetVelocity(j);
            mean+=bucket[i].GetVelocity(j);
        }
        bnd[0]=minval;bnd[1]=maxval;
        mean/=(Double_t)(end-start);
        return mean;
    }
    inline Double_t KDTree::BoundaryandMeanPhs(int j, Int_t start, Int_t end, Double_t *bnd,
        KDTreeOMPThreadPool &otp)
    {
        Double_t mean=bucket[start].GetPhase(j);
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
            if (bucket[i].GetPhase(j) < minval) minval = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > maxval) maxval = bucket[i].GetPhase(j);
            mean+=bucket[i].GetPhase(j);
        }
        bnd[0]=minval;bnd[1]=maxval;
        mean/=(Double_t)(end-start);
        return mean;
    }
    //@}
    /// \name Find the dispersion in a dimension (biased variance using 1/N as opposed to 1/(N-1) so that if N=2, doesn't crash)
    //@{
    inline Double_t KDTree::DispersionPos(int j, Int_t start, Int_t end, Double_t mean,
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
            disp+=(bucket[i].GetPosition(j)-mean)*(bucket[i].GetPosition(j)-mean);
        disp/=(Double_t)(end-start);
        return disp;
    }
    inline Double_t KDTree::DispersionVel(int j, Int_t start, Int_t end, Double_t mean,
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
        for (i = start ; i < end; i++)
            disp+=(bucket[i].GetVelocity(j)-mean)*(bucket[i].GetVelocity(j)-mean);
        disp/=(Double_t)(end-start);
        return disp;
    }
    inline Double_t KDTree::DispersionPhs(int j, Int_t start, Int_t end, Double_t mean,
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
            disp+=(bucket[i].GetPhase(j)-mean)*(bucket[i].GetPhase(j)-mean);
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
    inline Double_t KDTree::EntropyPos(int j, Int_t start, Int_t end,
        Double_t low, Double_t up, Double_t nbins, Double_t *nientropy,
        KDTreeOMPThreadPool &otp)
    {
        Int_t ibin, i;
        Double_t mtot=0.,entropy=0.;
        Double_t idx=static_cast<Double_t>(nbins)/(up-low);
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=static_cast<Int_t>(floor((bucket[i].GetPosition(j)-low)*idx));
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
    inline Double_t KDTree::EntropyVel(int j, Int_t start, Int_t end,
        Double_t low, Double_t up, Double_t nbins, Double_t *nientropy,
        KDTreeOMPThreadPool &otp)
    {
        Int_t ibin,i;
        Double_t mtot=0.,entropy=0.;
        Double_t idx=static_cast<Double_t>(nbins)/(up-low);
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=static_cast<Int_t>(floor((bucket[i].GetVelocity(j)-low)*idx));
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
    inline Double_t KDTree::EntropyPhs(int j, Int_t start, Int_t end,
        Double_t low, Double_t up, Double_t nbins, Double_t *nientropy,
        KDTreeOMPThreadPool &otp)
    {
        Int_t ibin,i;
        Double_t mtot=0.,entropy=0.;
        Double_t idx=static_cast<Double_t>(nbins)/(up-low);
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=static_cast<Int_t>(floor((bucket[i].GetPhase(j)-low)*idx));
            nientropy[ibin]+=bucket[i].GetMass();
        }
        mtot=1.0/mtot;
        for (Int_t i=0;i<nbins;i++) {
            if (nientropy[i]>0) {
                Double_t temp=nientropy[i]*mtot;
                entropy-=temp*log10(temp);
            }
        }
        return entropy/log10(nbins);
    }
    //@}

    /// \name Determine the median coordinates in some space
    //@{
    Double_t KDTree::MedianPos(int d, Int_t &k, Int_t start, Int_t end,
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
                x = bucket[k].GetPosition(d);
                swap(bucket[right],bucket[k]);
                pval = &bucket[k];
                i = left-1;
                j = right;
                while (1) {
                    while (i < j) if (bucket[++i].GetPosition(d) >= x) break;
                    while (i < j) if (bucket[--j].GetPosition(d) <= x) break;
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
            return bucket[k].GetPosition(d);
        }
        //requires that particle order is already balanced. Use with caution
        else
        {
            return bucket[k].GetPosition(d);
            //printf("Note yet implemented\n");
            //exit(9);
        }
    }
    Double_t KDTree::MedianVel(int d, Int_t &k, Int_t start, Int_t end,
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Int_t left = start;
        Int_t right = end - 1;
        Int_t i, j;
        Double_t x;
        Particle w;
        Particle *pval = NULL;

        if (balanced){
            while (left < right)
            {
                x = bucket[k].GetVelocity(d);
                swap(bucket[right],bucket[k]);
                pval = &bucket[k];
                i = left-1;
                j = right;
                while (1) {
                    while (i < j) if (bucket[++i].GetVelocity(d) >= x) break;
                    while (i < j) if (bucket[--j].GetVelocity(d) <= x) break;
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
            return bucket[k].GetVelocity(d);
        }
        //requires that particle order is already balanced. Use with caution
        else
        {
            return bucket[k].GetVelocity(d);
            //printf("Note yet implemented\n");
            //exit(9);
        }
    }
    Double_t KDTree::MedianPhs(int d, Int_t &k, Int_t start, Int_t end,
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Int_t left = start;
        Int_t right = end - 1;
        Int_t i, j;
        Double_t x;
        Particle w;
        Particle *pval = NULL;

        if (balanced){
            while (left < right)
            {
                x = bucket[k].GetPhase(d);
                swap(bucket[right],bucket[k]);
                pval = &bucket[k];
                i = left-1;
                j = right;
                while (1) {
                    while (i < j) {if (bucket[++i].GetPhase(d) >= x) break;}
                    while (i < j) {if (bucket[--j].GetPhase(d) <= x) break;}
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
            return bucket[k].GetPhase(d);
        }
        //requires that particle order is already balanced. Use with caution
        else
        {
            return bucket[k].GetPhase(d);
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
            if(splittingcriterion==1) {
                splitvalue[j] = (this->*spreadfunc)(j, start, end, bnd[j], otp)+1e-32;//addition incase lattice and no spread
                Double_t low, up;
                low=bnd[j][0]-2.0*(splitvalue[j])/(Double_t)(end-start);
                up=bnd[j][1]+2.0*(splitvalue[j])/(Double_t)(end-start);
                splitvalue[j] = (this->*entropyfunc)(j, start, end, low, up, nbins, entropybins.data(), otp);
            }
            else if (splittingcriterion==2) {
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
    ///adjust the sorted particle array so that the split is not at the median
    ///but approximative and splits at the particle with the largest distance
    ///between particles. This search is limited to a buffer region around
    ///the median index
    Double_t KDTree::AdjustMedianToMaximalDistancePos(int d,
        Int_t &splitindex, Int_t trueleft, Int_t trueright,
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Double_t splitvalue = MedianPos(d, splitindex, trueleft, trueright, otp, balanced);
        UInt_tree_t truesize = trueright - trueleft;
        UInt_tree_t bufferwidth = truesize * adaptivemedianfac;
        if (bufferwidth<minadaptivemedianregionsize) return splitvalue;
        UInt_tree_t left = splitindex - bufferwidth/2;
        UInt_tree_t right = splitindex + bufferwidth/2;
        UInt_tree_t size = right - left;
        vector<KDTreeForSorting> x(size);
        for (auto i=0;i<size;i++) {
            x[i].val = bucket[left+i].GetPosition(d);
            x[i].orgindex = left+i;
        }
        UInt_tree_t n=0;
        std::sort(x.begin(), x.end() , [](const KDTreeForSorting &a, const KDTreeForSorting &b) {
            return a.val < b.val;
        });
        UInt_tree_t newsplitindex = left;
        Double_t newsplitvalue;
        auto dist = std::abs(x[1].val - x[0].val);
        auto maxdist = dist;
        UInt_tree_t maxi = 0;
        for (UInt_tree_t i=1; i<size-1; i++)
        {
            dist = std::abs(x[i+1].val - x[i].val);
            if (dist > maxdist)
            {
                maxdist = dist;
                newsplitindex = i+x[i].orgindex;
                newsplitvalue = x[i].val;
                maxi = i;
            }
        }
        splitindex = newsplitindex;
        splitvalue = newsplitvalue;
        splitvalue = MedianPos(d, splitindex, trueleft, trueright, otp, balanced);
        return splitvalue;

        /*
        //sort buffer region
        std::sort(&bucket[left], &bucket[left] + size, [d](const Particle &a, const Particle &b) {
            return a.GetPosition(d) < b.GetPosition(d);
        });
        for (auto i=0;i<size;i++) x[i] = bucket[left+i].GetPosition(d);
        UInt_tree_t newsplitindex = left;
        Double_t newsplitvalue;
        auto dist = std::abs(x[1] - x[0]);
        auto maxdist = dist;
        UInt_tree_t maxi = 0;
        for (UInt_tree_t i=1; i<size-1; i++)
        {
            dist = std::abs(x[i+1] - x[i]);
            if (dist > maxdist)
            {
                maxdist = dist;
                newsplitindex = i+left;
                newsplitvalue = x[i];
                maxi = i;
            }
        }
        splitindex = newsplitindex;
        splitvalue = newsplitvalue;
        return splitvalue;
        */
    }
    Double_t KDTree::AdjustMedianToMaximalDistanceVel(int d,
        Int_t &splitindex, Int_t trueleft, Int_t trueright,
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Double_t splitvalue = MedianVel(d, splitindex, trueleft, trueright, otp, balanced);
        UInt_tree_t truesize = trueright - trueleft;
        UInt_tree_t bufferwidth = truesize * adaptivemedianfac;
        if (bufferwidth<minadaptivemedianregionsize) return splitvalue;
        UInt_tree_t left = splitindex - bufferwidth/2;
        UInt_tree_t right = splitindex + bufferwidth/2;
        UInt_tree_t size = right - left;
        vector<KDTreeForSorting> x(size);
        for (auto i=0;i<size;i++) {
            x[i].val = bucket[left+i].GetVelocity(d);
            x[i].orgindex = left+i;
        }
        UInt_tree_t n=0;
        std::sort(x.begin(), x.end() , [](const KDTreeForSorting &a, const KDTreeForSorting &b) {
            return a.val < b.val;
        });
        UInt_tree_t newsplitindex = left;
        Double_t newsplitvalue;
        auto dist = std::abs(x[1].val - x[0].val);
        auto maxdist = dist;
        UInt_tree_t maxi = 0;
        for (UInt_tree_t i=1; i<size-1; i++)
        {
            dist = std::abs(x[i+1].val - x[i].val);
            if (dist > maxdist)
            {
                maxdist = dist;
                newsplitindex = i+x[i].orgindex;
                newsplitvalue = x[i].val;
                maxi = i;
            }
        }
        splitindex = newsplitindex;
        splitvalue = newsplitvalue;
        splitvalue = MedianVel(d, splitindex, trueleft, trueright, otp, balanced);
        return splitvalue;
    }
    Double_t KDTree::AdjustMedianToMaximalDistancePhs(int d,
        Int_t &splitindex, Int_t trueleft, Int_t trueright,
        KDTreeOMPThreadPool &otp, bool balanced)
    {
        Double_t splitvalue = MedianPhs(d, splitindex, trueleft, trueright, otp, balanced);
        UInt_tree_t truesize = trueright - trueleft;
        UInt_tree_t bufferwidth = truesize * adaptivemedianfac;
        if (bufferwidth<minadaptivemedianregionsize) return splitvalue;
        UInt_tree_t left = splitindex - bufferwidth/2;
        UInt_tree_t right = splitindex + bufferwidth/2;
        UInt_tree_t size = right - left;
        vector<KDTreeForSorting> x(size);
        for (auto i=0;i<size;i++) {
            x[i].val = bucket[left+i].GetPhase(d);
            x[i].orgindex = left+i;
        }
        UInt_tree_t n=0;
        std::sort(x.begin(), x.end() , [](const KDTreeForSorting &a, const KDTreeForSorting &b) {
            return a.val < b.val;
        });
        UInt_tree_t newsplitindex = left;
        Double_t newsplitvalue;
        auto dist = std::abs(x[1].val - x[0].val);
        auto maxdist = dist;
        UInt_tree_t maxi = 0;
        for (UInt_tree_t i=1; i<size-1; i++)
        {
            dist = std::abs(x[i+1].val - x[i].val);
            if (dist > maxdist)
            {
                maxdist = dist;
                newsplitindex = i+x[i].orgindex;
                newsplitvalue = x[i].val;
                maxi = i;
            }
        }
        splitindex = newsplitindex;
        splitvalue = newsplitvalue;
        splitvalue = MedianPhs(d, splitindex, trueleft, trueright, otp, balanced);
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
                for(auto j=0;j<ND;j++) localcenter[j] += bucket[i].GetPhase(j);
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
                for(auto j=0;j<ND;j++) center[j] += bucket[i].GetPhase(j);
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
                for(auto j=0;j<ND;j++) pos[j] = bucket[i].GetPhase(j);
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
                for(auto j=0;j<ND;j++) pos[j] = bucket[i].GetPhase(j);
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
        // get center
        Double_t maxinterdist = 0.0;
        UInt_tree_t size = (localend - localstart);
        vector<KDTreeForSorting> x(size);
        for (auto i=0; i<size; i++) {
            x[i].val = bucket[i + localstart].GetPhase(splitdim);
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
                auto diff = pow(x[i+1].val - x[i].val,2.0);
                localmax = std::max(localmax, diff);
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
                auto diff = pow(x[i+1].val - x[i].val,2.0);
                maxinterdist = std::max(maxinterdist, diff);
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
        Double_t localfarthest;
        // if constructing adaptive tree where leaf nodes must be smaller than some size
        // calculate the farthest distance to the centre of the node
        if (rdist2adapt > 0) {
            center = DetermineCentreAndSmallestSphere(start, end, localfarthest, otp);
            // if checking that leaf nodes have interparticle spacings smaller than some value
            // then get interparticle spacing
            if (igetmaxinterparticlespacing) {
                //first get split dim
                splitdim = DetermineSplitDim(start, end, bnd, otp);
                //then get maximum interparticle spacing in split dimension
                auto maxinterdist = DetermineMaxInterParticleSpacing(start, end, splitdim, otp);
                isleafflag = ((size <= b && maxinterdist < rdist2adapt) || (size <= bmin));
            }
            // otherwise splitting criterion based on just farthest
    	    else {
                isleafflag = ((size <= b && localfarthest < rdist2adapt) || (size <= bmin));
            }
        }
        else {
            isleafflag = (size <= b);
        }
        if (isleafflag)
        {
            if (ibuildinparallel == false) numleafnodes++;
            for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j], otp);
            Node * leaf = new LeafNode(id, start, end,  bnd, ND);
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
            Double_t splitvalue = (this->*medianfunc)(splitdim, splitindex, start, end, otp, irearrangeandbalance);
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
                    left = BuildNodes(start, splitindex+1, newotp[0]);
                    if (rdist2adapt>0) DetermineCentreAndSmallestSphere(start, splitindex+1, left, newotp[0]);
                    #pragma omp task
                    right = BuildNodes(splitindex+1, end, newotp[1]);
                    if (rdist2adapt>0) DetermineCentreAndSmallestSphere(splitindex+1, end, right, newotp[1]);
                    #pragma omp taskwait
                }
#endif
                return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    left, right);
            }
            else {
                if (rdist2adapt > 0) {
                    Node *left, *right;
                    left = BuildNodes(start, splitindex+1, otp);
                    DetermineCentreAndSmallestSphere(start, splitindex+1, left, otp);
                    right = BuildNodes(splitindex+1, end, otp);
                    DetermineCentreAndSmallestSphere(splitindex+1, end, right, otp);
                    return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    left, right);
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
        if (treetype==TPHYS)
        {
            bmfunc=&NBody::KDTree::BoundaryandMeanPos;
            dispfunc=&NBody::KDTree::DispersionPos;
            spreadfunc=&NBody::KDTree::SpreadestPos;
            entropyfunc=&NBody::KDTree::EntropyPos;
            if (adaptivemedianfac > 0) medianfunc=&NBody::KDTree::AdjustMedianToMaximalDistancePos;
            else medianfunc=&NBody::KDTree::MedianPos;
        }
        else if (treetype==TVEL) {
            bmfunc=&NBody::KDTree::BoundaryandMeanVel;
            dispfunc=&NBody::KDTree::DispersionVel;
            spreadfunc=&NBody::KDTree::SpreadestVel;
            entropyfunc=&NBody::KDTree::EntropyVel;
            if (adaptivemedianfac > 0) medianfunc=&NBody::KDTree::AdjustMedianToMaximalDistanceVel;
            else medianfunc=&NBody::KDTree::MedianVel;
        }
        else if (treetype==TPHS) {
            bmfunc=&NBody::KDTree::BoundaryandMeanPhs;
            dispfunc=&NBody::KDTree::DispersionPhs;
            spreadfunc=&NBody::KDTree::SpreadestPhs;
            entropyfunc=&NBody::KDTree::EntropyPhs;
            if (adaptivemedianfac > 0) medianfunc=&NBody::KDTree::AdjustMedianToMaximalDistancePhs;
            else medianfunc=&NBody::KDTree::MedianPhs;
        }
        else if (treetype==TPROJ) {
            bmfunc=&NBody::KDTree::BoundaryandMeanPos;
            dispfunc=&NBody::KDTree::DispersionPos;
            spreadfunc=&NBody::KDTree::SpreadestPos;
            entropyfunc=&NBody::KDTree::EntropyPos;
            if (adaptivemedianfac > 0) medianfunc=&NBody::KDTree::AdjustMedianToMaximalDistancePos;
            else medianfunc=&NBody::KDTree::MedianPos;
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
        cout<<"At node "<<" "<<id<<" "<<node->GetLeaf()<<" "<<start<<" "<<end<<" "<<size<<" ";
        if (!node->GetLeaf()) cout<<((SplitNode*)node)->GetCutDim()<<" "<<((SplitNode*)node)->GetCutValue()<<" ";
        for (auto j=0;j<ND;j++)  cout<<"("<<node->GetBoundary(j,0)<<", "<<node->GetBoundary(j,1)<<")";
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
      double Rdistadapt,
      Double_t AdaptiveMedianFac,
      bool iGetMaxInterParticleSpacing
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
        numparts = nparts;
        numleafnodes=numnodes=0;
        bucket = p;
        b = bucket_size;
        bmin = std::max(static_cast<Int_t>(1),b/4);
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (Rdistadapt > 0) rdist2adapt = Rdistadapt*Rdistadapt;
        else rdist2adapt = -1;
        adaptivemedianfac = AdaptiveMedianFac;
        igetmaxinterparticlespacing = (Rdistadapt > 0 && iGetMaxInterParticleSpacing);
        if (Period!=NULL)
        {
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=Period[k];
        }
        else period=NULL;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            KDTreeOMPThreadPool otp = OMPInitThreadPool();
            root=BuildNodes(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
            if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_nested(inested);
#endif
    }

    KDTree::KDTree(System &s, Int_t bucket_size,
      int ttype, int smfunctype, int smres, int criterion, int aniso, int scale,
      Double_t **m,
      bool iBuildInParallel,
      bool iKeepInputOrder,
      double Rdistadapt,
      Double_t AdaptiveMedianFac,
      bool iGetMaxInterParticleSpacing
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
        bmin = std::max(static_cast<Int_t>(1),b/4);
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (Rdistadapt > 0) rdist2adapt = Rdistadapt*Rdistadapt;
        else rdist2adapt = -1;
        adaptivemedianfac = AdaptiveMedianFac;
        igetmaxinterparticlespacing = (Rdistadapt > 0 && iGetMaxInterParticleSpacing);
        if (s.GetPeriod()[0]>0&&s.GetPeriod()[1]>0&&s.GetPeriod()[2]>0){
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=s.GetPeriod()[k];
        }
        else period=NULL;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            KDTreeOMPThreadPool otp = OMPInitThreadPool();
            root=BuildNodes(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_nested(inested);
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
