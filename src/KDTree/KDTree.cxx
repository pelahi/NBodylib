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
        Double_t dx=(up-low)/nbins;
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=(Int_t)((bucket[i].GetPosition(j)-low)/dx);
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
        Double_t dx=(up-low)/nbins;
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=(Int_t)((bucket[i].GetVelocity(j)-low)/dx);
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
        Double_t dx=(up-low)/nbins;
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=(Int_t)((bucket[i].GetPhase(j)-low)/dx);
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
    inline Double_t KDTree::MedianPos(int d, Int_t k, Int_t start, Int_t end,
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
    inline Double_t KDTree::MedianVel(int d, Int_t k, Int_t start, Int_t end,
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
    inline Double_t KDTree::MedianPhs(int d, Int_t k, Int_t start, Int_t end,
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
    inline int KDTree::DetermineSplitDim(Int_t start, Int_t end, Double_t bnd[6][2],
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
            entropybins.resize(nbins);
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

    inline void KDTree::js_qsort(int js_start, int js_end, int js_dim)
    {
	    int js_ind = (js_start + js_end)/2;
	    int js_i = js_start;
	    int js_j = js_end;
	    Double_t js_xx = bucket[js_ind].GetPhase(js_dim);
	    while(1)
	    {
		    while(bucket[js_i].GetPhase(js_dim) < js_xx) js_i++;
		    while(bucket[js_j].GetPhase(js_dim) > js_xx) js_j--;

		    if(js_i >= js_j) break;
		    swap(bucket[js_i],bucket[js_j]);
		    js_i++;
		    js_j--;
	    }
	    if(js_start < js_i - 1) js_qsort(js_start, js_i-1, js_dim);
	    if(js_j +1 < js_end) js_qsort(js_j + 1, js_end, js_dim);
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
            for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j], otp);

	    LeafNode *js_LN;
	    js_LN = new LeafNode(id ,start, end,  bnd, ND);
	    js_LN->SetLeaf(1);

	    return js_LN;
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

    ///--JS--
    /// Build OMP Domains by using Adaptive KDTree. OMP domains are decomposed by considering
    /// interparticle distances to avoid an overhead when linking across OMP domains
    Node *KDTree::BuildNodes_OMP(Int_t start, Int_t end, KDTreeOMPThreadPool &otp)
    {
	    Double_t bnd[6][2];
	    Int_t size = end - start, k;
	    Int_tree_t id = 0;
	    int splitdim;
	    Double_t splitvalue;
	    int js_ind0 = start;
	    int js_ind1 = end-1;
	    int js_ompskip=-1;

        
	    //if not building in parallel can set ids here and update number of nodes
	    //otherwise, must set after construction
	    if (ibuildinparallel == false) {
		    id = numnodes;
		    numnodes++;
	    }

	    //See whether this OMP region is further decomposed or not.
	    //If the maximum interparticle distance is quite larger than the linking length,
	    // decompose the domain further.
	    b=200000000;
	    if (size <= b){
		    js_ompskip=1;
		    if(size > 1000000){ // Minimum size of OMP domain (arbitrary chosen)
			    double js_dx=0., js_dx2;
			    int js_nn = (end - start) / 8; // Buffer for splitting
			    for(int js_dim=0; js_dim<ND; js_dim++){
				    js_qsort(js_ind0, js_ind1, js_dim);
				    for(int js_ind=start + js_nn; js_ind<end - js_nn; js_ind++){
					    js_dx2 = abs(bucket[js_ind+1].GetPhase(js_dim) - bucket[js_ind].GetPhase(js_dim));
					    if(js_dx2 > js_dx) js_dx = js_dx2;
				    }

				    if(js_dx > 2.0*js_rdist){js_ompskip=-1; break;}
				    js_dx = 0.;
			    }
		    }
	    }

	    //Leaf Node Construction
	    if (js_ompskip>0){
		    if (ibuildinparallel == false) numleafnodes++;
		    for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j], otp);

		    LeafNode *js_LN;
		    js_LN = new LeafNode(id, start, end, bnd, ND);
		    js_LN->SetLeaf(1);

		    return js_LN;
	    }
	    else
	    {
		    bool irearrangeandbalance=true;
		    if (ikeepinputorder) irearrangeandbalance=false;

		    splitdim = DetermineSplitDim(start, end, bnd, otp);
		    js_qsort(js_ind0, js_ind1, splitdim);

		    double js_dx=0., js_dx2;
		    int js_nn = (end - start) / 8;

		    for(int js_ind=start + js_nn; js_ind<end - js_nn; js_ind++){
			    js_dx2 = abs(bucket[js_ind+1].GetPhase(splitdim) - bucket[js_ind].GetPhase(splitdim));
			    if(js_dx2 > js_dx){js_dx=js_dx2; k=js_ind; splitvalue=bucket[k].GetPhase(splitdim);}
		    }
	    }

	    //Now Split the node
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
			    left = BuildNodes_OMP(start, k+1, newotp[0]);
			    #pragma omp task
			    right = BuildNodes_OMP(k+1, end, newotp[1]);
			    #pragma omp taskwait
		    }


                    SplitNode *end_node = new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);
                    left->SetSibling(right);
                    right->SetSibling(left);
                    left->SetParent(end_node);
                    right->SetParent(end_node);
                    return end_node;

                    //return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    //    left, right);
#endif
	    }
	    else {
		    Node *left, *right;

		    left = BuildNodes_OMP(start, k+1, otp);
		    right = BuildNodes_OMP(k+1, end, otp);

                    SplitNode *end_node = new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);
                    left->SetSibling(right);
                    right->SetSibling(left);
                    left->SetParent(end_node);
                    right->SetParent(end_node);
                    return end_node;
                    //return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND,
                    //    left, right);
	    }	    
    }
    
    ///--JS--
    /// Build Normal Adaptive KDTree
    Node *KDTree::BuildNodes_ADT(Int_t start, Int_t end, KDTreeOMPThreadPool &otp)
    {
	    Double_t bnd[6][2];
	    Int_t size = end - start, k;
	    Int_tree_t id = 0;
	    int splitdim;
	    Double_t splitvalue;
	    int js_ind0 = start;
	    int js_ind1 = end-1;

	    //if not building in parallel can set ids here and update number of nodes
	    //otherwise, must set after construction
	    if (ibuildinparallel == false) {
		    id = numnodes;
		    numnodes++;
	    }

	    //Leaf Node Construction
	    if (size <= b){
		    if (ibuildinparallel == false) numleafnodes++;
		    for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j], otp);

		    LeafNode *js_LN;
		    js_LN = new LeafNode(id, start, end, bnd, ND);
		    js_LN->SetLeaf(1);

		    return js_LN;
	    }
	    else
	    {
		    bool irearrangeandbalance=true;
		    if (ikeepinputorder) irearrangeandbalance=false;

		    ///
		    //double js_dx=0., js_dx2;
		    //int js_skip = 0, js_dim2;
		    //int js_nn = (end - start) / 8;
		    //for(int js_dim=0; js_dim<ND; js_dim++){
		    //    js_dim2 = js_dim;
		    //	js_qsort(js_ind0, js_ind1, js_dim2);
		    //    for(int js_ind=start + js_nn; js_ind<end - js_nn; js_ind++){
		    //    	js_dx2 = abs(bucket[js_ind+1].GetPhase(js_dim) - bucket[js_ind].GetPhase(js_dim));
		    //    	if(js_dx2 > js_dx){
		    //    		splitdim = js_dim;
		    //    		k = js_ind;
		    //    		splitvalue=bucket[k].GetPhase(splitdim);
		    //    		js_dx = js_dx2;
		    //    	}
		    //    }
		    //}
		    //js_qsort(js_ind0, js_ind1, splitdim);
		    ///

		    splitdim = DetermineSplitDim(start, end, bnd, otp);
		    js_qsort(js_ind0, js_ind1, splitdim);

		    double js_dx=0., js_dx2;
		    int js_nn = (end - start) / 8;
		    if (js_nn == 0) js_nn = 1;

		    for(int js_ind=start + js_nn; js_ind<end - js_nn; js_ind++){
		            js_dx2 = abs(bucket[js_ind+1].GetPhase(splitdim) - bucket[js_ind].GetPhase(splitdim));
		            if(js_dx2 > js_dx){js_dx=js_dx2; k=js_ind; splitvalue=bucket[k].GetPhase(splitdim);}
		    }
		    

		    //k = start + (size - 1) / 2;
		    //splitvalue = (this->*medianfunc)(splitdim, k, start, end, otp, irearrangeandbalance);
	    }

	    //Now Split the node
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
			    left = BuildNodes_ADT(start, k+1, newotp[0]);
			    #pragma omp task
			    right = BuildNodes_ADT(k+1, end, newotp[1]);
			    #pragma omp taskwait
		    }

		    //Now save the largest distance from the center of this node
		    //Here, the center is just defined as the mean coordinates of the included particles
		    //TO DO LIST
		    //	Define center by using the 'finding the smallest sphere' algorithm
		    //	(https://en.wikipedia.org/wiki/Smallest-circle_problem)

		    //Left
		    Double_t js_center[ND], js_centertmp, js_pos[ND];
		    Double_t js_dd=-1, js_dd2;

		    //Left - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(int js_i=start; js_i<k+1; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (k+1 - start));
		    }
		    for(int js_i=0; js_i<ND; js_i++) left->SetCenter(js_center[js_i],js_i);

		    //Left - largest distance
		    for(int js_i=start; js_i<k+1; js_i++){
			    for(int js_j=0; js_j<ND; js_j++) js_pos[js_j] = bucket[js_i].GetPhase(js_j);
			    js_dd2 = DistanceSqd(js_pos, js_center, ND);
			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    left->SetFarthest(js_dd);

		    //Right
		    js_dd = -1.;
		    //Right - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(int js_i=k+1; js_i<end; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (end - (k+1)));
		    }
		    for(int js_i=0; js_i<ND; js_i++) right->SetCenter(js_center[js_i],js_i);

		    //Right - largest distance
		    for(int js_i=k+1; js_i<end; js_i++){
			    for(int js_j=0; js_j<ND; js_j++) js_pos[js_j] = bucket[js_i].GetPhase(js_j);
			    js_dd2 = DistanceSqd(js_pos, js_center, ND);
			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    right->SetFarthest(js_dd);

		    SplitNode *end_node = new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);
		    left->SetSibling(right);
		    right->SetSibling(left);
		    left->SetParent(end_node);
		    right->SetParent(end_node);
		    return end_node;
		    //return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);

#endif
	    }
	    else {
		    Node *left, *right;

		    left = BuildNodes_ADT(start, k+1, otp);
		    right = BuildNodes_ADT(k+1, end, otp);

		    //Now save the largest distance from the center of this node
		    //Here, the center is just defined as the mean coordinates of the included particles
		    //TO DO LIST
		    //	Define center by using the 'finding the smallest sphere' algorithm
		    //	(https://en.wikipedia.org/wiki/Smallest-circle_problem)

		    //Left
		    Double_t js_center[ND], js_centertmp, js_pos[ND];
		    Double_t js_dd=-1, js_dd2;

		    //Left - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(int js_i=start; js_i<k+1; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (k+1 - start));
		    }
		    for(int js_i=0; js_i<ND; js_i++) left->SetCenter(js_center[js_i],js_i);

		    //Left - largest distance
		    for(int js_i=start; js_i<k+1; js_i++){
			    for(int js_j=0; js_j<ND; js_j++) js_pos[js_j] = bucket[js_i].GetPhase(js_j);
			    js_dd2 = DistanceSqd(js_pos, js_center, ND);
			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    left->SetFarthest(js_dd);

		    //Right
		    js_dd = -1.;
		    //Right - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(int js_i=k+1; js_i<end; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (end - (k+1)));
		    }
		    for(int js_i=0; js_i<ND; js_i++) right->SetCenter(js_center[js_i],js_i);

		    //Right - largest distance
		    for(int js_i=k+1; js_i<end; js_i++){
			    for(int js_j=0; js_j<ND; js_j++) js_pos[js_j] = bucket[js_i].GetPhase(js_j);
			    js_dd2 = DistanceSqd(js_pos, js_center, ND);
			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    right->SetFarthest(js_dd);

		    SplitNode *end_node = new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);
		    left->SetSibling(right);
		    right->SetSibling(left);
		    left->SetParent(end_node);
		    right->SetParent(end_node);
		    return end_node;
		    //return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);
	    }	    
    }
    
    ///--JS--
    /// Build Normal Adaptive KDTree for FOFSearchCriterion
    Node *KDTree::BuildNodes_CRIT(Int_t start, Int_t end, KDTreeOMPThreadPool &otp, Double_t *param)
    {
	    Double_t bnd[6][2];
	    Int_t size = end - start, k;
	    Int_tree_t id = 0;
	    int splitdim;
	    Double_t splitvalue;
	    int js_ind0 = start;
	    int js_ind1 = end-1;

        
	    //if not building in parallel can set ids here and update number of nodes
	    //otherwise, must set after construction
	    if (ibuildinparallel == false) {
		    id = numnodes;
		    numnodes++;
	    }

	    //Leaf Node Construction
	    if (size <= b){
		    if (ibuildinparallel == false) numleafnodes++;
		    for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j], otp);

		    LeafNode *js_LN;
		    js_LN = new LeafNode(id, start, end, bnd, ND);
		    js_LN->SetLeaf(1);

		    return js_LN;
	    }
	    else
	    {
		    bool irearrangeandbalance=true;
		    if (ikeepinputorder) irearrangeandbalance=false;

		    //--JS--
		    //Should scale phase-space coordinates first otherwise the split dimension
		    // incorrectly chosen, which results in very poor performance in the core 
		    // search.
		    Double_t js_xc = 1.0/sqrt(param[1]);
		    Double_t js_vc = 1.0/sqrt(param[2]);

		    for(Int_t js_i=start; js_i<end; js_i++) bucket[js_i].ScalePhase(js_xc, js_vc);
		    splitdim = DetermineSplitDim(start, end, bnd, otp);

		    js_xc = 1.0/js_xc;
		    js_vc = 1.0/js_vc;
		    for(Int_t js_i=start; js_i<end; js_i++) bucket[js_i].ScalePhase(js_xc, js_vc);
		    js_qsort(js_ind0, js_ind1, splitdim);

		    double js_dx=0., js_dx2;
		    int js_nn = (end - start) / 8;

		    for(int js_ind=start + js_nn; js_ind<end - js_nn; js_ind++){
			    js_dx2 = abs(bucket[js_ind+1].GetPhase(splitdim) - bucket[js_ind].GetPhase(splitdim));
			    if(js_dx2 > js_dx){js_dx=js_dx2; k=js_ind; splitvalue=bucket[k].GetPhase(splitdim);}
		    }

		    //Adaptive Tree based on background particles only
		    //k = start + (size - 1) / 2;
		    //splitvalue = (this->*medianfunc)(splitdim, k, start, end, otp, irearrangeandbalance);

		    //double js2_dx0, js2_dx1, js2_dx2=-1.;
		    //int js2_nn=0;
		    //for(int js_ind=start + js_nn; js_ind<end - js_nn; js_ind++){
		    //        if(bucket[js_ind].GetPotential()<param[9]){
		    //    	    js2_dx1 = bucket[js_ind].GetPhase(splitdim);
		    //    	    if(js2_nn==0){
		    //    		    js2_nn++;
		    //    		    js2_dx0=js2_dx1;
		    //    	    	    k= js_ind;
		    //    		    splitvalue=bucket[k].GetPhase(splitdim);
		    //    	    }
		    //    	    else
		    //    	    {
		    //    		    if(js2_dx1 - js2_dx0 > js2_dx2){
		    //    			    k = js_ind;
		    //    			    splitvalue=bucket[k].GetPhase(splitdim);
		    //    			    js2_dx2 = js2_dx1 - js2_dx0;
		    //    		    }
		    //    		    js2_dx0=js2_dx1;
		    //    	    }
		    //        }
		    //}

		    //k = start + (size - 1) / 2;
		    //splitvalue = (this->*medianfunc)(splitdim, k, start, end, otp, irearrangeandbalance);

	    }

	    //Now Split the node
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
			    left = BuildNodes_CRIT(start, k+1, newotp[0], param);
			    #pragma omp task
			    right = BuildNodes_CRIT(k+1, end, newotp[1], param);
			    #pragma omp taskwait
		    }

		    //Now save the largest distance from the center of this node
		    //Here, the center is just defined as the mean coordinates of the included particles
		    //TO DO LIST
		    //	Define center by using the 'finding the smallest sphere' algorithm
		    //	(https://en.wikipedia.org/wiki/Smallest-circle_problem)

		    //Left
		    Double_t js_center[ND], js_centertmp, js_pos[3], js_vel[3];
		    Double_t js_cenPos[3], js_cenVel[3];
		    Double_t js_dd=-1, js_dd2;

		    //Left - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(Int_t js_i=start; js_i<k+1; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (k+1 - start));
		    }

		    for(int js_i=0; js_i<ND; js_i++){
			    left->SetCenter(js_center[js_i], js_i);
			    js_cenPos[js_i] = js_center[js_i];
			    if(js_i>=3) js_cenVel[js_i-3] = js_center[js_i];
		    }

		    //Left - Largest Distance
		    for(Int_t js_i=start; js_i<k+1; js_i++){
			    js_dd2 = 0.;
			    for(int js_j=0; js_j<3; js_j++) js_pos[js_j] = bucket[js_i].GetPosition(js_j);
			    if(ND==6) for(int js_j=0; js_j<3; js_j++) js_vel[js_j] = bucket[js_i].GetVelocity(js_j);
			    js_dd2 += DistanceSqd(js_pos, js_cenPos, 3)/param[1];
			    if(ND==6) js_dd2 += DistanceSqd(js_vel, js_cenVel, 3)/param[2];

			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    left->SetFarthest(js_dd);

		    //Right
		    js_dd = -1.;
		    //Right - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(Int_t js_i=k+1; js_i<end; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (end - (k+1)));
		    }

		    for(int js_i=0; js_i<ND; js_i++){
			    right->SetCenter(js_center[js_i], js_i);
			    js_cenPos[js_i] = js_center[js_i];
			    if(js_i>=3) js_cenVel[js_i-3] = js_center[js_i];
		    }

		    //Right - largest distance
		    for(Int_t js_i=k+1; js_i<end; js_i++){
			    js_dd2 = 0.;
			    for(int js_j=0; js_j<3; js_j++) js_pos[js_j] = bucket[js_i].GetPosition(js_j);
			    if(ND==6) for(int js_j=0; js_j<3; js_j++) js_vel[js_j] = bucket[js_i].GetVelocity(js_j);

			    js_dd2 += DistanceSqd(js_pos, js_cenPos, 3)/param[1];
			    if(ND==6) js_dd2 += DistanceSqd(js_vel, js_cenVel, 3)/param[2];

			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    right->SetFarthest(js_dd);

		    return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);

#endif
	    }
	    else {
		    Node *left, *right;

		    left = BuildNodes_CRIT(start, k+1, otp, param);
		    right = BuildNodes_CRIT(k+1, end, otp, param);

		    //Now save the largest distance from the center of this node
		    //Here, the center is just defined as the mean coordinates of the included particles
		    //TO DO LIST
		    //	Define center by using the 'finding the smallest sphere' algorithm
		    //	(https://en.wikipedia.org/wiki/Smallest-circle_problem)

		    //Left
		    Double_t js_center[ND], js_centertmp, js_pos[3], js_vel[3];
		    Double_t js_cenPos[3], js_cenVel[3];
		    Double_t js_dd=-1, js_dd2;

		    //Left - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(Int_t js_i=start; js_i<k+1; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (k+1 - start));
		    }

		    for(int js_i=0; js_i<ND; js_i++){
			    left->SetCenter(js_center[js_i], js_i);
			    js_cenPos[js_i] = js_center[js_i];
			    if(ND==6) js_cenVel[js_i-3] = js_center[js_i];
		    }

		    //Left - Largest Distance
		    for(Int_t js_i=start; js_i<k+1; js_i++){
			    js_dd2 = 0.;
			    for(int js_j=0; js_j<3; js_j++) js_pos[js_j] = bucket[js_i].GetPosition(js_j);
			    if(ND==6) for(int js_j=0; js_j<3; js_j++) js_vel[js_j] = bucket[js_i].GetVelocity(js_j);
			    js_dd2 += DistanceSqd(js_pos, js_cenPos, 3)/param[1];
			    if(ND==6) js_dd2 += DistanceSqd(js_vel, js_cenVel, 3)/param[2];

			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    left->SetFarthest(js_dd);

		    //Right
		    js_dd = -1.;
		    //Right - Define Center
		    for(int js_j=0; js_j<ND; js_j++){
			    js_centertmp=0.;
			    for(Int_t js_i=k+1; js_i<end; js_i++) js_centertmp += bucket[js_i].GetPhase(js_j);
			    js_center[js_j] = js_centertmp / (Double_t (end - (k+1)));
		    }

		    for(int js_i=0; js_i<ND; js_i++){
			    right->SetCenter(js_center[js_i], js_i);
			    js_cenPos[js_i] = js_center[js_i];
			    if(ND==6) js_cenVel[js_i-3] = js_center[js_i];
		    }

		    //Right - largest distance
		    for(Int_t js_i=k+1; js_i<end; js_i++){
			    js_dd2 = 0.;
			    for(int js_j=0; js_j<3; js_j++) js_pos[js_j] = bucket[js_i].GetPosition(js_j);
			    if(ND==6) for(int js_j=0; js_j<3; js_j++) js_vel[js_j] = bucket[js_i].GetVelocity(js_j);

			    js_dd2 += DistanceSqd(js_pos, js_cenPos, 3)/param[1];
			    if(ND==6) js_dd2 += DistanceSqd(js_vel, js_cenVel, 3)/param[2];

			    if(js_dd2 > js_dd) js_dd = js_dd2;
		    }
		    right->SetFarthest(js_dd);


		    return new SplitNode(id, splitdim, splitvalue, size, bnd, start, end, ND, left, right);
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
            medianfunc=&NBody::KDTree::MedianPos;
        }
        else if (treetype==TVEL) {
            bmfunc=&NBody::KDTree::BoundaryandMeanVel;
            dispfunc=&NBody::KDTree::DispersionVel;
            spreadfunc=&NBody::KDTree::SpreadestVel;
            entropyfunc=&NBody::KDTree::EntropyVel;
            medianfunc=&NBody::KDTree::MedianVel;
        }
        else if (treetype==TPHS) {
            bmfunc=&NBody::KDTree::BoundaryandMeanPhs;
            dispfunc=&NBody::KDTree::DispersionPhs;
            spreadfunc=&NBody::KDTree::SpreadestPhs;
            entropyfunc=&NBody::KDTree::EntropyPhs;
            medianfunc=&NBody::KDTree::MedianPhs;
        }
        else if (treetype==TPROJ) {
            bmfunc=&NBody::KDTree::BoundaryandMeanPos;
            dispfunc=&NBody::KDTree::DispersionPos;
            spreadfunc=&NBody::KDTree::SpreadestPos;
            entropyfunc=&NBody::KDTree::EntropyPos;
            medianfunc=&NBody::KDTree::MedianPos;
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
        //if (node->GetCount() <= b) {
	if(node->GetLeaf() > 0){
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
        //if (node->GetCount() > b) {
	if(node->GetLeaf() < 0){
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
    //KDTree for OMP building (specified by the first Double_t)
    KDTree::KDTree(Double_t rdist, Particle *p, Int_t nparts, Int_t bucket_size,
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
	js_rdist = rdist;
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
            root=BuildNodes_OMP(0,numparts, otp);
            if (ibuildinparallel) BuildNodeIDs();
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
            if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_nested(inested);
#endif
    }

    //KDTree for normal adaptive KDTree building (specified by the first interger)
    KDTree::KDTree(Int_t js_int, Particle *p, Int_t nparts, Int_t bucket_size,
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
            root=BuildNodes_ADT(0,numparts, otp);
	    root->SetFarthest(1e31);
	    for(int js_i=0; js_i<ND; js_i++) root->SetCenter(0., js_i);

            if (ibuildinparallel) BuildNodeIDs();
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
            if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
#ifdef USEOPENMP
        omp_set_nested(inested);
#endif
    }

    //KDTree for the FOFSearchCriterion routine (specified by the first two Double_t)
    KDTree::KDTree(Double_t rdist, Double_t *param, Particle *p, Int_t nparts, Int_t bucket_size,
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
	js_rdist = rdist;

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
            root=BuildNodes_CRIT(0,numparts, otp, param);
	    root->SetFarthest(1e31);
	    for(int js_i=0; js_i<ND; js_i++) root->SetCenter(0., js_i);

            if (ibuildinparallel) BuildNodeIDs();
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
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
            if (iresetorder) qsort(bucket, numparts, sizeof(Particle), IDCompare);
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
        for (auto i=0;i<ompthreadpool.nactivethreads;i++) ompthreadpool.activethreadids[i]=i;
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
            for (auto i=0;i<newthreadpool[0].nactivethreads;i++)
                newthreadpool[0].activethreadids[i] = ompthreadpool.activethreadids[i];
            for (auto i=0;i<newthreadpool[1].nactivethreads;i++)
                newthreadpool[1].activethreadids[i] = ompthreadpool.activethreadids[i+newthreadpool[0].nactivethreads];
        }
#endif
        return newthreadpool;
    }
}
