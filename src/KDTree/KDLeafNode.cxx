/*! \file KDLeafNode.cxx
 *  \brief This file contains implementation of Leaf Node functions

*/

#include <iostream>
#include <KDNode.h>

namespace NBody
{
    ///\name Leaf Node Functions
    //@{
    ///\name Non-periodic calls
    //@{
    void LeafNode::FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
            if (dist2 < pq->TopPriority() && dist2 > 0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = VelDistSqd(bucket[target].GetVelocity(),bucket[i].GetVelocity(), dim);
            if (dist2 < pq->TopPriority() && dist2 >0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = PhaseDistSqd(bucket[target].GetPosition(),bucket[i].GetPosition(),
                bucket[target].GetVelocity(),bucket[i].GetVelocity());
            if (dist2 < pq->TopPriority() && dist2 > 0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, Double_t *metric)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = MetricDistSqd(bucket[target].GetPosition(),bucket[i].GetPosition(),
                bucket[target].GetVelocity(),bucket[i].GetVelocity(),metric);
            if (dist2 < pq->TopPriority() && dist2 > 0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = MetricwithTensorDistSqd(bucket[target].GetPosition(),bucket[i].GetPosition(),
                bucket[target].GetVelocity(),bucket[i].GetVelocity(),m0,m1,gm);
            if (dist2 < pq->TopPriority() && dist2 > 0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
            int cmpval=cmp(bucket[target],bucket[i],params);
            if (dist2 < pq->TopPriority() && dist2 > 0 && cmpval==1)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestCheck(Double_t rd, FOFcheckfunc check, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
            int checkval=check(bucket[i],params);
            if (dist2 < pq->TopPriority() && dist2 > 0 && checkval==0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }

    void LeafNode::FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = DistanceSqd(x,bucket[i].GetPosition(), dim);
            if (dist2 < pq->TopPriority())
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
        }
    }
    void LeafNode::FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *v, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = VelDistSqd(v,bucket[i].GetVelocity(), dim);
            if (dist2 < pq->TopPriority())
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
        }
    }
    void LeafNode::FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = PhaseDistSqd(x,bucket[i].GetPosition(),v,bucket[i].GetVelocity());
            if (dist2 < pq->TopPriority())
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
        }
    }
    void LeafNode::FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, int dim)
    {
        FindNearestPos(rd,bucket,pq,off,x.GetCoord(),dim);
    }
    void LeafNode::FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate v, int dim)
    {
        FindNearestVel(rd,bucket,pq,off,v.GetCoord(),dim);
    }
    void LeafNode::FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v)
    {
        FindNearestPhase(rd,bucket,pq,off,x.GetCoord(),v.GetCoord());
    }

    void LeafNode::FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *metric)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = MetricDistSqd(x,bucket[i].GetPosition(),
                v,bucket[i].GetVelocity(),metric);
            if (dist2 < pq->TopPriority() && dist2 > 0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
        }
    }
    void LeafNode::FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = MetricwithTensorDistSqd(x,bucket[i].GetPosition(),
                v,bucket[i].GetVelocity(),m0,m1,gm);
            if (dist2 < pq->TopPriority() && dist2 > 0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
        }
    }
    void LeafNode::FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *metric)
    {
        FindNearestMetric(rd,bucket,pq,off,x.GetCoord(),v.GetCoord(),metric);
    }
    void LeafNode::FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        FindNearestMetricwithTensor(rd,bucket,pq,off,x.GetCoord(),v.GetCoord(),m0,m1,gm);
    }
    void LeafNode::FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Particle &p, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (!(bucket[i]==p)){
            Double_t dist2 = DistanceSqd(p.GetPosition(),bucket[i].GetPosition(), dim);
            int cmpval=cmp(p,bucket[i],params);
            if (dist2 < pq->TopPriority() && dist2 > 0 && cmpval==1)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }
    void LeafNode::FindNearestCheck(Double_t rd, FOFcheckfunc check, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Particle &p, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (!(bucket[i]==p)){
            Double_t dist2 = DistanceSqd(p.GetPosition(),bucket[i].GetPosition(), dim);
            int checkval=check(bucket[i],params);
            if (dist2 < pq->TopPriority() && dist2 > 0 && checkval==0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
            }
        }
    }

    void LeafNode::FindNearestCheck(Double_t rd, FOFcheckfunc check, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate &x, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = DistanceSqd(x.GetCoord(),bucket[i].GetPosition(), dim);
            int checkval=check(bucket[i],params);
            if (dist2 < pq->TopPriority() && dist2 > 0 && checkval==0)
            {
                pq->Pop();
                pq->Push(i, dist2);
            }
        }
    }

    void LeafNode::SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t target, int dim)
    {
        //first check to see if entire node lies wihtin search distance
        Double_t maxr0=0.,maxr1=0.;
        for (int j=0;j<dim;j++){
            maxr0+=(bucket[target].GetPosition(j)-xbnd[j][0])*(bucket[target].GetPosition(j)-xbnd[j][0]);
            maxr1+=(bucket[target].GetPosition(j)-xbnd[j][1])*(bucket[target].GetPosition(j)-xbnd[j][1]);
        }
        if (maxr0<fdist2&&maxr1<fdist2)
            for (Int_t i = bucket_start; i < bucket_end; i++)
            {
                Int_t id=bucket[i].GetID();
                Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
                Group[id]=iGroup;
                pdist2[id]=dist2;
            }
        else
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
            if (dist2 < fdist2)
            {
                Int_t id=bucket[i].GetID();
                Group[id]=iGroup;
                pdist2[id]=dist2;
            }
            }
        }
    }

    void LeafNode::SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *x, int dim)
    {
        //first check to see if entire node lies wihtin search distance
        Double_t maxr0=0.,maxr1=0.;
        for (int j=0;j<dim;j++){
            maxr0+=(x[j]-xbnd[j][0])*(x[j]-xbnd[j][0]);
            maxr1+=(x[j]-xbnd[j][1])*(x[j]-xbnd[j][1]);
        }
        if (maxr0<fdist2&&maxr1<fdist2)
            for (Int_t i = bucket_start; i < bucket_end; i++)
            {
                Int_t id=bucket[i].GetID();
                Double_t dist2 = DistanceSqd(x,bucket[i].GetPosition(), dim);
                Group[id]=iGroup;
                pdist2[id]=dist2;
            }
        else
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = DistanceSqd(x,bucket[i].GetPosition(), dim);
            if (dist2 < fdist2)
            {
                Int_t id=bucket[i].GetID();
                Group[id]=iGroup;
                pdist2[id]=dist2;
            }
        }
    }
    void LeafNode::SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Coordinate x, int dim)
    {
        SearchBallPos(rd,fdist2,iGroup,bucket,Group,pdist2,off,x.GetCoord(),dim);
    }

    void LeafNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Int_t target, Int_t &nt, int dim)
    {
        //first check to see if entire node lies wihtin search distance
        Double_t maxr0=0.,maxr1=0.;
        for (int j=0;j<dim;j++){
            maxr0+=(bucket[target].GetPosition(j)-xbnd[j][0])*(bucket[target].GetPosition(j)-xbnd[j][0]);
            maxr1+=(bucket[target].GetPosition(j)-xbnd[j][1])*(bucket[target].GetPosition(j)-xbnd[j][1]);
        }
        if (maxr0<fdist2&&maxr1<fdist2)
            for (Int_t i = bucket_start; i < bucket_end; i++)
                tagged[nt++]=i;
        else
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
            if (dist2 < fdist2) tagged[nt++]=i;
            }
        }
    }

    void LeafNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *x, Int_t &nt, int dim)
    {
        //first check to see if entire node lies wihtin search distance
        Double_t maxr0=0.,maxr1=0.;
        for (int j=0;j<dim;j++){
            maxr0+=(x[j]-xbnd[j][0])*(x[j]-xbnd[j][0]);
            maxr1+=(x[j]-xbnd[j][1])*(x[j]-xbnd[j][1]);
        }
        if (maxr0<fdist2&&maxr1<fdist2)
            for (Int_t i = bucket_start; i < bucket_end; i++)
                tagged[nt++]=i;
        else
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = DistanceSqd(x,bucket[i].GetPosition(), dim);
            if (dist2 < fdist2)
            {
                tagged[nt++]=i;
            }
        }
    }
    void LeafNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Coordinate x, Int_t &nt, int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x.GetCoord(),nt,dim);
    }

    void LeafNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, vector<Int_t> &tagged, Double_t* off, Int_t target, int dim)
    {
        //first check to see if entire node lies wihtin search distance
        Double_t maxr0=0.,maxr1=0.;
        for (int j=0;j<dim;j++){
            maxr0+=(bucket[target].GetPosition(j)-xbnd[j][0])*(bucket[target].GetPosition(j)-xbnd[j][0]);
            maxr1+=(bucket[target].GetPosition(j)-xbnd[j][1])*(bucket[target].GetPosition(j)-xbnd[j][1]);
        }
        if (maxr0<fdist2&&maxr1<fdist2)
            for (Int_t i = bucket_start; i < bucket_end; i++)
                tagged.push_back(i);
        else
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
            if (dist2 < fdist2) tagged.push_back(i);
            }
        }
    }

    void LeafNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, vector<Int_t> &tagged, Double_t* off, Double_t *x, int dim)
    {
        //first check to see if entire node lies wihtin search distance
        Double_t maxr0=0.,maxr1=0.;
        for (int j=0;j<dim;j++){
            maxr0+=(x[j]-xbnd[j][0])*(x[j]-xbnd[j][0]);
            maxr1+=(x[j]-xbnd[j][1])*(x[j]-xbnd[j][1]);
        }
        if (maxr0<fdist2&&maxr1<fdist2)
            for (Int_t i = bucket_start; i < bucket_end; i++)
                tagged.push_back(i);
        else
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            Double_t dist2 = DistanceSqd(x,bucket[i].GetPosition(), dim);
            if (dist2 < fdist2)
            {
                tagged.push_back(i);
            }
        }
    }
    void LeafNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, vector<Int_t> &tagged, Double_t* off, Coordinate x, int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x.GetCoord(),dim);
    }

    void LeafNode::SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target&&(Group[bucket[i].GetID()]>iGroup||Group[bucket[i].GetID()]==0)){
//            if (i!=target&&Group[bucket[i].GetID()]!=iGroup){
            if (cmp(bucket[target],bucket[i],params))
            {
                Double_t dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition(), dim);
                Int_t id=bucket[i].GetID();
                Group[id]=iGroup;
                pdist2[id]=dist2;
            }
            }
        }
    }
    void LeafNode::SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Particle &target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (!(bucket[i]==target)&&(Group[bucket[i].GetID()]>iGroup||Group[bucket[i].GetID()]==0)){
            if (cmp(target,bucket[i],params))
            {
                Double_t dist2 = DistanceSqd(target.GetPosition(),bucket[i].GetPosition(), dim);
                Int_t id=bucket[i].GetID();
                Group[id]=iGroup;
                pdist2[id]=dist2;
            }
            }
        }
    }
    void LeafNode::SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged,Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            if (cmp(bucket[target],bucket[i],params))
            {
                tagged[nt++]=i;
            }
            }
        }
    }
    void LeafNode::SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Particle &target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (!(bucket[i]==target)){
            if (cmp(target,bucket[i],params))
            {
                tagged[nt++]=i;
            }
            }
        }
    }
    void LeafNode::SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, vector<Int_t> &tagged,Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target){
            if (cmp(bucket[target],bucket[i],params))
            {
                tagged.push_back(i);
            }
            }
        }
    }
    void LeafNode::SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, vector<Int_t> &tagged, Double_t* off, Particle &target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (!(bucket[i]==target)){
            if (cmp(target,bucket[i],params))
            {
                tagged.push_back(i);
            }
            }
        }
    }

    void LeafNode::SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Int_t target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (i!=target) if (Group[bucket[i].GetID()]>-1)
            if ((Group[bucket[i].GetID()]>iGroup||Group[bucket[i].GetID()]==0)){
            if (cmp(bucket[target],bucket[i],params))
//            if (i!=target&&(Group[bucket[i].GetID()]>iGroup||Group[bucket[i].GetID()]==0)){
//            if (cmp(bucket[target],bucket[i],params))
            {
                Int_t id=bucket[i].GetID();
                Group[id]=iGroup;
            }
            }
        }
    }
    void LeafNode::SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Particle &target, int dim)
    {
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (!(bucket[i]==target)) if (Group[bucket[i].GetID()]>-1)
            if ((Group[bucket[i].GetID()]>iGroup||Group[bucket[i].GetID()]==0)){
            if (cmp(target,bucket[i],params))
//            if (!(bucket[i]==target)&&(Group[bucket[i].GetID()]>iGroup||Group[bucket[i].GetID()]==0)){
//            if (cmp(target,bucket[i],params))
            {
                Int_t id=bucket[i].GetID();
                Group[id]=iGroup;
            }
            }
        }
    }

    void LeafNode::FOFSearchBall(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target)
    {
        //if bucket already linked and particle already part of group, do nothing.
        //if(BucketFlag[nid]&&Head[target]==Head[bucket_start])return;
	if(BucketFlag[nid])return;
        //this flag is initialized to !=0 and if entire bucket searched and all particles already linked,
        //then BucketFlag[nid]=1
        //int flag=Head[bucket_start];
	int flag=1;

	//--JS--
	//Node Skip by using trigonometric inequalities only acts when numdim==3
	Double_t js_pos[6], js_dist, js_rr;
	for(int js_j=0; js_j<numdim; js_j++) js_pos[js_j] = bucket[target].GetPhase(js_j);
	js_dist = DistanceSqd(js_pos, js_center, numdim);
	js_rr = js_farthest;

	if(sqrt(js_dist) >= sqrt(js_rr) + sqrt(fdist2)){
		flag=0;
	}
	else if(sqrt(js_dist) <= abs(sqrt(js_rr) - sqrt(fdist2)) && fdist2 > js_rr){
		//The the entire node lies within search distance
		Int_t id;
		for (Int_t i = bucket_start; i < bucket_end; i++){
			id=bucket[i].GetID();
			if (Group[id]) continue;
			Group[id]=iGroup;
			Fifo[iTail++]=i;
			Len[iGroup]++;

			Next[Tail[Head[target]]]=Head[i];
			Tail[Head[target]]=Tail[Head[i]];
			Head[i]=Head[target];

			if(iTail==nActive)iTail=0;
		}
	}
	else{
                Int_t id;
                Double_t dist2;
                for (Int_t i = bucket_start; i < bucket_end; i++)
                {
                    //if (flag!=Head[i])flag=0;
                    id=bucket[i].GetID();
                    if (Group[id]) continue;
                    dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition());
                    if (numdim==6) dist2+=DistanceSqd(bucket[target].GetVelocity(),bucket[i].GetVelocity());

                    //if (flag!=Head[i])flag=0;
		    //flag=0;
			if(Group[id]==0)flag=0;

                    if (dist2 < fdist2) {
                        Group[id]=iGroup;
                        Fifo[iTail++]=i;
                        Len[iGroup]++;

                        Next[Tail[Head[target]]]=Head[i];
                        Tail[Head[target]]=Tail[Head[i]];
                        Head[i]=Head[target];

                        if(iTail==nActive)iTail=0;
                        flag=0;
                    }
                }
	}
		//Otherwise check each particle individually

        //Double_t maxr0=0.,maxr1=0.;
        //for (int j=0;j<numdim;j++){
        //    maxr0+=(bucket[target].GetPhase(j)-xbnd[j][0])*(bucket[target].GetPhase(j)-xbnd[j][0]);
        //    maxr1+=(bucket[target].GetPhase(j)-xbnd[j][1])*(bucket[target].GetPhase(j)-xbnd[j][1]);
        //}
        //first check to see if entire node lies wihtin search distance
        //if (maxr0<fdist2&&maxr1<fdist2){
        //    Int_t id;
        //    for (Int_t i = bucket_start; i < bucket_end; i++){
        //        id=bucket[i].GetID();
        //        if (Group[id]) continue;
        //        Group[id]=iGroup;
        //        Fifo[iTail++]=i;
        //        Len[iGroup]++;

        //        Next[Tail[Head[target]]]=Head[i];
        //        Tail[Head[target]]=Tail[Head[i]];
        //        Head[i]=Head[target];

        //        if(iTail==nActive)iTail=0;
        //    }
        //}
        ////otherwise check each particle individually
        //else {
        //    Int_t id;
        //    Double_t dist2;
        //    for (Int_t i = bucket_start; i < bucket_end; i++)
        //    {
        //        if (flag!=Head[i])flag=0;
        //        id=bucket[i].GetID();
        //        if (Group[id]) continue;
        //        dist2 = DistanceSqd(bucket[target].GetPosition(),bucket[i].GetPosition());
        //        if (numdim==6) dist2+=DistanceSqd(bucket[target].GetVelocity(),bucket[i].GetVelocity());
        //        if (dist2 < fdist2) {
        //            Group[id]=iGroup;
        //            Fifo[iTail++]=i;
        //            Len[iGroup]++;

        //            Next[Tail[Head[target]]]=Head[i];
        //            Tail[Head[target]]=Tail[Head[i]];
        //            Head[i]=Head[target];

        //            if(iTail==nActive)iTail=0;
        //            flag=0;
        //        }
        //    }
        //}
        if (flag) {
		BucketFlag[nid]=1;
		if(BucketFlag[sibling->GetID()]==1)BucketFlag[parent->GetID()]=1;
	}
    }
    void LeafNode::FOFSearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target)
    {
        //if bucket already linked and particle already part of group, do nothing.
        //if(BucketFlag[nid]&&Head[target]==Head[bucket_start])return;
	if(BucketFlag[nid])return;
        //this flag is initialized to !=0 and if entire bucket searched and all particles already linked,
        //then BucketFlag[nid]=1
	//int flag=Head[bucket_start];
	int flag=1;


	Double_t js_pos[3], js_vel[3], js_dist=0., js_rr;
	Double_t js_posCen[3], js_velCen[3];
	for(int js_j=0; js_j<3; js_j++) {js_pos[js_j] = bucket[target].GetPosition(js_j); js_posCen[js_j] = js_center[js_j];}
	if(numdim==6) for(int js_j=3; js_j<6; js_j++) {js_vel[js_j-3] = bucket[target].GetVelocity(js_j-3); js_velCen[js_j-3] = js_center[js_j];}

	js_dist += DistanceSqd(js_pos, js_posCen, 3)/params[1];
	if(numdim==6) js_dist += DistanceSqd(js_vel, js_velCen, 3)/params[2];
	js_rr = js_farthest;

	if(sqrt(js_dist) >= sqrt(js_rr) + 1.0){
		flag=0;
	}
	else if(sqrt(js_dist) <= abs(sqrt(js_rr) - 1.0) && 1.0 > js_rr){
		for(Int_t i=bucket_start; i < bucket_end; i++){
			Int_t id = bucket[i].GetID();
			//if(Group[id]==iGroup) continue;
			if(Group[id]<0) continue;
			if(Group[id]) continue;

			Group[id]=iGroup;
			Fifo[iTail++]=i;
			Len[iGroup]++;

			Next[Tail[Head[target]]]=Head[i];
			Tail[Head[target]]=Tail[Head[i]];
			Head[i]=Head[target];
			if(iTail==nActive)iTail=0;
		}
	}
	else{
        	for (Int_t i = bucket_start; i < bucket_end; i++)
        	{
        	    //if (flag!=Head[i])flag=0;
        	    Int_t id=bucket[i].GetID();
        	    //if already linked don't do anything
        	    //if (Group[id]==iGroup) continue;
		    if (Group[id]) continue;
        	    //if tag below zero then don't do anything
        	    if (Group[id]<0) continue;

        	    //if (flag!=Head[i])flag=0;
		    flag=0;

        	    if (cmp(bucket[target],bucket[i],params)) {
        	        Group[id]=iGroup;
        	        Fifo[iTail++]=i;
        	        Len[iGroup]++;

        	        Next[Tail[Head[target]]]=Head[i];
        	        Tail[Head[target]]=Tail[Head[i]];
        	        Head[i]=Head[target];
        	        if(iTail==nActive)iTail=0;
        	        flag=0;
        	    }
        	}
	}
	
        //for (Int_t i = bucket_start; i < bucket_end; i++)
        //{
        //    if (flag!=Head[i])flag=0;
        //    Int_t id=bucket[i].GetID();
        //    //if already linked don't do anything
        //    //if (Group[id]==iGroup) continue;
	//    if (Group[id]) continue;
        //    //if tag below zero then don't do anything
        //    if (Group[id]<0) continue;
        //    if (cmp(bucket[target],bucket[i],params)) {
        //        Group[id]=iGroup;
        //        Fifo[iTail++]=i;
        //        Len[iGroup]++;

        //        Next[Tail[Head[target]]]=Head[i];
        //        Tail[Head[target]]=Tail[Head[i]];
        //        Head[i]=Head[target];
        //        if(iTail==nActive)iTail=0;
        //        flag=0;
        //    }
        //}
        if (flag) BucketFlag[nid]=1;
    }
    void LeafNode::FOFSearchCriterionSetBasisForLinks(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target)
    {
        //if bucket already linked and particle already part of group, do nothing.
        if(BucketFlag[nid]&&Head[target]==Head[bucket_start])return;
        //this flag is initialized to !=0 and if entire bucket searched and all particles already linked,
        //then BucketFlag[nid]=1
        int flag=Head[bucket_start];
        for (Int_t i = bucket_start; i < bucket_end; i++)
        {
            if (flag!=Head[i])flag=0;
            Int_t id=bucket[i].GetID();
            //if already linked don't do anything
            //if (Group[id]==iGroup) continue;
	    if (Group[id]) continue;
            //if tag below zero then don't do anything
            if (Group[id]<0) continue;
            if (cmp(bucket[target],bucket[i],params)) {
                //also possible particle is tagged in another group and cannot be used for generating new links so this link invalid
                //so check if particle can be used for links by using check function and if result <0 then do nothing
                if (Group[id]!=iGroup && Group[id]>0) continue;

                Group[id]=iGroup;
                Fifo[iTail++]=i;
                Len[iGroup]++;

                Next[Tail[Head[target]]]=Head[i];
                Tail[Head[target]]=Tail[Head[i]];
                Head[i]=Head[target];
                if(iTail==nActive)iTail=0;
                flag=0;
            }
        }
        if (flag) BucketFlag[nid]=1;
    }
    //@}

    ///\name Periodic
    //@{
    //As leaf nodes should never be head node and split nodes account for periodicity, these are the same as the non-periodic case.
    void LeafNode::FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        FindNearestPos(rd,bucket,pq,off,target,dim);
    }
    void LeafNode::FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        FindNearestVel(rd,bucket,pq,off,target,dim);
    }
    void LeafNode::FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target)
    {
        Coordinate x0,v;
        x0=Coordinate(bucket[target].GetPosition());
        v=Coordinate(bucket[target].GetVelocity());
        FindNearestPhase(rd,bucket,pq,off,x0,v);
    }
    void LeafNode::FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, Double_t *metric)
    {
        Coordinate x0,v;
        x0=Coordinate(bucket[target].GetPosition());
        v=Coordinate(bucket[target].GetVelocity());
        FindNearestMetric(rd,bucket,pq,off,x0,v,metric);
    }
    void LeafNode::FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        Coordinate x0,v;
        x0=Coordinate(bucket[target].GetPosition());
        v=Coordinate(bucket[target].GetVelocity());
        FindNearestMetricwithTensor(rd,bucket,pq,off,x0,v,m0,m1,gm);
    }
    void LeafNode::FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle p0;
        p0=bucket[target];
        FindNearestCriterion(rd,cmp,params,bucket,pq,off,p0,dim);
    }
    void LeafNode::FindNearestCheckPeriodic(Double_t rd, FOFcheckfunc check, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle p0;
        p0=bucket[target];
        FindNearestCheck(rd,check,params,bucket,pq,off,p0,dim);
    }

    void LeafNode::FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, int dim)
    {
        FindNearestPos(rd,bucket,pq,off,x,dim);
    }
    void LeafNode::FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *v, int dim)
    {
        FindNearestVel(rd,bucket,pq,off,v,dim);
    }
    void LeafNode::FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, Double_t *v)
    {
        FindNearestPhase(rd,bucket,pq,off,x,v);
    }
    void LeafNode::FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, int dim)
    {
        FindNearestPosPeriodic(rd,bucket,pq,off,p,x.GetCoord(),dim);
    }
    void LeafNode::FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate v, int dim)
    {
        FindNearestVelPeriodic(rd,bucket,pq,off,p,v.GetCoord(),dim);
    }
    void LeafNode::FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, Coordinate v)
    {
        FindNearestPhasePeriodic(rd,bucket,pq,off,p,x.GetCoord(),v.GetCoord());
    }

    void LeafNode::FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, Double_t *v, Double_t *metric)
    {
        FindNearestMetric(rd,bucket,pq,off,x,v,metric);
    }
    void LeafNode::FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        FindNearestMetricwithTensor(rd,bucket,pq,off,x,v,m0,m1,gm);
    }
    void LeafNode::FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, Coordinate v, Double_t *metric)
    {
        FindNearestMetricPeriodic(rd,bucket,pq,off,p,x.GetCoord(),v.GetCoord(),metric);
    }
    void LeafNode::FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        FindNearestMetricwithTensorPeriodic(rd,bucket,pq,off,p,x.GetCoord(),v.GetCoord(),m0,m1,gm);
    }
    void LeafNode::FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        FindNearestCriterion(rd,cmp,params,bucket,pq,off,p0,dim);
    }
    void LeafNode::FindNearestCheckPeriodic(Double_t rd, FOFcheckfunc check, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        FindNearestCheck(rd,check,params,bucket,pq,off,p0,dim);
    }
    void LeafNode::FindNearestCheckPeriodic(Double_t rd, FOFcheckfunc check, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate &x0, int dim)
    {
        FindNearestCheck(rd,check,params,bucket,pq,off,x0,dim);
    }

    void LeafNode::SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Coordinate x0;
        x0=Coordinate(bucket[target].GetPosition());
        SearchBallPos(rd,fdist2,iGroup,bucket,Group,pdist2,off,x0,dim);
    }
    void LeafNode::SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *p, Double_t *x, int dim)
    {
        SearchBallPos(rd,fdist2,iGroup,bucket,Group,pdist2,off,x,dim);
    }
    void LeafNode::SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *p, Coordinate x, int dim)
    {
        SearchBallPosPeriodic(rd,fdist2,iGroup,bucket,Group,pdist2,off,p,x.GetCoord(),dim);
    }

    void LeafNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *p, Int_t target, Int_t &nt, int dim)
    {
        Coordinate x0;
        x0=Coordinate(bucket[target].GetPosition());
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x0,nt,dim);
    }
    void LeafNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2,  Particle *bucket, Int_t *tagged, Double_t *off, Double_t *p, Double_t *x, Int_t &nt, int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x,nt,dim);
    }
    void LeafNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *p, Coordinate x, Int_t &nt, int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x.GetCoord(),nt,dim);
    }

    void LeafNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, vector<Int_t> &tagged, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Coordinate x0;
        x0=Coordinate(bucket[target].GetPosition());
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x0,dim);
    }
    void LeafNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2,  Particle *bucket, vector<Int_t> &tagged, Double_t *off, Double_t *p, Double_t *x, int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x,dim);
    }
    void LeafNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, vector<Int_t> &tagged, Double_t *off, Double_t *p, Coordinate x, int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x.GetCoord(),dim);
    }

    void LeafNode::SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterion(rd,cmp,params,iGroup,bucket,Group,pdist2,off,x0,dim);
    }
    void LeafNode::SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,x0,dim);
    }
    void LeafNode::SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t *off, Double_t *p, Particle &target, int dim)
    {
        SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target,dim);
    }

    void LeafNode::SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, vector<Int_t> &tagged, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterionTagged(rd,cmp,params,bucket,tagged,off,x0,dim);
    }
    void LeafNode::SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, vector<Int_t> &tagged, Double_t *off, Double_t *p, Particle &target, int dim)
    {
        SearchCriterionTagged(rd,cmp,params,bucket,tagged,off,target,dim);
    }

    void LeafNode::SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,x0,dim);
    }
    void LeafNode::SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        SearchCriterion(rd,cmp,params,iGroup,bucket,Group,pdist2,off,p0,dim);
    }
    void LeafNode::SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,p0,dim);
    }
    void LeafNode::FOFSearchBallPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t *off, Double_t *p, Int_t target)
    {
        FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
    }

    void LeafNode::FOFSearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t *off, Double_t *p, Int_t target)
    {
        FOFSearchCriterion(rd, cmp, params, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
    }

    void LeafNode::FOFSearchCriterionSetBasisForLinksPeriodic(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t *off, Double_t *p, Int_t target)
    {
        FOFSearchCriterionSetBasisForLinks(rd, cmp, check, params, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
    }
    //@}

    //@}

}
