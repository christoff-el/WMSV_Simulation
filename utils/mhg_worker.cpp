#ifndef MHG_WORKER_CPP
#define MHG_WORKER_CPP 1

/* MEX function [s,coef]=mhg([MAX,K],alpha,p,q,x,y)
 computes the truncated hypergeometric function pFq ^alpha(p;q;x;y)
 of one or two matrix arguments.
 
 The sum is only over those partitions kappa such that 
 a) |kappa|<=MAX 
 b) kappa has no more than n=length(x) parts
 c) kappa_1<=K (only if K is specified, if K is omitted, this restriction 
    is not applied)
 
 p and q are arrays, so mhg(30,9,[3 4],[5 6 7],[0.5 0.6],[0.8,0.9]) is 
 2F3^9([3 4],[5 6 7];[0.5 0.6], [0.8 0.9]) summed over all kappa with 
 |kappa|<=30

 K and y may be omitted. If y is omitted, the hypergeometric function of one 
 matrix argument is computed.

 Copyright, Plamen Koev, Massachusetts Institute of Technology
 Written: May 2004. 
 
 Updated:
 
 December 2005: introduced coef -- array from 0 to MAX, such that 
                coef[i]=sum over all partitions of size i. 
                Then c=sum_{i=0}^MAX coef[i].
 
 May 2007:      The Jacks are now Schur normalized leading to 
                   easier updates, i.e., Sx_kappa=Jx_kappa/H^*_kappa
                   where H^_kappa is the product of upper hook lengths
                Introduced Raymond Kan's faster Q_kappa and 
                   beta_{lambda mu} updates; 
                Removed all recursions; 
                Introduced a faster update for partitions with 
                   exactly n parts from Kaneko (1993) eq. (5.3);

 June 2007 ideas to consider: 
        1. malloc not more than choose(n+K,n) space if K is specified
        2. can we get J^(1/alpha) from J^(alpha), Stanley 89, Cor. 3.5<
        
 April 2014 Ported to C++ by Christopher Davis.
*/


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector>

using namespace std;

template<typename SCAL>
SCAL
mhg_worker(const int MAX, double alpha,
			const int np, const SCAL *p, 
			const int nq, const SCAL *q, 
			const int n, const SCAL *x) {

	alpha=1.0;
   	int    			nmu, i, j, k;
   	vector<int>		f(MAX+2);
	vector<SCAL>	xn((n+1)*(MAX+2));
	vector<SCAL>	prodx(n+1);

   	SCAL 		 c, s, zn, dn, t, q1, q2;
   	vector<SCAL> z(n+1,1);
   	vector<SCAL> coef(MAX+1,0);
   	vector<SCAL> kt(n+1);
   	vector<SCAL> mt(n+1);
   	vector<SCAL> blm(n+1);


   	int 		lg, w, slm, nhstrip, gz;
   	vector<int> l(n+1);
   	vector<int> mu(n+1);
   	vector<int> d(n);
   	vector<int> g(n+1);
   	vector<int> ww(n+1,1);
	vector<int> lmd(n+1);
	
	//for (i=1;i<=MAX;i++) coef[i]=0;  /* set to zero, these are the coefficients of the polynomial */ 
  	coef[0]=1; /* free term equals one */
  
  	w = 0; /* index of the zero partition, currently l*/

  	/* figure out the number of partitions |kappa|<= MAX with at most n parts */

  	for (i=1;i<=MAX+1;i++) f[i]=i;
	
  	for (i=2;i<n;i++) for (j=i+1;j<=MAX+1;j++) f[j]+=f[j-i];

  	w=f[MAX+1];
  	
  	vector<int> 	D(w+1);
  	vector<SCAL> 	Sx(n*(w+1));

  	prodx[1]=x[0];
  	for (i=2;i<=n;i++) prodx[i]=prodx[i-1]*x[i-1];
  	for (i=1; i<=n; i++) {
    	Sx[n+i-1]=1;
    	xn[(MAX+2)*i+1]=1;
    	for (j=2;j<=MAX+1;j++) xn[(MAX+2)*i+j]=xn[(MAX+2)*i+j-1]*x[i-1];
  	}

  
     /* l[0] is the second element of MAX, when MAX has 2 elements */
	l[0]=MAX; 
  	/* this is what limits l[1] by the second element of MAX if needed and 
    	 allows for the check l[i]<l[i-1] to be OK even for i=1 */
  
  	//for (i=1;i<=n;i++) z[i]=1;
  
  	for (i=1;i<=n;i++) kt[i]=-i;

  	//for (i=1;i<=n;i++) ww[i]=1;

  	int heap = MAX+2;
  	SCAL cc=1;
  	int h=1;
  	int sl=1;  /* sl= sum(l) */
  
  	while (h>0) {
      	if ((l[h]<l[h-1]) && (MAX>=sl) && (abs(z[h])!=0)) {
          
          	l[h]++;
          
          	if ((l[h]==1) && (h>1) && (h<n)) {
              	D[ww[h]]=heap;
              	ww[h]=heap;
              	k=MAX-sl+l[h];
              	if (k>l[h-1]) k=l[h-1];
              	heap+=k;
          	}
          	else ww[h]++;
          	w=ww[h];
              
          	/* Update Q */
          	c=(1-h)/alpha+l[h]-1;
          	zn=alpha;
          	dn=kt[h]+(double)h+1.0;
          	for (j=0;j<np;j++)  zn*=p[j]+c;
          	for (j=0;j<nq;j++)  dn*=q[j]+c;
	
          	kt[h]+=alpha;
          	for (j=1;j<h;j++) {
              	t=kt[j]-kt[h];
              	zn*=t;
              	dn*=t+1.0;
          	}
          	z[h]*=zn/dn;
          
          	/* Working hard only when l has less than n parts */
          
          
          	if (h<n) {
              	t=h+1.0-alpha; cc=1; for (j=1; j<=h;j++) cc*=(t+kt[j])/((double)h+kt[j]); 
              
              	/* computing the index of l-ones(1,h) */
              	nmu=l[1]; k=2; while ((k<=h)&&(l[k]>1)) nmu=D[nmu]+l[k++]-2;
                
              	Sx[w*n+h-1]=cc*prodx[h]*Sx[nmu*n+h-1];

	            cc=1; /* this way we can just update from 1 in the h=n case*/

              	d[h-1]--; /* technically this has to execute only when h>1 
                	           but is OK if it is always executed; d[0] will 
                    	       end up being -MAX at the end of the code */

              	d[h]=l[h];  /* for (k=1;k<h;k++) d[k]=l[k]-l[k+1]; 
                	             this happens automatically now via updates */
              
              	lg=0; for (k=1;k<=h;k++) if (d[k]>0) {lg++; g[lg]=k;}
              	slm=1; /* this is sum(l-mu) */
              	nhstrip=1; for (k=1;k<=lg;k++) nhstrip*=d[g[k]]+1; nhstrip--;

              	memcpy(&mu[1],&l[1],sizeof(int)*h);
              	memcpy(&mt[1],&kt[1],sizeof(SCAL)*h);
              	for (k=1;k<=lg;k++) { blm[k]=1; lmd[k]=l[g[k]]-d[g[k]]; }
              
              	for (i=1;i<=nhstrip;i++) {
                  	j=lg;
                  	gz=g[lg];
                  	while (mu[gz]==lmd[j]) {
                      	mu[gz]=l[gz];
                      	mt[gz]=kt[gz];
                      	slm-=d[gz];
                      	j--;
                      	gz=g[j];
                  	}
                  	t=kt[gz]-mt[gz];
                  
                  	zn=1.0+t;
                  	dn=t+alpha;
                  	for (k=1; k<gz; k++) {
                      	q1=mt[k]-mt[gz];
                      	q2=kt[k]-mt[gz];
                      	zn*=(alpha-1.0+q1)*(1.0+q2);
                      	dn*=q1*(alpha+q2); 
                  	}
                  	blm[j]*=zn/dn;
                                    
                  	mu[gz]--;
                  	mt[gz]-=alpha;
                  	slm++;

                  	for (k=j+1;k<=lg;k++) blm[k]=blm[j];    
                  
                  	/* next, find the index of mu */
                  	nmu=mu[1]+1; for (k=2;k<=h-(mu[h]==0);k++) nmu=D[nmu]+mu[k]-1;
                  
                  	for (k=h+1; k<=n;k++) 
                      	Sx[w*n+k-1]+=blm[j]*Sx[nmu*n+k-2]*xn[k*(MAX+2)+slm];
                  
             	}
             
             	for (k=h; k<n; k++) Sx[w*n+k]+=Sx[w*n+k-1];
              	coef[sl]+=z[h]*Sx[w*n+n-1];

          	} /* of "if h<n" */
          	else {
              /* computing the index of the partition l-l[n]*ones(1,n) */
              	nmu=l[1]-l[n]+1;
              	k=2; while ((k<n)&&(l[k]>l[n])) nmu=D[nmu]+l[k++]-1-l[n];
              	/* cc is 1 if l[n]==1, (guaranteed by the h<n case); 
                	 we then update from the previous */ 
            
              	cc*=(1.0/alpha+l[n]-1.0)*complex<double>(real(prodx[n])/l[n],imag(prodx[n])/l[n]);
                for (k=1;k<n;k++) cc*=(1.0+kt[k]-kt[n])/(alpha+kt[k]-kt[n]);
                coef[sl]+=z[n]*cc*Sx[nmu*n+n-1];
          	}
          	if (h<n) {
              	z[h+1]=z[h];
              	h++;
              	ww[h]=w;
          	}
          	sl++;
      	}
      	else { 
          	sl-=l[h];
          	l[h]=0;
          	kt[h]=-h;
          	h--;
      	}
  	} /* of while h>0 */
  
  	s=0; 
  	for (i=0;i<MAX+1;i++) {
  		s+=coef[i];
  	}
	
	/*delete[] lmd;
	delete[] blm;
	delete[] mt;
	delete[] g;
	delete[] d;
	delete[] ww;
	delete[] kt;
	delete[] mu;
	delete[] z;
	delete[] l;
	delete[] prodx;
	delete[] xn;
	delete[] Sx;
	delete[] D;
	delete[] coef;*/
	
	return s;
	
}

#endif	//MHG_WORKER_CPP