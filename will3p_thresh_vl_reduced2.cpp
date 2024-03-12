#include <Rcpp.h>
//[[Rcpp::depends(dqrng, BH, sitmo)]]
#include <xoshiro.h>
#include <dqrng_distribution.h>
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]
#include <omp.h>

using namespace Rcpp;

//[[Rcpp::depends(dqrng)]]

//[[Rcpp::export]]
NumericMatrix ct_func2_cpp_thresh(const NumericVector& a, const NumericVector& b, const NumericVector& c, const NumericVector& l, const NumericVector& thresh,
                            const NumericVector& p, const int nthreads, const NumericVector& pseeds, const int vg, const int vk, const int vt, double stoch) 
{
  dqrng::normal_distribution distA(a[0],a[1]); 
  dqrng::normal_distribution distB(b[0],b[1]);
  dqrng::normal_distribution distC(c[0],c[1]);
  dqrng::normal_distribution distL2(l[0],l[1]);
  dqrng::normal_distribution distP(thresh[0],thresh[1]);
  dqrng::uniform_distribution distT(-stoch, stoch);
  
  //double lower;
  //double upper;
  //if(stoch<0.1){
  //  lower=-stoch;
  //  upper=0;
  //} 
  //else{
  //  lower=-stoch;
  //  upper=stoch;
  //} 
  
  /* Set-up and some bound checking */
  double lag = thresh[2];
  double flat = l[0];
  
  double vli = p[1];
  if(vli<0) {
    Rprintf("Max CT is %lg, now reduced to 40\n",vli);
    vli=0;
  }
  double vl_max = p[2];
  if(vl_max>15) {
    Rprintf("Max VL is %lg, now decreased to 1\n",vl_max);
    vl_max=15;
  }
  if(vl_max<=vli) {
    Rprintf("Max vL is <= min VL, now adjusted\n",vl_max);
    vl_max=15;
    vli=0;
  }
  int ct_n=(int) 2*(vl_max-vli)+1;
  double t_max = p[0];
  if(t_max>55) {
    Rprintf("Max time is %lg, now reduced to 55\n",t_max);
    t_max=55;
  }		
  int t_n=(int) (t_max+1);
  
  double a_lower = 0.5;
  double b_lower = 0.25;
  
  int ncores=nthreads; /* need a copy to be able to change it */
  if(ncores>32) {
    Rprintf("ncores is %i, now reduced to 32\n",ncores);
    ncores=32;
  }	
  int test_pop = (int) p[3];
  if(test_pop<100) {
    Rprintf("Popiulation size is %i, now increased to 100\n",test_pop);
    test_pop=100;		
  }
  int n_s = test_pop/ncores; /* pop to simulate per thread */
  
  /* initialise thread specific random number generators with the seeds provided
   - neeed to *8 to avoid false-sharing (cache line conflicts */
  dqrng::xoroshiro128plus* rng = new dqrng::xoroshiro128plus[ncores*8];
  for (int i = 0; i < ncores*8; i++) rng[i] = dqrng::xoroshiro128plus(pseeds(i/8));
  
  omp_set_num_threads(ncores); /* set number of threads */
  
  /* static thread specific frequency matrix - need to zero */	 
  double ft[32][41][51];  /* max 32 threads, 41 CT categories, 32 timesteps */
  for(int tn=0;tn<ncores;tn++)
    for(int i=0;i<ct_n;i++)
      for(int j=0;j<t_n;j++)
        ft[tn][i][j]=0;
  
  /* The one parallel loop which does all the work */
#pragma omp parallel for schedule(static,1)
  for(int tn=0;tn<ncores;tn++) {
    int l=tn*n_s;
    int m=l+n_s;
    for(int k=l; k<m; k++){
      double asv;
      double bsv;
      double lsv;
      double tsv;
      double tstoch = distT(rng[tn*8]);
      
      if(flat==0) lsv=0;
      else lsv = exp(distL2(rng[tn*8]));
      
      //do {
      tsv = distP(rng[tn*8]);
      //} while (tsv<vli);
      
      double csv = distC(rng[tn*8]);
      double cv;
      
      if(vg==0){
        do {
          asv = exp(distA(rng[tn*8]));
        } while (asv<(csv-vli)/14);
        do {
          bsv = exp(distB(rng[tn*8]));
        } while (bsv<(csv-vli)/28);
      }
      else if(vg==1){
        csv = distC(rng[tn*8]);
        do {
          asv = exp(distA(rng[tn*8]));
          bsv = exp(distB(rng[tn*8]))/asv;
        } while (asv<0.5 && bsv<0.25);
      }
      else if(vg==2){
        do {
          asv = exp(distA(rng[tn*8]));
          csv = distC(rng[tn*8])*asv;
          bsv = exp(distB(rng[tn*8]))/asv;
        } while (asv<(csv-vli)/14 && bsv<(csv-vli)/28);
      }
      else if(vg==3){
        csv = distC(rng[tn*8]);
        asv = exp(distA(rng[tn*8]));
        bsv = exp(distB(rng[tn*8]));
      }
      else if(vg==4){
        csv = distC(rng[tn*8]);
        asv = exp(distA(rng[tn*8]));
        bsv = exp(distB(rng[tn*8]))/asv;
      }
      else if(vg==5){
        asv = exp(distA(rng[tn*8]));
        csv = distC(rng[tn*8])*asv;
        bsv = exp(distB(rng[tn*8]))/asv;
      }
      else if(vg==6){
        csv = distC(rng[tn*8]);
        asv = a_lower+exp(distA(rng[tn*8]));
        bsv = b_lower+exp(distB(rng[tn*8]));
      }
      
      if(vk==1){
        double t1=(csv-tsv)/asv-lag;
        //double vt=(tsv-vli)/asv;
        for(int j=0;j<t_n;j++) { // 7 should be 0
          int t=(double) j;
            if(tsv > csv) cv = vli;
            else{
              if(t <= t1+tstoch) cv = tsv+asv*(lag+t-tstoch);
              //else if(j <= t1+tstoch-vt-lag+lsv) cv = csv;
              else cv = csv-bsv*(t-tstoch-t1);
            }
            
          if(cv<vli) 
            cv=vli; //+1;
          else if(cv>vl_max)
            cv=vl_max;				
          ft[tn][(int) (2*vl_max-2*cv)][j]++;
        }
      }
    }
  }
  
  /* ensure total simulated pop accounts for any rounding after division by ncores */
  double test_pop_d=(double) (n_s*ncores);
  
  NumericMatrix freq(t_n,ct_n); /* final freq matrix in Rcpp matrix format */
  
  /* sum thread specific frequencies and divide by total pop */
  for(int j=0;j<t_n;j++)
    for(int i=0;i<ct_n;i++) {
      double s=0;
      for(int tn=0;tn<ncores;tn++)
        s+=ft[tn][i][j];
      freq(j,i)=s/test_pop_d;
    }
    
  return freq;
  
}


