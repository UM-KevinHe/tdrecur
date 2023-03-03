//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' created on 2/14/2022: this one use stratified model

// [[Rcpp::export]]
List compute_facility_idx(arma::colvec &facility, int num_facility){
            std::vector<std::vector<int> > facility_idx(num_facility);//C++ list
            for(int i=0; i<facility.size(); i++){
                facility_idx[facility[i]].push_back(i);
            }
            return wrap(facility_idx);
}

//' this function calculates Sm (1 by p), sum of z in event set; O, the vector of total events in each 
//' facility; num_events, number of all events; events_per_day_facility, number of events each day in each facility
//'@param hosp_begin matrix of hospitalization begin date
//'@param max_d the largest number of hospitalization in each segment
//'@param D number of distinct failure time number of days)
//'@param z time-independent variables
// [[Rcpp::export]]
List compute_d12(arma::mat &hosp_begin, arma::colvec &max_d, int D,
arma::mat &z, arma::cube &Z_tv, arma::cube &t_tv, arma::mat &max_v,
	arma::colvec &facility, int num_facility){
	  int N = hosp_begin.n_rows;
	  int p1 = z.n_cols;	
    int p2 = Z_tv.n_slices; 
    int p = p1 + p2;
    arma::rowvec Sm(p); //Sm vector
    arma::colvec O(num_facility);
    int num_events = 0;
    arma::mat events_per_day_facility(num_facility, D); 
    arma::colvec events_per_day(D);
    arma::mat tempz, tempt;
    // i for patients, s for parameters
    for(int i = 0; i < N; i ++){
    	if(max_d(i) > 0){
      	O(facility(i)) += max_d(i);
      	num_events += max_d(i);
      	for (int s = 0; s < p1; s ++){
      		  Sm(s) += z(i, s) * max_d(i);
	  	  }
    		
    	  for (int j = 0; j < max_d(i); j ++ ){
    	  	  events_per_day_facility(facility(i), hosp_begin(i, j) - 1) += 1;
    	  	  events_per_day(hosp_begin(i, j) - 1) += 1;
            for (int s = 0; s < p2; s++) {
              tempz = Z_tv.slice(s);
              tempt = t_tv.slice(s);
              
              for (int t = max_v(i, s) - 1; t >= 0; t--) {
                if (tempt(i, t) <= hosp_begin(i, j)) {
                  Sm(s + p1) += tempz(i, t + 1);
                  break;
                }
              }
              if (tempt(i, 0) > hosp_begin(i, j)) {
                Sm(s + p1) += tempz(i, 0);
              }
            }
    	  }
    	}
    }
    return List::create(Named("Sm") = Sm,
    	Named("events_per_day_facility") = events_per_day_facility,
    	Named("events_per_day") = events_per_day,
    	Named("O") = O,
    	Named("num_events") = num_events);
}

//'@param l log partial likelihood
//'@param S0 D+1 vector, we will do day 0 (begining of day 1) to day D (end of day D)
//'@param t_start a vector of segment start date, starting from 1, minus 1 when using it.
//'@param t_end a vector of segment end date, starting from 1.
//'@param Sm 1 by p vector
//'@param ym number of changes on each day
// [[Rcpp::export]]
void ddloglik_cpp6(double l, double lchange, arma::colvec &L1, arma::mat &L2,
arma::colvec &t_start, arma::colvec &t_end, arma::rowvec &Sm, 
arma::mat events_per_day_facility, 
arma::mat& z, arma::cube &Z_tv, arma::cube &t_tv,
                          arma::mat &max_v, 
arma::colvec &beta, int N, List facility_idx,
                   arma::mat &hosp_begin, arma::colvec &max_d){
  L1.fill(0);
  L2.fill(0);
  lchange = l; // previous l
	l = 0;
	//int N = t_start.n_rows;
	int D = events_per_day_facility.n_cols; //number of distinct failure time (number of days)
	int p1 = z.n_cols;
	int p2 = Z_tv.n_slices;
  int p = p1 + p2;
  
  arma::mat zt = arma::zeros<arma::mat>(N, p2);
  arma::mat tempz, tempt;
  	
  double S0 = 0.0;
  arma::rowvec S1(p);
  arma::mat S2(p, p);
  double exp_z_beta0;
  arma::colvec zvi_beta =
    z * beta.rows(0, p1 - 1); // zbeta for variables that are time-independent
  arma::rowvec zi(p);
  arma::mat zti_beta(1, 1); // zbeta for time t, patient i

	// i for patients, j, s for covariates, k for time
  for (int k = 0; k < D; k++) { // for time
  	//cout<<"k = "<<k<<endl;

    zt.fill(0);
    for (int s = 0; s < p2; s++) { // for covariates
      tempz = Z_tv.slice(s); // now it's array, we can get one row
      tempt = t_tv.slice(s);
      for (int i = 0; i < N; i++) { // for patient
        // is it possible to combine these two for-loops and make it faster?
        // some steps below are not needed if all change times are between
        // t_start and t_end;
        // todo ???
        if ((k >= (t_start(i) - 1)) & (k < t_end(i))) {
          // only process when patient i is at risk on day k
          zt(i, s) = tempz(i, 0); // initialize the vector as 0s
          for (int flag = 0; flag < max_v(i, s); flag++) {
            if (tempt(i, flag) - 1 <= k) {
              // update if change happens before k
              zt(i, s) = tempz(i, flag + 1);
            }
          }
        }
      }
    }
      
  	for (int f=0; f<facility_idx.size(); f++)
  	{
      SEXP f_idx_temp = facility_idx[f];
      IntegerVector f_idx(f_idx_temp); //conversion
      int nj = f_idx.size();
  	  S0 = 0.0;
      S1.fill(0);
      S2.fill(0);
      
  		for(int i = 0; i < nj; i ++ ) {
    		if ((k >= (t_start(f_idx[i]) - 1)) & (k < t_end(f_idx[i]))) {
          zti_beta = zt.row(f_idx[i]) * beta.rows(p1, p - 1);
          exp_z_beta0 = exp(zvi_beta(f_idx[i]) + zti_beta(0, 0));
          S0 += exp_z_beta0;
          zi = join_rows(z.row(f_idx[i]), zt.row(f_idx[i]));
          S1 += zi * exp_z_beta0;
          S2 += zi.t() * zi * exp_z_beta0;
        }
  		 
    	}

      if(S0 == 0) // in case no one is at risk
      			{
      				S0 = 1e-6;
      				//cout<<"events_per_day_facility("<<f<<", "<<k<<") = "<<events_per_day_facility(f, k)<<endl;
      			}

      L2 += events_per_day_facility(f, k) / S0 * (S2 - S1.t() * S1 / S0);
      for (int s = 0; s < p; s++) {
        L1[s] -=
          events_per_day_facility(f, k) / S0 * S1(s);
      }      
      l -= events_per_day_facility(f, k) * log(S0);
  	}
  	
  }
  arma::mat Zbeta = Sm * beta;
  l += Zbeta(0, 0);
  lchange -= l;
	for(int j = 0; j < p1 + p2; j ++ ){ 
    	  L1[j] += Sm(j);
  }   
}

//'@param alpha facility effect for each person, length-N vector
//'@param alpha_star facility effect for each facility, length-num_facility vector
//'@param num_events total number of events
//'@param O number of events by facility
// [[Rcpp::export]]
void update_beta(double l, double lchange, arma::colvec &L1, arma::mat &L2,
arma::rowvec &Sm, 
arma::colvec &t_start, arma::colvec &t_end,
arma::mat events_per_day_facility, 
arma::mat& z, arma::cube &Z_tv, arma::cube &t_tv,
                   arma::mat &max_v,arma::colvec &beta,
arma::colvec &update, List facility_idx,
	arma::mat &hosp_begin, arma::colvec &max_d){
		int p1 = z.n_cols;
  int p2 = Z_tv.n_slices;
  int N = t_start.n_rows;
    
        // first step: update beta
        ddloglik_cpp6(l, lchange, L1, L2, t_start, t_end, Sm, events_per_day_facility, 
        z, Z_tv, t_tv,max_v, beta, N, facility_idx,
                  hosp_begin, max_d);
        
        update = arma::solve(L2 + arma::eye(p1 + p2, p1 + p2) * 1e-6, L1);
        
        beta += update;
        
        for(int i = 0; i < p1 + p2; i++){
        	if(beta(i) > 5){
        		update(i) -= (beta(i) - 5);
        		 beta(i) = 5;
        		 }
        		else 
        			{
        				if(beta(i) < -5){
        					update(i) -= (beta(i) + 5);
        					beta(i) = -5;
        					}
        			}
        }
	    	/*for(int i = 0; i< N; i++){
    	  	exp_z_beta0(i) = exp(z_beta(i) + cnh_beta(i)); //exp_z_beta0 has no facility information, length N;
    	  }*/
    	 
        //this loop is one iteration for beta
        
        
        
 
}

//'Estimate beta, lambda, and alpha
//'
//'@param max_d the largest number of hospitalization in each segment
//'@param delta indicator vector of event
//'@param z design matrix
//'@param tol tolerance of precision, when highest change in beta is lower than or equal to tol, iteration stops
//'@param facility 1 to num_facility
//'
//'@return estimated beta, lambda, alpha
//'
//'
//'@export
// [[Rcpp::export(rng = false)]]
List cox_breslow_dep(arma::mat &hosp_begin, arma::colvec &max_d, arma::colvec &t_start, arma::colvec &t_end, 
arma::mat &z, arma::cube &Z_tv, arma::cube &t_tv,
arma::mat &max_v,  double tolb, double tola, int maxiter, arma::colvec &facility, 
int num_facility, int D, arma::colvec beta){
    cout<<"this is shr14"<<endl;
    int p1 = z.n_cols;
    int p2 = Z_tv.n_slices;
    int p = p1 + p2;
    int N = hosp_begin.n_rows;
    //arma::colvec beta = arma::zeros<arma::colvec>(p); //beta not change back to 0 each iteration;
    	
    List result = compute_d12(hosp_begin, max_d, D, z, Z_tv, t_tv, max_v,
    	 facility, num_facility);
    List facility_idx = compute_facility_idx(facility, num_facility);

    arma::rowvec Sm = result("Sm");
    arma::mat events_per_day_facility = result("events_per_day_facility");
    arma::colvec events_per_day = result("events_per_day");
    arma::colvec O = result("O");
    int num_events = result("num_events");
    

    arma::colvec update = arma::ones<arma::colvec>(p);
    arma::colvec update1 = arma::ones<arma::colvec>(p);
    arma::colvec update_rel = arma::ones<arma::colvec>(p);
    
    double l;
    double l1;
    double lchange = 1.0;
    double lchange_rel = 1.0;
    
    arma::colvec L1(p);
    arma::mat L2(p, p);
    arma::colvec S0(D + 1);
    
    arma::colvec alpha(N, arma::fill::zeros);
    arma::colvec diffalpha(num_facility, arma::fill::value(10.0));
    arma::colvec diffalpha_rel(num_facility, arma::fill::value(10.0));
    arma::colvec alpha_star(num_facility);
    arma::colvec alpha_ori(num_facility);
    
    arma::colvec lambda(D, arma::fill::value(1e-3));
    arma::mat Epre(num_facility, D + 1);
    arma::mat exp_z_beta0(N, D);
    
    int iter = 0;
      
    while ((abs(update_rel).max() > tolb) & (abs(lchange_rel) > tolb) & (iter < maxiter)) {

        update_beta(l, lchange, L1, L2, 
        Sm, t_start, t_end, events_per_day_facility, 
        z, Z_tv, t_tv, max_v, beta,
        update, facility_idx,
                  hosp_begin, max_d
        );   	
        iter += 1;
        if(iter == 1){
        	update1 = update;
        	l1 = l;
        }
        update_rel = update / update1;
        if(iter > 1){
        	lchange_rel = lchange / (l - l1);
        }
    }
    if(abs(update_rel).max() > tolb){
    	cout<<"Warning: beta not converged!"<<endl;
    }
    cout<<"number of iterations for beta: "<<iter<<endl;
    
    //arma::mat Ui = ddloglik_Ui(L1, L2, S0, t_start, t_end, Sm, events_per_day_facility, 
    //    exp_z_beta0, z, cnh, beta, facility_idx, hosp_begin, max_d);
    //arma::colvec varbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));
    //arma::colvec varbeta = diagvec(pinv((L2) + arma::eye(p, p) * 1e-6));
    //arma::colvec sanvarbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6) * (Ui * Ui.t()) * pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));
    

    
    
    return List::create(Named("beta") = wrap(beta),
                        Named("update") = wrap(update), 
                        //Named("sebeta") = wrap(sqrt(varbeta)),
                        //Named("sansebeta") = wrap(sqrt(sanvarbeta)),
                        //Named("exp_z_beta0") = wrap(exp_z_beta0),
                        //Named("z_value") = wrap(z_value),
                        Named("L1") = L1,
                        Named("Info") = L2);//,
                        //Named("lambda") = wrap(lambda),
                        //Named("alpha_star") = wrap(alpha_star),
                        //Named("Epre") = wrap(Epre));
}
