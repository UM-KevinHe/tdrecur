//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' created on 12/11/2021: this one use stratified model

// [[Rcpp::export]]
List compute_facility_idx_ind(IntegerVector &facility, int num_facility){
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
List compute_d(NumericMatrix &hosp_begin, IntegerVector &max_d, int D,
NumericMatrix &z,
IntegerVector &facility, int num_facility){
	  int N = hosp_begin.nrow();
	  int p = z.ncol();
    NumericVector Sm(p); //Sm vector
    IntegerVector O(num_facility);
    int num_events = 0;
    IntegerMatrix events_per_day_facility(num_facility, D);
    IntegerVector dm(D);

    // i for patients, s for parameters
    for(int i = 0; i < N; i ++){
    	if(max_d(i) > 0){
    	  for (int j = 0; j < max_d(i); j ++ ){
    	  	  events_per_day_facility(facility(i), hosp_begin(i, j) - 1) += 1;
    	  	  dm(hosp_begin(i, j) - 1) += 1;
          	O(facility(i)) += 1;
          	num_events += 1;
          	for (int s = 0; s < p; s ++){
          		  Sm(s) += z(i, s);
    	  	  }
    	  }
    	}
    }
    return List::create(Named("Sm") = Sm,
    	Named("events_per_day_facility") = events_per_day_facility,
    	Named("dm") = dm,
    	Named("O") = O,
    	Named("num_events") = num_events);
}

//'@param l log partial likelihood
//'@param S0 D+1 vector, we will do day 0 (begining of day 1) to day D (end of day D)
//'@param t_start a vector of segment start date, starting from 1, minus 1 when using it.
//'@param t_end a vector of segment end date, starting from 1.
//'@param Sm p by 1 vector
//'@param ym number of changes on each day
//'@param exp_z_beta vector of exp(zi * beta + alphai), length N, with repeated KECC_ID
// [[Rcpp::export]]
void ddloglik_cpp(NumericVector &L1, NumericMatrix &L2,
NumericVector &S0, IntegerVector &t_start, IntegerVector &t_end, NumericVector &Sm,
IntegerMatrix events_per_day_facility,
NumericVector& exp_z_beta, NumericMatrix& z,
arma::colvec &beta, List facility_idx){
	std::fill(L1.begin(), L1.end(), 0);
	std::fill(L2.begin(), L2.end(), 0);
	//l = 0;
	//int N = t_start.length();
	int D = S0.length() - 1; //number of distinct failure time (number of days)
	int p = z.ncol();
  NumericMatrix S1(D + 1, p);
    NumericMatrix V(D + 1, p);
	// i for patients, j for covariates
	for (int f=0; f<facility_idx.size(); f++)
	{
    SEXP f_idx_temp = facility_idx[f];
    IntegerVector f_idx(f_idx_temp); //conversion
    int nj = f_idx.size();
	  std::fill(S0.begin(), S0.end(), 0);
	  std::fill(S1.begin(), S1.end(), 0);

		for(int i = 0; i < nj; i ++ ) {
		 	  S0[t_start[f_idx[i]] - 1] += exp_z_beta(f_idx[i]);
		    S0[t_end[f_idx[i]]] -= exp_z_beta(f_idx[i]);
    	  for(int j = 0; j < p; j ++){
    	   	  S1(t_start[f_idx[i]] - 1, j) += z(f_idx[i], j) * exp_z_beta(f_idx[i]);
    	   	  S1(t_end[f_idx[i]], j) -= z(f_idx[i], j) * exp_z_beta(f_idx[i]);
    	  }

  	}
  	// k for days
    for(int k = 1; k < D + 1; k ++){
    		S0(k) += S0(k - 1);
    		if (S0(k - 1) == 0)
    			{
    				S0(k - 1) = 1e-6;
    			}
    		for(int j = 0; j < p; j ++){
    		   S1(k, j) += S1(k - 1, j);
    		}
    }
    if(S0(D) == 0)
    			{
    				S0(D) = 1e-6;
    			}


    //i for time, j, k for parameters, s for patients

    //j = 0 to p: z;
    for(int j = 0; j < p; j ++){
	    std::fill(V.begin(), V.end(), 0);
    	for (int s=0; s< nj; s++){
        	  //k = 0 to p: z;
            for(int k=0; k< p; k++){
            	V(t_end[f_idx[s]], k) -= z(f_idx[s], k) * exp_z_beta(f_idx[s]) * z(f_idx[s], j);
            	V(t_start[f_idx[s]] - 1, k) += z(f_idx[s], k) * exp_z_beta(f_idx[s]) * z(f_idx[s], j);
            }

    	}
      for(int k = 0; k < p; k ++){
            for(int i=1; i< D + 1; i ++){
                V(i, k) += V(i - 1, k);
            }
      } //V is S2(, j) in this loop

        for(int k=0; k<p; k++){
            for(int i=0; i< D + 1;i++){
                V(i,k) = V(i,k)/S0[i]-S1(i,k)*S1(i,j)/(S0[i]*S0[i]);
            }
        } //V is I(, j) in this loop
        for(int i=0;i<D;i++){
            for(int k=0; k< p; k++){
                L2(k,j) += V(i,k) * events_per_day_facility(f, i);
            }

        }
    }


    //i for time, j for parameters
    for(int j = 0; j < p; j ++ ){
        for(int i = 0; i < D; i ++){
    	  	L1[j] -= events_per_day_facility(f, i) * S1(i, j) / S0(i); //need to match the id here
    	  }
    }
	}
	for(int j = 0; j < p; j ++ ){
    	  L1[j] += Sm(j);
  }
}

//'@param alpha facility effect for each person, length-N vector
//'@param alpha_star facility effect for each facility, length-num_facility vector
//'@param num_events total number of events
//'@param O number of events by facility
// [[Rcpp::export]]
void update_beta_ind(NumericVector &L1, NumericMatrix &L2,
NumericVector &S0, NumericVector &Sm,
IntegerVector &t_start, IntegerVector &t_end,
IntegerMatrix events_per_day_facility,
NumericVector& exp_z_beta0, NumericMatrix &z, arma::colvec &beta,
arma::colvec &update, arma::colvec &z_beta,
IntegerVector &facility, int num_facility,
int D, int N, List facility_idx){
		int p = z.ncol();
	    	for(int i = 0; i< N; i++){
    	  	exp_z_beta0(i) = exp(z_beta(i)); //exp_z_beta0 has no facility information, length N;
    	  }

        // first step: update beta
        ddloglik_cpp(L1, L2, S0, t_start, t_end, Sm, events_per_day_facility,
        exp_z_beta0, z, beta, facility_idx);

        update = arma::solve(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6, as<arma::colvec>(L1),
                             arma::solve_opts::allow_ugly);
        beta += update;
        /*
        for(int i = 0; i < p; i++){
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
        }*/
        z_beta = as<arma::mat>(z) * beta;
        //this loop is one iteration for beta




}

//'@param l log partial likelihood
//'@param S0 D+1 vector, we will do day 0 (begining of day 1) to day D (end of day D)
//'@param t_start a vector of segment start date, starting from 1, minus 1 when using it.
//'@param t_end a vector of segment end date, starting from 1.
//'@param Sm p by 1 vector
//'@param ym number of changes on each day
//'@param exp_z_beta vector of exp(zi * beta + alphai), length N, with repeated KECC_ID
// [[Rcpp::export]]
arma::mat ddloglik_Ui(NumericVector &L1, NumericMatrix &L2,
NumericVector &S0, IntegerVector &t_start, IntegerVector &t_end, NumericVector &Sm,
IntegerMatrix events_per_day_facility,
NumericVector& exp_z_beta, NumericMatrix& z,
arma::colvec &beta, List facility_idx, NumericMatrix &hosp_begin, IntegerVector &max_d){
	std::fill(L1.begin(), L1.end(), 0);
	std::fill(L2.begin(), L2.end(), 0);
	//l = 0;
	//int N = t_start.length();
	int D = S0.length() - 1; //number of distinct failure time (number of days)
	int p = z.ncol();
  NumericMatrix S1(D + 1, p);
  NumericMatrix V(D + 1, p);
  int N = t_start.length();
  arma::mat Ui(p, N, arma::fill::zeros);
  int d = 0; // indicator of event for patient i
	// i for patients, j for covariates
	for (int f=0; f<facility_idx.size(); f++)
	{
    SEXP f_idx_temp = facility_idx[f];
    IntegerVector f_idx(f_idx_temp); //conversion
    int nj = f_idx.size();
	  std::fill(S0.begin(), S0.end(), 0);
	  std::fill(S1.begin(), S1.end(), 0);
	  std::fill(V.begin(), V.end(), 0);

		for(int i = 0; i < nj; i ++ ) {
		 	  S0[t_start[f_idx[i]] - 1] += exp_z_beta(f_idx[i]);
		    S0[t_end[f_idx[i]]] -= exp_z_beta(f_idx[i]);
    	  for(int j = 0; j < p; j ++){
    	   	  S1(t_start[f_idx[i]] - 1, j) += z(f_idx[i], j) * exp_z_beta(f_idx[i]);
    	   	  S1(t_end[f_idx[i]], j) -= z(f_idx[i], j) * exp_z_beta(f_idx[i]);
    	  }

  	}
  	// k for days
    for(int k = 1; k < D + 1; k ++){
    		S0(k) += S0(k - 1);
    		if (S0(k - 1) == 0)
    			{
    				S0(k - 1) = 1e-6;
    			}
    		for(int j = 0; j < p; j ++){
    		   S1(k, j) += S1(k - 1, j);
    		}
    }
    if(S0(D) == 0)
    			{
    				S0(D) = 1e-6;
    			}



    //i for patient; t for hospital flag
    for(int i = 0; i < nj; i ++){
    	int t = 0;
	  	  for (int j = t_start(f_idx[i]) - 1; j < t_end(f_idx[i]); j ++){
	  	  	  if ((t < max_d(f_idx[i])) & (hosp_begin(f_idx[i], t) - 1 == j))
          			{
          				d = 1;
          				t ++ ;
          			}
          	for (int s = 0; s < p; s ++){
          		  Ui(s, f_idx[i]) += (z(f_idx[i], s) - S1(j, s) / S0(j)) * (d - events_per_day_facility(f, j) / S0(j) * exp_z_beta(f_idx[i]));
    	  	  }
    	  	  d = 0;
    	  }
    }
    //i for time, j, k for parameters, s for patients

    //j = 0 to p: z;
    for(int j = 0; j < p; j ++){
	    std::fill(V.begin(), V.end(), 0);
    	for (int s=0; s< nj; s++){
        	  //k = 0 to p: z;
            for(int k=0; k< p; k++){
            	V(t_end[f_idx[s]], k) -= z(f_idx[s], k) * exp_z_beta(f_idx[s]) * z(f_idx[s], j);
            	V(t_start[f_idx[s]] - 1, k) += z(f_idx[s], k) * exp_z_beta(f_idx[s]) * z(f_idx[s], j);
            }

    	}
      for(int k = 0; k < p; k ++){
            for(int i=1; i< D + 1; i ++){
                V(i, k) += V(i - 1, k);
            }
      } //V is S2(, j) in this loop

        for(int k=0; k<p; k++){
            for(int i=0; i< D + 1;i++){
                V(i,k) = V(i,k)/S0[i]-S1(i,k)*S1(i,j)/(S0[i]*S0[i]);
            }
        } //V is I(, j) in this loop
        for(int i=0;i<D;i++){
            for(int k=0; k< p; k++){
                L2(k,j) += V(i,k) * events_per_day_facility(f, i);
            }

        }
    }


    //i for time, j for parameters
    for(int j = 0; j < p; j ++ ){
        for(int i = 0; i < D; i ++){
    	  	L1[j] -= events_per_day_facility(f, i) * S1(i, j) / S0(i); //need to match the id here
    	  }
    }
	}
	for(int j = 0; j < p; j ++ ){
    	  L1[j] += Sm(j);
  }
  return Ui;
}
//'@param l log partial likelihood
//'@param S0 D+1 vector, we will do day 0 (begining of day 1) to day D (end of day D)
//'@param t_start a vector of segment start date, starting from 1, minus 1 when using it.
//'@param t_end a vector of segment end date, starting from 1.
//'@param Sm p by 1 vector
//'@param ym number of changes on each day
//'@param exp_z_beta vector of exp(zi * beta + alphai), length N, with repeated KECC_ID
// [[Rcpp::export]]
void ddloglik_cpp11(
NumericVector &S0, IntegerVector &t_start, IntegerVector &t_end,
NumericVector& exp_z_beta){
	std::fill(S0.begin(), S0.end(), 0);
	//l = 0;
	int N = t_start.length();
	int D = S0.length() - 1; //number of distinct failure time (number of days)
	// i for patients, j for covariates
	for(int i = 0; i < N; i ++ ) {
		 	  S0[t_start[i] - 1] += exp_z_beta(i);
		    S0[t_end[i]] -= exp_z_beta(i);
	}
	// k for days
    for(int k = 1; k < D + 1; k ++){
    		S0(k) += S0(k - 1);
    		if (S0(k - 1) == 0)
    			{
    				S0(k - 1) = 1e-6;
    			}
    }
    if(S0(D) == 0)
    			{
    				S0(D) = 1e-6;
    			}
}

//'@param alpha facility effect for each person, length-N vector
//'@param alpha_star facility effect for each facility, length-num_facility vector
//'@param num_events total number of events
//'@param O number of events by facility
// [[Rcpp::export]]
void update_alpha11(NumericVector &alpha,
NumericVector &S0,
IntegerVector &t_start, IntegerVector &t_end,
IntegerVector dm,
NumericVector& exp_z_beta, NumericVector& exp_z_beta0,
NumericVector &diffalpha, NumericVector &alpha_star,
IntegerVector &facility, int num_facility,
int D, NumericVector &lambda, int num_events, NumericVector &O, int N,
NumericMatrix &Epre){
	    	for(int i = 0; i< N; i++){
    	  	alpha(i) = alpha_star(facility(i));
    	  	exp_z_beta(i) = exp_z_beta0(i) * exp(alpha(i));
    	  }

        // first step: update S0
        ddloglik_cpp11(S0, t_start, t_end,
        exp_z_beta);

        // second step: update lambda
        for (int i = 0; i < D; i ++ ){
            lambda(i) = dm(i) / S0(i);
        }

        // third step: calculate alpha
        double Cnum = 0;
	      std::fill(Epre.begin(), Epre.end(), 0);
        // i for person, j for time
        for (int i = 0; i < N; i++){


          	if(t_end(i) < D){
          	  for(int j = t_start(i) - 1; j < t_end(i); j++){
                Epre(facility(i), j) += lambda(j) * exp_z_beta0(i);
              }
            }else{for(int j = t_start(i) - 1; j < D; j++){
                Epre(facility(i), j) += lambda(j) * exp_z_beta0(i);
              }
            }

        }
        NumericVector E = Epre(_, 0);
        for (int f=0; f<num_facility; f++)
        {for (int j = 1; j < D; j ++){
        	E(f) += Epre(f, j);
        }
        		Cnum += E(f);
        }
        //11/26/2020: this step can be faster if I have a list of sum(z) for each facility for each day, but that requires larger memory:n_f * D * p;


        double C = Cnum/num_events;
        for (int f = 0; f < num_facility; f ++){
        	if (O(f) == 0){
        		diffalpha(f) = alpha_star(f) + 10;
        		alpha_star(f) = -10;
        	}else{
        		diffalpha(f) = alpha_star(f) - log(C * O(f) / E(f));
            alpha_star(f) = log(C * O(f) / E(f));
          }
        }

}

//'Estimate beta, lambda, and alpha
//'
//'@param hosp_begin matrix of event time for each participant
//'@param max_d number of events for each row (participant)
//'@param delta indicator vector of event
//'@param z design matrix
//'@param tol tolerance of precision, when highest change in beta is lower than or equal to tol, iteration stops
//'@param facility 1 to num_facility
//'@param beta initial value for beta
//'
//'@return estimated beta, lambda, alpha
//'
//'
//'@export
// [[Rcpp::export(rng = false)]]
List cox_breslow_ind(NumericMatrix &hosp_begin, IntegerVector &max_d, IntegerVector &t_start, IntegerVector &t_end,
NumericMatrix &z, double tolb, double tola, int maxiter, IntegerVector &facility,
int num_facility, int D, arma::colvec beta){
    int p = z.ncol();
    int N = hosp_begin.nrow();
    //arma::colvec beta = arma::zeros<arma::colvec>(p); //beta not change back to 0 each iteration;

    List result = compute_d(hosp_begin, max_d, D, z, facility, num_facility);
    List facility_idx = compute_facility_idx_ind(facility, num_facility);

    NumericVector Sm = result("Sm");
    IntegerMatrix events_per_day_facility = result("events_per_day_facility");
    IntegerVector dm = result("dm");
    NumericVector O = result("O");
    int num_events = result("num_events");

    arma::colvec z_beta = as<arma::mat>(z) * beta;

    arma::colvec update = arma::ones<arma::colvec>(p);


    NumericVector L1(p);
    NumericMatrix L2(p, p);
    NumericVector S0(D + 1);

    NumericVector alpha(N);
    NumericVector diffalpha(num_facility, 10.0);
    NumericVector alpha_star(num_facility);
    NumericVector alpha_ori(num_facility);

    NumericVector exp_z_beta(N);
    NumericVector exp_z_beta0(N, 1.0);
    NumericVector lambda(D);
    int iter = 0;
    NumericMatrix Epre(num_facility, D + 1);

    //Rcout<<"this is shr12"<<endl;
    while ((abs(update).max() > tolb) & (iter < maxiter)) {

        update_beta_ind(L1, L2,
        S0, Sm, t_start, t_end, events_per_day_facility,
        exp_z_beta0, z, beta,
        update, z_beta, facility, num_facility,
        D, N, facility_idx);
        iter += 1;
    }
    if(iter == maxiter){
    	Rcout<<"Maximum number of iterations reached for beta!"<<endl;
    }
    //Rcout<<"number of iterations for beta: "<<iter<<endl;

    arma::mat Ui = ddloglik_Ui(L1, L2, S0, t_start, t_end, Sm, events_per_day_facility,
        exp_z_beta0, z, beta, facility_idx, hosp_begin, max_d);
    //arma::colvec varbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));
    arma::colvec varbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));
    arma::colvec sanvarbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6) * (Ui * Ui.t()) * pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));

	    	for(int i = 0; i< N; i++){
    	  	exp_z_beta0(i) = exp(z_beta(i)); //exp_z_beta0 has no facility information, length N;
    	  }
    iter = 0;
    while ((max(abs(diffalpha)) > tola) & (iter < maxiter)) {

        update_alpha11(alpha,
        S0, t_start, t_end, dm,
        exp_z_beta, exp_z_beta0,
        diffalpha, alpha_star, facility, num_facility,
        D, lambda, num_events, O,  N, Epre);
        iter += 1;
    }
    if(iter == maxiter){
    	Rcout<<"Maximum number of iterations reached for alpha!"<<endl;
    }
    //Rcout<<"number of iterations for alpha: "<<iter<<endl;

    return List::create(Named("beta") = wrap(beta),
                        Named("update") = wrap(update),
                        Named("sebeta") = wrap(sqrt(varbeta)),
                        Named("sansebeta") = wrap(sqrt(sanvarbeta)),
                        Named("L1") = L1,
                        Named("Info") = L2,
                        Named("lambda") = wrap(lambda),
                        Named("alpha_star") = wrap(alpha_star),
                        Named("Epre") = wrap(Epre));
}
