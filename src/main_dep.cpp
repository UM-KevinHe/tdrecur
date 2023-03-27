//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' created on 3/7/2022: this file saves common functions used
#include "shr_dep.h"
#include "shr_omp.h"

// [[Rcpp::export]]
List compute_facility_idx(arma::colvec &facility, int num_facility){
            std::vector<std::vector<int> > facility_idx(num_facility);//C++ list
            for(unsigned int i=0; i<facility.size(); i++){
                facility_idx[facility[i]].push_back(i);
            }
            return wrap(facility_idx);
}


//'@param alpha facility effect for each person, length-N vector
//'@param alpha_star facility effect for each facility, length-num_facility vector
//'@param num_events total number of events
//'@param events_per_facility number of events by facility
// [[Rcpp::export]]
void update_beta(double l, double lchange, arma::colvec &L1, arma::mat &L2,
arma::rowvec &Sm,
arma::colvec &t_start, arma::colvec &t_end,
arma::mat events_per_day_facility,
arma::mat& z, arma::cube &Z_tv, arma::cube &t_tv,
                   arma::mat &max_v,arma::colvec &beta,
arma::colvec &update, arma::colvec &update_rel, List facility_idx,
	arma::mat &hosp_begin, arma::colvec &max_d, bool parallel = true, const unsigned int &nthreads=1){
		int p1 = z.n_cols;
  int p2 = Z_tv.n_slices;
  int N = t_start.n_rows;

        // first step: update beta
        if(parallel){
        	ddloglik_omp(L1, L2, t_start, t_end, Sm, events_per_day_facility,
        z, Z_tv, t_tv,max_v, beta, N, facility_idx,
                  hosp_begin, max_d, nthreads);
        }
        else{
        ddloglik_cpp6(l, lchange, L1, L2, t_start, t_end, Sm, events_per_day_facility,
        z, Z_tv, t_tv,max_v, beta, N, facility_idx,
                  hosp_begin, max_d);
                }

        update = arma::solve(L2 + arma::eye(p1 + p2, p1 + p2) * 1e-6, L1);
        update_rel = update / beta;
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


//'@param alpha facility effect for each person, length-N vector
//'@param alpha_star facility effect for each facility, length-num_facility vector
//'@param num_events total number of events
//'@param events_per_facility number of events by facility
// [[Rcpp::export]]
void update_alpha11(arma::colvec &alpha,
arma::colvec &S0,
arma::colvec &t_start, arma::colvec &t_end,
arma::colvec &events_per_day,
arma::mat& exp_z_beta0,
arma::colvec &diffalpha, arma::colvec &alpha_star,
arma::mat &z, 	arma::cube &Z_tv, arma::cube &t_tv, arma::mat &max_v,
   arma::colvec &beta,
arma::colvec &facility, int num_facility,
int D, arma::colvec &lambda, int num_events, arma::colvec &events_per_facility, int N,
arma::mat &Epre, arma::mat &hosp_begin, arma::colvec &max_d, bool parallel = true, const unsigned int &nthreads=1){


        // first step: update S0
        if(parallel){
        ddloglik_all_omp(alpha,
                   t_start, t_end,
                   Z_tv, t_tv,
                   max_v, z, beta,
                   D,
                   N, S0, exp_z_beta0,
                  hosp_begin, max_d, nthreads);
        }
        else{
        	ddloglik_cpp14(alpha,
                   t_start, t_end,
                   Z_tv, t_tv,
                   max_v, z, beta,
                   D,
                   N, S0, exp_z_beta0);
      }

        // second step: update lambda
        for (int i = 0; i < D; i ++ ){
            lambda(i) = events_per_day(i) / S0(i);
        }

        // third step: calculate alpha
        double Cnum = 0;
	      std::fill(Epre.begin(), Epre.end(), 0);
        // i for person, j for time
        for (int i = 0; i < N; i++){


          	if(t_end(i) < D){
          	  for(int j = t_start(i) - 1; j < t_end(i); j++){
                Epre(facility(i), j) += lambda(j) * exp_z_beta0(i, j);
              }
            }else{for(int j = t_start(i) - 1; j < D; j++){
                Epre(facility(i), j) += lambda(j) * exp_z_beta0(i, j);
              }
            }

        }
        arma::colvec E = Epre.col(0);
        for (int f=0; f<num_facility; f++)
        {for (int j = 1; j < D; j ++){
        	E(f) += Epre(f, j);
        }
        		Cnum += E(f);
        }
        //11/26/2020: this step can be faster if I have a list of sum(z) for each facility for each day, but that requires larger memory:n_f * D * p;


        double C = Cnum/num_events;
        for (int f = 0; f < num_facility; f ++){
        	if (events_per_facility(f) == 0){
        		diffalpha(f) = alpha_star(f) + 10;
        		alpha_star(f) = -10;
        	}else{
        		diffalpha(f) = alpha_star(f) - log(C * events_per_facility(f) / E(f));
            alpha_star(f) = log(C * events_per_facility(f) / E(f));
          }
        }

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
int num_facility, int D, arma::colvec beta,
bool parallel = true, const unsigned int &nthreads=1){

    int p1 = z.n_cols;
    int p2 = Z_tv.n_slices;
    int p = p1 + p2;
    int N = hosp_begin.n_rows;
    //arma::colvec beta = arma::zeros<arma::colvec>(p); //beta not change back to 0 each iteration;
    List result;
    if (parallel){
    	result = compute_d_omp(hosp_begin, max_d, D, z, Z_tv, t_tv, max_v,
    	 facility, num_facility, nthreads);
    }
    else{
    	result = compute_d12(hosp_begin, max_d, D, z, Z_tv, t_tv, max_v,
    	 facility, num_facility);
    }
    List facility_idx = compute_facility_idx(facility, num_facility);

    arma::rowvec Sm = result("Sm");
    arma::mat events_per_day_facility = result("events_per_day_facility");
    arma::colvec events_per_day = result("events_per_day");
    arma::colvec events_per_facility = result("events_per_facility");
    int num_events = result("num_events");


    arma::colvec update = arma::ones<arma::colvec>(p);
    arma::colvec update1 = arma::ones<arma::colvec>(p);
    arma::colvec update_rel = arma::ones<arma::colvec>(p);

    double l = 0.0;
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
        update, update_rel, facility_idx,
                  hosp_begin, max_d, parallel, nthreads
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
    	Rcout<<"Warning: beta not converged!"<<endl;
    }
    Rcout<<"number of iterations for beta: "<<iter<<endl;

    //arma::mat Ui = ddloglik_Ui(L1, L2, S0, t_start, t_end, Sm, events_per_day_facility,
    //    exp_z_beta0, z, cnh, beta, facility_idx, hosp_begin, max_d);
    //arma::colvec varbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));
    //arma::colvec varbeta = diagvec(pinv((L2) + arma::eye(p, p) * 1e-6));
    //arma::colvec sanvarbeta = diagvec(pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6) * (Ui * Ui.t()) * pinv(as<arma::mat>(L2) + arma::eye(p, p) * 1e-6));


    iter = 0;
    while ((max(abs(diffalpha)) > tola) & (iter < maxiter)) {

        update_alpha11(alpha,
        S0, t_start, t_end, events_per_day,
        exp_z_beta0,
        diffalpha, alpha_star, z, Z_tv, t_tv, max_v,
    beta,
   	facility, num_facility,
        D, lambda, num_events, events_per_facility,  N, Epre,
        hosp_begin, max_d, parallel, nthreads
        );
        iter += 1;
    }
    Rcout<<"number of iterations for alpha: "<<iter<<endl;

    return List::create(Named("beta") = wrap(beta),
                        Named("update") = wrap(update),
                        //Named("sebeta") = wrap(sqrt(varbeta)),
                        //Named("sansebeta") = wrap(sqrt(sanvarbeta)),
                        //Named("exp_z_beta0") = wrap(exp_z_beta0),
                        //Named("z_value") = wrap(z_value),
                        Named("L1") = L1,
                        Named("Info") = L2,
                        Named("lambda") = wrap(lambda),
                        Named("alpha_star") = wrap(alpha_star),
                        Named("Epre") = wrap(Epre));
}
