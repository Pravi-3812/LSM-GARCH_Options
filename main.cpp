#include <boost/format.hpp>
#include <boost/date_time.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cassert>
#include <algorithm>
#include <chrono>
#include "PolyfitBoost.hpp"
 // Include the chrono library for time measurement


using namespace std;
using namespace std::chrono;

class AmericanOptionsLSMC {
private:
    string option_type;
    double S0, strike, T, r, div, sigma;
    int M, simulations;
    double time_unit, discount;
    vector<vector<double>> MCprice_matrix;

public:
    // Constructor
    AmericanOptionsLSMC(string option_type, double S0, double strike, double T, int M, double r, double div, double sigma, int simulations) {
        assert(option_type == "call" || option_type == "put");
        assert(S0 > 0 && strike > 0 && T > 0 && M > 0 && r >= 0 && div >= 0 && sigma > 0 && simulations > 0);

        this->option_type = option_type;
        this->S0 = S0;
        this->strike = strike;
        this->T = T;
        this->M = M;
        this->r = r;
        this->div = div;
        this->sigma = sigma;
        this->simulations = simulations;

        this->time_unit = T / double(M);
        this->discount = exp(-r * time_unit);

        // Initialize the MC price matrix
        MCprice_matrix.resize(M + 1, vector<double>(simulations, 0.0));

    }

    // Monte Carlo Price Matrix
    void generate_MCprice_matrix() {
        random_device rd;
        mt19937 gen(rd());
        normal_distribution<> d(0, 1);

        MCprice_matrix[0].assign(simulations, S0);
        //for(int i = 0; i<simulations ; i++)
        //cout<<MCprice_matrix[0][i] ;

        for (int t = 1; t <= M; ++t) {
            vector<double> brownian(simulations / 2);
            for (int i = 0; i < simulations / 2; ++i) {
                brownian[i] = d(gen);
            }
        /*    vector<double> brownians2(simulations / 2);
            for(int i = 0 ; i<simulations ; i++) {
                brownians2[i] = -1*brownians2[i] ;
            }
            vector<double> brownian(simulations) ;
            for(int i = 0 ; i<simulations/2 ; i++)
                brownian[i] = brownians1[i] ;
            for(int i = simulations/2 ; i<simulations ; i++)
                brownian[i] = brownian[i-sumulations/2] ;
           for(int i = 0; i<simulations/2 ; ++i){
                MCprice_matrix[t][i] = MCprice_matrix[t-1][i]*exp(r-sigma*sigma/2)*time_unit + sigma*brownian[i]*sqrt(time_unit) ;
                MCprice_matrix[t][i+simulations/2] = MCprice_matrix[t-1][i]*exp(r-sigma*sigma/2)*time_unit + sigma*brownian[i]*sqrt(time_unit) ;
            } */


            // Brownian motion
            for (int i = 0; i < simulations / 2; ++i) {
                MCprice_matrix[t][i] = MCprice_matrix[t - 1][i] * exp((r - sigma * sigma / 2) * time_unit + sigma * brownian[i] * sqrt(time_unit));
                MCprice_matrix[t][simulations / 2 + i] = MCprice_matrix[t - 1][simulations / 2 + i] * exp((r - sigma * sigma / 2) * time_unit - sigma * brownian[i] * sqrt(time_unit));
            }
        }
    }

    // Payoff Matrix
    vector<vector<double>> get_MCpayoff() {
        //for(int i = 0 ; i < M+1 ; i++)
            //for(int j = 0 ; j<simulations; j++)
              //cout<<MCprice_matrix[i][j] ;

        vector<vector<double>> payoff(M + 1, vector<double>(simulations, 0.0));

        if (option_type == "call") {
            for (int t = 0; t <= M; ++t) {
                for (int i = 0; i < simulations; ++i) {
                    payoff[t][i] = max(MCprice_matrix[t][i] - strike, 0.0);
                }
            }
        } else if (option_type == "put") {
            for (int t = 0; t <= M; ++t) {
                for (int i = 0; i < simulations; ++i) {
                    payoff[t][i] = max(strike - MCprice_matrix[t][i], 0.0);
                }
            }
        }

        return payoff;
    }

    // Value Vector Calculation
  /*  vector<double> get_value_vector() {
        vector<vector<double>> payoff = get_MCpayoff();
        vector<vector<double>> value_matrix(M + 1, vector<double>(simulations, 0.0));
        value_matrix[M] = payoff[M];

        for (int t = M - 1; t > 0; --t) {
            vector<double> regression(simulations);
            vector<double> continuation_value(simulations);

            // Linear regression using simple least squares (polyfit equivalent in C++)
            for (int i = 0; i < simulations; ++i) {
                regression[i] = MCprice_matrix[t][i];
            }

            // Regression
            for (int i = 0; i < simulations; ++i) {
                continuation_value[i] = regression[i] * discount;
            }

            for (int i = 0; i < simulations; ++i) {
                value_matrix[t][i] = max(payoff[t][i], continuation_value[i]);
            }
        }

        vector<double> value_vector(simulations);
        for (int i = 0; i < simulations; ++i) {
            value_vector[i] = value_matrix[1][i] * discount;
        }

        return value_vector;
    }
*/

    vector<double> get_value_vector() {
        vector<vector<double>> payoff = get_MCpayoff();
        vector<vector<double>> value_matrix(M + 1, vector<double>(simulations, 0.0));
        value_matrix[M] = payoff[M];
        vector<double> value_vector = payoff[M] ;

        for (int t = M - 1; t > 0; --t) {
        vector<double> btm = value_matrix[t+1] ;
        for(int i = 0; i <simulations ; i++)
            btm[i] = discount*btm[i] ;
        vector<double> coeff = polyfit_boost(MCprice_matrix[t], btm, 5) ;
        vector<double> continuation_values = polyval(coeff , MCprice_matrix[t] ) ;
        for(int i = 0 ; i< simulations ; i++){
            if( payoff[t][i] > continuation_values[i] )
                value_matrix[t][i] = payoff[t][i] ;
            else
                value_matrix[t][i] = btm[i] ;
        }
        }
        for(int i = 0 ; i<simulations ; i++ )
            value_vector[i] = discount * value_matrix[1][i] ;

        return value_vector;
    }

    // Option Price
    double get_price() {
        vector<double> value_vector = get_value_vector();
        double sum = 0.0;

        for (int i = 0; i < simulations; ++i) {
            sum += value_vector[i];
        }

        return sum / simulations;
    }

    double delta()
    {
       double diff = S0 * 0.01 ;
       //AmericanOptionsLSMC americanOption("put", S0 , 40.0, 1, 50, 0.06, 0.06, 0.2, 1500);
                AmericanOptionsLSMC americanoption1(option_type, S0 + diff, strike, T, M, r, div, sigma, simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption1.generate_MCprice_matrix();

                // Get the option price
                double price1 = americanoption1.get_price();
                AmericanOptionsLSMC americanoption2(option_type, S0 - diff, strike, T, M, r, div, sigma, simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption2.generate_MCprice_matrix();

                // Get the option price
                double price2 = americanoption1.get_price();
                double delk = (price1 - price2) / (2*diff) ;
                return delk ;
    }
    double gamma()
    {
      double diff = S0 * 0.01 ;
       //AmericanOptionsLSMC americanOption("put", S0 , 40.0, 1, 50, 0.06, 0.06, 0.2, 1500);
                AmericanOptionsLSMC americanoption1(option_type, S0 + diff, strike, T, M, r, div, sigma, simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption1.generate_MCprice_matrix();

                // Get the option price
                double price1 = americanoption1.delta();
                AmericanOptionsLSMC americanoption2(option_type, S0 - diff, strike, T, M, r, div, sigma, simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption2.generate_MCprice_matrix();

                // Get the option price
                double price2 = americanoption1.delta();
                double delk = (price1 - price2) / (2*diff) ;
                return delk ;

    }

    double vega()
    {
               double diff = sigma * 0.01 ;
       //AmericanOptionsLSMC americanOption("put", S0 , 40.0, 1, 50, 0.06, 0.06, 0.2, 1500);
                AmericanOptionsLSMC americanoption1(option_type, S0 , strike, T, M, r, div, sigma + diff, simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption1.generate_MCprice_matrix();

                // Get the option price
                double price1 = americanoption1.get_price();
                AmericanOptionsLSMC americanoption2(option_type, S0 , strike, T, M, r, div, sigma - diff, simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption2.generate_MCprice_matrix();

                // Get the option price
                double price2 = americanoption1.get_price();
                double delk = (price1 - price2) / (2*diff) ;
                return delk ;
    }

    double rho()
    {
        double diff = r*0.01 ;
                AmericanOptionsLSMC americanoption1(option_type, S0 , strike, T, M, r + diff, div, sigma , simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption1.generate_MCprice_matrix();

                // Get the option price
                double price1 = americanoption1.get_price();
                AmericanOptionsLSMC americanoption2(option_type, S0 , strike, T, M, r - diff , div, sigma , simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption2.generate_MCprice_matrix();

                // Get the option price
                double price2 = americanoption1.get_price();
                double delk = (price1 - price2) / (2*diff) ;
                return delk ;
    }

    double theta()
    {
                double Timedelta = 1/252 ;
                AmericanOptionsLSMC americanoption1(option_type, S0 , strike, T + Timedelta, M, r, div, sigma , simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption1.generate_MCprice_matrix();

                // Get the option price
                double price1 = americanoption1.get_price();
                AmericanOptionsLSMC americanoption2(option_type, S0 , strike, T - Timedelta , M , r , div, sigma , simulations) ;
                // Generate the Monte Carlo price matrix
                americanoption2.generate_MCprice_matrix();

                // Get the option price
                double price2 = americanoption1.get_price();
                double delk = (price1 - price2) / (2*Timedelta) ;
                return delk ;
    }



};
double volatilityfutt(double errort , double volatility , double c0 = 0.4599 , double c1 = 0.378 , double c2 = 0.622)
{
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(0.5, 1);
    return d(gen) *(c0 + c1*errort*errort + c2*volatility) ;
}
double assetvalue(double Asset , double r , double volatility, double Asset_premium , double t_error )
{
 double volatility_fut = volatilityfutt(t_error , volatility ) ;
 double asset3 = Asset * exp(r - 0.5*volatility_fut +  t_error + Asset_premium*sqrt(volatility_fut)) ;
 return asset3 ;
}

int main() {
    // Record the start time

/*
    // Iterate over multiple values of S0, volatility, and time to maturity
    for (double S0 : {36.0, 38.0, 40.0, 42.0, 44.0}) {  // Initial stock price values
        for (double vol : {0.2, 0.4}) {  // Volatility values
            for (double T : {1.0, 2.0}) {  // Times to maturity
                // Create the American option
                AmericanOptionsLSMC americanOption("put", S0, 40.0, T, 50, 0.06, 0.06, vol, 4);

                // Generate the Monte Carlo price matrix
                americanOption.generate_MCprice_matrix();

                // Get the option price
                double price = americanOption.get_price();

                // Print the result
                cout << "Initial price: " << S0 << ", Sigma: " << vol << ", Expire: " << T
                     << " --> Option Value: " << price << endl;
            }
        }
    } */
    double volatility = volatilityfutt(0.7 , 0.06) , stock = assetvalue(36 , 0.06 , 0.06 , 0.06 , 0.7 ) ;
//AmericanOptionsLSMC(string option_type, double S0, double strike, double T, int M, double r, double div, double sigma, int simulations)
                AmericanOptionsLSMC americanOption("put", stock , 40.0, 1, 50, 0.06, 0.06, volatility , 1500);

                // Generate the Monte Carlo price matrix
                americanOption.generate_MCprice_matrix();

                // Get the option price
                double price = americanOption.get_price();
                cout<<price ;
    



    return 0;
}
