// Copyright (c) 2007 Intel Corp.
// Modernized by J. Daniel Garcia 2021
// Universidad Carlos III de Madrid
//
// Black-Scholes
// Analytical method for calculating European Options
//
//
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice
// Hall, John C. Hull,

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

constexpr int NUM_RUNS = 100;

template <typename fptype>
struct OptionData {
  fptype s;          // spot price
  fptype strike;     // strike price
  fptype r;          // risk-free interest rate
  fptype divq;       // dividend rate
  fptype v;          // volatility
  fptype t;          // time to maturity or option expiration in years
  //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
  char OptionType;   // Option type.  "P"=PUT, "C"=CALL
  fptype divs;       // dividend vals (not used in this test)
  fptype DGrefval;   // DerivaGem Reference Value
};

template <typename fptype>
OptionData<fptype> *data;

template<typename Number>
std::istream & operator>>(std::istream & is, OptionData<Number> & option)
{
  is >> option.s
     >> option.strike
     >> option.r
     >> option.divq
     >> option.v
     >> option.t
     >> option.OptionType
     >> option.divs
     >> option.DGrefval;
}

template <typename fptype> fptype *prices;
int numOptions;

int    * otype;
template <typename fptype> fptype * sptprice;
template <typename fptype> fptype * strike;
template <typename fptype> fptype * rate;
template <typename fptype> fptype * volatility;
template <typename fptype> fptype * otime;
int numError = 0;
int nThreads;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
template <typename T>
constexpr T inv_sqrt_2xPI = 0.39894228040143270286;

template<typename Number>
Number cdf_normal(Number x)
{
  // Check for negative value of x
  int sign;
  if (x < 0.0) {
    x = -x;
    sign = 1;
  }
  else {
    sign = 0;
  }

  Number xK2 = 1 / (1 + Number(0.2316419) * x);
  Number xK2_2 = xK2 * xK2;
  Number xK2_3 = xK2_2 * xK2;
  Number xK2_4 = xK2_3 * xK2;
  Number xK2_5 = xK2_4 * xK2;

  Number xLocal_1 = xK2 * 0.319381530;
  Number xLocal_2 = xK2_2 * (-0.356563782);
  Number xLocal_3 = xK2_3 * 1.781477937;
  xLocal_2 = xLocal_2 + xLocal_3;
  xLocal_3 = xK2_4 * (-1.821255978);
  xLocal_2 = xLocal_2 + xLocal_3;
  xLocal_3 = xK2_5 * 1.330274429;
  xLocal_2 = xLocal_2 + xLocal_3;
  xLocal_1 = xLocal_2 + xLocal_1;

  // Compute NPrimeX term common to both four & six decimal accuracy calcs
  Number xNPrimeofX = std::exp(-0.5f * x * x) * inv_sqrt_2xPI<Number>;
  Number xLocal = 1.0 - xLocal_1 * xNPrimeofX;

  if (sign) {
    return 1.0 - xLocal;
  }
  else {
    return xLocal;
  }
}

template <typename Number>
Number BlkSchlsEqEuroNoDiv(Number spot_price,
    Number strike, Number rate, Number volatility,
    Number time, int option_type)
{
  Number xPowerTerm = volatility * volatility / 2;
  Number d1 = (rate + xPowerTerm) * time + std::log(spot_price / strike);

  Number den = volatility * std::sqrt(time);
  d1 = d1 / den;
  Number d2 = d1 - den;

  Number NofXd1 = cdf_normal(d1);
  Number NofXd2 = cdf_normal(d2);

  Number FutureValueX = strike * std::exp(-(rate)*(time));
  if (option_type == 0) {
    return spot_price * NofXd1 - FutureValueX * NofXd2;
  } else {
    return FutureValueX * (1 - NofXd2) - spot_price * (1 - NofXd1);
  }
}


template <typename Number>
void compute_values() {
  for (int j=0; j<NUM_RUNS; j++) {
    for (int i=0; i<numOptions; i++) {
      /* Calling main function to calculate option value based on
       * Black & Scholes's equation.
       */
      Number price = BlkSchlsEqEuroNoDiv(sptprice<Number>[i], strike<Number>[i],
          rate<Number>[i], volatility<Number>[i], otime<Number>[i],
          otype[i]);
      prices<Number>[i] = price;
    }
  }
}

int main (int argc, char **argv)
{
  using fptype = float;

  std::cout << "PARSEC Benchmark Suite" << std::endl;

  if (argc != 4)
  {
    printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
    exit(1);
  }
  nThreads = std::stoi(argv[1]);
  std::string inputFile = argv[2];
  std::string outputFile = argv[3];

  //Read input data from file
  std::ifstream input{inputFile};
  if(!input) {
    std::cerr << "ERROR: Unable to open file `" << inputFile << "'.\n";
    exit(1);
  }
  input >> numOptions;
  if(!input) {
    std::cerr << "ERROR: Unable to read from file `" << inputFile << "'.\n";
    exit(1);
  }
  if(nThreads > numOptions) {
    std::cerr << "WARNING: Not enough work,"
      << "reducing number of threads to match number of options.\n";
    nThreads = numOptions;
  }

  if(nThreads != 1) {
    std::cerr << "Error: <nthreads> must be 1 (serial version)\n";
    exit(1);
  }

  // alloc spaces for the option data
  data<fptype> = (OptionData<fptype>*)malloc(numOptions*sizeof(OptionData<fptype>));
  prices<fptype> = (fptype*)malloc(numOptions*sizeof(fptype));
  for (int loopnum = 0; loopnum < numOptions; ++ loopnum )
  {
    input >> data<fptype>[loopnum];
  }
  if(!input) {
    std::cerr << "ERROR: Unable to read from file `" << inputFile << "'.\n";
    exit(1);
  }

  std::cout << "Num of Options: " << numOptions << "\n";
  std::cout << "Num of Runs: " << NUM_RUNS << "\n";

  constexpr int PAD = 256;
  constexpr int LINESIZE = 64;

  fptype * buffer = (fptype *) malloc(5 * numOptions * sizeof(fptype) + PAD);
  sptprice<fptype> = (fptype *) (((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
  strike<fptype> = sptprice<fptype> + numOptions;
  rate<fptype> = strike<fptype> + numOptions;
  volatility<fptype> = rate<fptype> + numOptions;
  otime<fptype> = volatility<fptype> + numOptions;

  int * buffer2 = (int *) malloc(numOptions * sizeof(fptype) + PAD);
  otype = (int *) (((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

  for (int i=0; i<numOptions; i++) {
    otype[i]      = (data<fptype>[i].OptionType == 'P') ? 1 : 0;
    sptprice<fptype>[i]   = data<fptype>[i].s;
    strike<fptype>[i]     = data<fptype>[i].strike;
    rate<fptype>[i]       = data<fptype>[i].r;
    volatility<fptype>[i] = data<fptype>[i].v;
    otime<fptype>[i]      = data<fptype>[i].t;
  }

  std::cout << "Size of data: "
    << numOptions * (sizeof(OptionData<fptype>) + sizeof(int))
    << "\n";

  //serial version
  int tid=0;
  compute_values<fptype>();

  //Write prices to output file
  std::ofstream output{outputFile};
  if(!output) {
    std::cerr << "ERROR: Unable to open file `" << outputFile << "'.\n";
    exit(1);
  }
  output << numOptions << "\n";
  if(!output) {
    std::cerr << "ERROR: Unable to write to file `" << outputFile << "'.\n";
    exit(1);
  }
  output.precision(18);
  output.setf(std::ios::fixed );
  for(int i=0; i<numOptions; i++) {
    output << prices<fptype>[i] << "\n";
  }
  if(!output) {
    std::cerr << "ERROR: Unable to write to file `" << outputFile << "'.\n";
    exit(1);
  }

  free(data<fptype>);
  free(prices<fptype>);

  return 0;
}

