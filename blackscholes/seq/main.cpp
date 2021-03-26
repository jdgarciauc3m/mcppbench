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

template<typename fptype>
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

template<typename Number>
std::istream & operator>>(std::istream & is, OptionData<Number> & option)
{
  is >> option.s >> option.strike >> option.r >> option.divq >> option.v
     >> option.t >> option.OptionType >> option.divs >> option.DGrefval;
  return is;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244

template<typename Number>
Number cdf_normal(Number x)
{
  constexpr Number inv_sqrt_2xPI = 0.39894228040143270286;
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
  Number xNPrimeofX = std::exp(-0.5f * x * x) * inv_sqrt_2xPI;
  Number xLocal = 1.0 - xLocal_1 * xNPrimeofX;

  if (sign) {
    return 1.0 - xLocal;
  }
  else {
    return xLocal;
  }
}

template<typename Number>
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

  Number FutureValueX = strike * std::exp(-(rate) * (time));
  if (option_type == 0) {
    return spot_price * NofXd1 - FutureValueX * NofXd2;
  }
  else {
    return FutureValueX * (1 - NofXd2) - spot_price * (1 - NofXd1);
  }
}


template<typename Number>
void compute_values(
    int n_options,
    Number *p_sptprice,
    Number *p_strike,
    Number *p_rate,
    Number *p_volatility,
    Number *p_otime,
    int *p_otype,
    Number * p_prices
)
{
  for (int j = 0; j < NUM_RUNS; j++) {
    for (int i = 0; i < n_options; i++) {
      /* Calling main function to calculate option value based on
       * Black & Scholes's equation.
       */
      Number price = BlkSchlsEqEuroNoDiv(p_sptprice[i], p_strike[i],
          p_rate[i], p_volatility[i], p_otime[i],
          p_otype[i]);
      p_prices[i] = price;
    }
  }
}

int main(int argc, char **argv)
{
  using fptype = float;

  std::cout << "PARSEC Benchmark Suite" << std::endl;

  if (argc != 4) {
    printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
    exit(1);
  }
  int nThreads = std::stoi(argv[1]);
  std::string inputFile = argv[2];
  std::string outputFile = argv[3];

  //Read input data from file
  std::ifstream input{inputFile};
  if (!input) {
    std::cerr << "ERROR: Unable to open file `" << inputFile << "'.\n";
    exit(1);
  }
  int numOptions;
  input >> numOptions;
  if (!input) {
    std::cerr << "ERROR: Unable to read from file `" << inputFile << "'.\n";
    exit(1);
  }
  if (nThreads > numOptions) {
    std::cerr << "WARNING: Not enough work,"
              << "reducing number of threads to match number of options.\n";
    nThreads = numOptions;
  }

  if (nThreads != 1) {
    std::cerr << "Error: <nthreads> must be 1 (serial version)\n";
    exit(1);
  }

  // alloc spaces for the option data
  auto data = (OptionData<fptype> *) malloc(numOptions * sizeof(OptionData<fptype>));
  auto prices = (fptype *) malloc(numOptions * sizeof(fptype));
  for (int loopnum = 0; loopnum < numOptions; ++loopnum) {
    input >> data[loopnum];
  }
  if (!input) {
    std::cerr << "ERROR: Unable to read from file `" << inputFile << "'.\n";
    exit(1);
  }

  std::cout << "Num of Options: " << numOptions << "\n";
  std::cout << "Num of Runs: " << NUM_RUNS << "\n";

  constexpr int PAD = 256;
  constexpr int LINESIZE = 64;

  auto buffer = (fptype *) malloc(5 * numOptions * sizeof(fptype) + PAD);
  auto sptprice = (fptype *) (((unsigned long long) buffer + PAD) & ~(LINESIZE - 1));
  auto strike = sptprice + numOptions;
  auto rate = strike + numOptions;
  auto volatility = rate + numOptions;
  auto otime = volatility + numOptions;

  auto buffer2 = (int *) malloc(numOptions * sizeof(fptype) + PAD);
  auto otype = (int *) (((unsigned long long) buffer2 + PAD) & ~(LINESIZE - 1));

  for (int i = 0; i < numOptions; i++) {
    otype[i] = (data[i].OptionType == 'P') ? 1 : 0;
    sptprice[i] = data[i].s;
    strike[i] = data[i].strike;
    rate[i] = data[i].r;
    volatility[i] = data[i].v;
    otime[i] = data[i].t;
  }

  std::cout << "Size of data: "
            << numOptions * (sizeof(OptionData<fptype>) + sizeof(int))
            << "\n";

  //serial version
  compute_values(numOptions, sptprice, strike, rate,
      volatility, otime, otype, prices);
  //Write prices to output file
  std::ofstream output{outputFile};
  if (!output) {
    std::cerr << "ERROR: Unable to open file `" << outputFile << "'.\n";
    exit(1);
  }
  output << numOptions << "\n";
  if (!output) {
    std::cerr << "ERROR: Unable to write to file `" << outputFile << "'.\n";
    exit(1);
  }
  output.precision(18);
  output.setf(std::ios::fixed);
  for (int i = 0; i < numOptions; i++) {
    output << prices[i] << "\n";
  }
  if (!output) {
    std::cerr << "ERROR: Unable to write to file `" << outputFile << "'.\n";
    exit(1);
  }

  free(data);
  free(prices);

  return 0;
}

