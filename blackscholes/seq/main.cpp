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
#include <utility>
#include <vector>
#include <chrono>
#include <fmt/ostream.h>

constexpr int NUM_RUNS = 100;

class invalid_io {
public:
  explicit invalid_io(std::string file_name)
      : file_name_{std::move(file_name)} {}

  [[nodiscard]] std::string what() const  {
    return "ERROR: Unable to access file `" + file_name_ + "'.";
  }

private:
  std::string file_name_;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244

template<typename Number>
Number cdf_normal(Number x)
{
  constexpr double inv_sqrt_2xPI = 0.39894228040143270286;
  // Check for negative value of x
  int sign;
  if (x < 0.0) {
    x = -x;
    sign = 1;
  }
  else {
    sign = 0;
  }

  Number xK2 = 0.2316419 * x;
  xK2 = 1.0 + xK2;
  xK2 = 1.0 / xK2;
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
  Number expValues = std::exp(-0.5f * x * x);
  Number xNPrimeofX = expValues;
  xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

  Number xLocal = xLocal_1 * xNPrimeofX;
  xLocal = 1.0 - xLocal;

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

template<typename fptype>
struct option {
  fptype spot_price;
  fptype strike_price;
  fptype rate;          // risk-free interest rate
  fptype dividend_rate;
  fptype volatility;
  fptype time;          // time to maturity or option expiration in years
  //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
  char option_type;   // Option type.  "P"=PUT, "C"=CALL
  fptype divs;       // dividend vals (not used in this test)
  fptype DGrefval;   // DerivaGem Reference Value
};

template<typename Number>
std::istream & operator>>(std::istream & is, option<Number> & option)
{
  is >> option.spot_price >> option.strike_price >> option.rate
     >> option.dividend_rate >> option.volatility
     >> option.time >> option.option_type >> option.divs >> option.DGrefval;
  return is;
}

template<typename Name>
class input_option_portfolio;

template<typename Number>
std::istream & operator>>(std::istream & is, input_option_portfolio<Number> & pf);

template<typename Number>
class input_option_portfolio {
public:
  [[nodiscard]] std::size_t size() const { return options.size(); }

  [[nodiscard]] std::size_t data_size() const
  {
    return options.size() * sizeof(option<Number>);
  }

  const option<Number> & operator[](std::size_t i) const { return options[i]; }

  friend std::istream & operator>><>(std::istream & is, input_option_portfolio<Number> & pf);

private:
  std::vector<option<Number>> options;
};

template<typename Number>
std::istream & operator>>(std::istream & is, input_option_portfolio<Number> & pf)
{
  int num_options;  // NOLINT
  is >> num_options;
  if (!is) {
    return is;
  }

  pf.options.reserve(num_options);

  for (int i = 0; i < num_options; ++i) {
    option<Number> o; // NOLINT
    is >> o;
    //if (!is) { return is; }
    pf.options.push_back(o);
  }
  return is;
}

template<typename Number>
input_option_portfolio<Number> read_input_portfolio(const std::string & file_name)
{
  std::ifstream input{file_name};
  if (!input) { throw invalid_io{file_name}; }
  input_option_portfolio<Number> input_portfolio;
  input >> input_portfolio;
  if (!input) { throw invalid_io{file_name}; }
  return input_portfolio;
}


template<class Number>
class compute_option_portfolio {
public:
  explicit compute_option_portfolio(const input_option_portfolio<Number> & pf)
  {
    const std::size_t num_options = pf.size();
    reserve(num_options);
    for (std::size_t i = 0; i < num_options; ++i) {
      spot_price.push_back(pf[i].spot_price);
      strike_price.push_back(pf[i].strike_price);
      rate.push_back(pf[i].rate);
      volatility.push_back(pf[i].volatility);
      time.push_back(pf[i].time);
      option_type.push_back((pf[i].option_type == 'P') ? 1 : 0);
    }
  }

  [[nodiscard]] std::vector<Number> compute_prices() const
  {
    const std::size_t size = spot_price.size();
    std::vector<Number> result;
    result.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
      result.push_back(BlkSchlsEqEuroNoDiv(spot_price[i], strike_price[i],
          rate[i], volatility[i], time[i], option_type[i]));
    }
    return result;
  }

private:
  void reserve(std::size_t n)
  {
    spot_price.reserve(n);
    strike_price.reserve(n);
    rate.reserve(n);
    volatility.reserve(n);
    time.reserve(n);
    option_type.reserve(n);
  }

private:
  std::vector<Number> spot_price;
  std::vector<Number> strike_price;
  std::vector<Number> rate;
  std::vector<Number> volatility;
  std::vector<Number> time;
  std::vector<int8_t> option_type;
};

template<typename Number>
std::vector<Number> compute_values(compute_option_portfolio<Number> & options) {
  std::vector<Number> result;
  for (int j = 0; j < NUM_RUNS; j++) {
    result = options.compute_prices();
  }
  return result;
}

template<typename Number>
void write_prices(const std::string & outputFile, std::vector<Number> & prices) {
  std::ofstream output{outputFile};
  if (!output) { throw invalid_io{outputFile}; }
  fmt::print(output,"{}\n",prices.size());
  if (!output) { throw invalid_io{outputFile}; }
  for (const auto & p : prices) {
    fmt::print(output,"{:.18f}\n",p);
  }
  if (!output) { throw invalid_io{outputFile}; }
}

int main(int argc, char *argv[])
try
{
  using namespace std::chrono;

  auto timer1 = high_resolution_clock::now();

  std::ios_base::sync_with_stdio(false);

  using fptype = float;

  std::cout << "PARSEC Benchmark Suite" << std::endl;

  if (argc != 4) {
    std::cerr << "Usage:\n\t" << argv[0] // NOLINT
        << " <nthreads> <inputFile> <outputFile>\n";
    return 1;
  }
  int nThreads = std::stoi(argv[1]); // NOLINT
  if (nThreads != 1) {
    std::cerr << "Error: <nthreads> must be 1 (serial version)\n";
    exit(1);
  }
  std::string inputFile = argv[2]; // NOLINT
  std::string outputFile = argv[3]; // NOLINT

  //Read input data from file
  auto input_portfolio = read_input_portfolio<fptype>(inputFile);
  std::cout << "Num of Options: " << input_portfolio.size() << "\n";
  std::cout << "Num of Runs: " << NUM_RUNS << "\n";

  auto timer2 = high_resolution_clock::now();

  // Compute portfolio prices
  compute_option_portfolio<fptype> compute_porfolio{input_portfolio};
  std::cout << "Size of data: " << input_portfolio.data_size() << "\n";
  auto prices = compute_values(compute_porfolio);

  auto timer3 = high_resolution_clock::now();

  //Write prices to output file
  write_prices(outputFile, prices);

  auto timer4 = high_resolution_clock::now();

  auto d1 = duration_cast<microseconds>(timer2 - timer1);
  auto d2 = duration_cast<microseconds>(timer3 - timer2);
  auto d3 = duration_cast<microseconds>(timer4 - timer2);

  std::cout << "Reading time: " << static_cast<double>(d1.count()) / 1000.0 << "\n";
  std::cout << "Processing time: " << static_cast<double>(d2.count()) / 1000.0 << "\n";
  std::cout << "Writing time: " << static_cast<double>(d3.count()) / 1000.0 << "\n";

  return 0;
}
catch (invalid_io & e) {
  std::cerr << e.what() << "\n";
}
catch (...) {
  std::cerr << "Unexpected exception\n";
}
