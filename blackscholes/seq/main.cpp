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

#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <chrono>
#include <fmt/os.h>
#include <fmt/chrono.h>
#include <fmt/ostream.h>

constexpr int NUM_RUNS = 100;

class invalid_io {
public:
  explicit invalid_io(std::string file_name)
      : file_name_{std::move(file_name)} {}

  [[nodiscard]] std::string what() const
  {
    return "ERROR: Unable to access file `" + file_name_ + "'.";
  }

private:
  std::string file_name_;
};

class invalid_file_io {
};

template<typename T>
class dynamic_array {
public:

  dynamic_array() = default;

  explicit dynamic_array(std::size_t n)
      :
      size_{n}, buffer_{new T[n]} {}

  [[nodiscard]] std::size_t size() const { return size_; }

  T & operator[](std::size_t i) { return buffer_[i]; }

  const T & operator[](std::size_t i) const { return buffer_[i]; }

  auto begin() { return buffer_.get(); }
  auto begin() const { return buffer_.get(); }
  auto end() { return buffer_.get() +  size_; }
  auto end() const { return buffer_.get() +  size_; }

private:
  std::size_t size_ = 0;
  std::unique_ptr<T[]> buffer_ = nullptr;
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

  void read(std::FILE *file);
};

template<typename Number>
void option<Number>::read(FILE *file)
{
  int ret = fscanf(file, "%f %f %f %f %f %f %c %f %f",
      &spot_price, &strike_price, &rate, &dividend_rate, &volatility, &time,
      &option_type, &divs, &DGrefval);
  if (ret != 9) {
    throw invalid_file_io{};
  }
}

template<typename Name>
class input_option_portfolio;

template<typename Number>
std::istream & operator>>(std::istream & is, input_option_portfolio<Number> & pf);

template<typename Number>
class input_option_portfolio {
public:
  [[nodiscard]] std::size_t size() const { return size_; }

  [[nodiscard]] std::size_t data_size() const {
    return size_ * (sizeof(Number) + sizeof(option<Number>));
  }

  const option<Number> & operator[](std::size_t i) const { return options_[i]; }

  void read(FILE *file);

private:
  std::size_t size_;
  std::unique_ptr<option<Number>[]> options_;
};

template<typename Number>
void input_option_portfolio<Number>::read(FILE *file)
{
  int ret = fscanf(file, "%lu", &size_);
  if (ret != 1) {
    throw invalid_file_io{};
  }
  options_.reset(new option<Number>[size_]);
  std::for_each_n(options_.get(), size_, [&](auto & o) { o.read(file); });
}

template<typename Number>
input_option_portfolio<Number> read_input_portfolio_f(const std::string & file_name)
{
  try {
    FILE *file = fopen(file_name.c_str(), "r");
    if (!file) { throw invalid_io{file_name}; }
    input_option_portfolio<Number> input_portfolio;
    input_portfolio.read(file);
    return input_portfolio;
  }
  catch (invalid_file_io) {
    throw invalid_io{file_name};
  }
}

template<class Number>
class compute_option_portfolio {
public:
  explicit compute_option_portfolio(const input_option_portfolio<Number> & pf) :
    num_options{pf.size()},
    spot_price{new Number[num_options]},
    strike_price{new Number[num_options]},
    rate{new Number[num_options]},
    volatility{new Number[num_options]},
    time{new Number[num_options]},
    option_type{new std::int8_t[num_options]}
{
    for (std::size_t i = 0; i < num_options; ++i) {
      spot_price[i] = pf[i].spot_price;
      strike_price[i] = pf[i].strike_price;
      rate[i] = pf[i].rate;
      volatility[i] = pf[i].volatility;
      time[i] = pf[i].time;
      option_type[i] = pf[i].option_type == 'P' ? 1 : 0;
    }
  }

  [[nodiscard]] dynamic_array<Number> compute_prices() const
  {
    dynamic_array<Number> result{num_options};
    for (std::size_t i = 0; i < num_options; ++i) {
      result[i] = BlkSchlsEqEuroNoDiv(spot_price[i], strike_price[i],
          rate[i], volatility[i], time[i], option_type[i]);
    }
    return result;
  }

private:

private:
  std::size_t num_options;
  std::unique_ptr<Number[]> spot_price;
  std::unique_ptr<Number[]> strike_price;
  std::unique_ptr<Number[]> rate;
  std::unique_ptr<Number[]> volatility;
  std::unique_ptr<Number[]> time;
  std::unique_ptr<int8_t[]> option_type;
};

template<typename Number>
dynamic_array<Number> compute_values(compute_option_portfolio<Number> & options)
{
  dynamic_array<Number> result;
  for (int j = 0; j < NUM_RUNS; j++) {
    result = options.compute_prices();
  }
  return result;
}

template<typename Number>
void write_prices(const std::string & outputFile,
    const dynamic_array<Number> & prices)
{
  auto output = fmt::output_file(outputFile);
  output.print("{}\n", prices.size());
  for (const auto & p : prices) {
    output.print("{:.18f}\n", p);
  }
}

int main(int argc, char *argv[])
try
{
  using namespace std::chrono;

  using fptype = float;

  fmt::print("PARSEC Benchmark Suite\n");

  if (argc != 4) {
    fmt::print(std::cerr, "Usage:\n\t{} <nthreads> <inputFile> <outputFile>\n", argv[0]);
    return 1;
  }
  int nThreads = std::stoi(argv[1]); // NOLINT
  if (nThreads != 1) {
    fmt::print(std::cerr, "Error: <nthreads> must be 1 (serial version)\n");
    exit(1);
  }
  std::string inputFile = argv[2]; // NOLINT
  std::string outputFile = argv[3]; // NOLINT

  auto timer1 = high_resolution_clock::now();

  //Read input data from file
  auto input_portfolio = read_input_portfolio_f<fptype>(inputFile);
  fmt::print("Num of Options: {}\n", input_portfolio.size());
  fmt::print("Num of Runs: {}\n", NUM_RUNS);

  auto timer2 = high_resolution_clock::now();

  // Compute portfolio prices
  compute_option_portfolio<fptype> compute_porfolio{input_portfolio};
  fmt::print("Size of data: {}\n", input_portfolio.data_size());
  auto prices = compute_values(compute_porfolio);

  auto timer3 = high_resolution_clock::now();

  //Write prices to output file
  write_prices(outputFile, prices);

  auto timer4 = high_resolution_clock::now();

  auto d1 = duration_cast<microseconds>(timer2 - timer1);
  auto d2 = duration_cast<microseconds>(timer3 - timer2);
  auto d3 = duration_cast<microseconds>(timer4 - timer3);
  auto d4 = duration_cast<microseconds>(timer4 - timer1);

  fmt::print("Reading time: {} ms\n", d1.count() / 1000.0);
  fmt::print("Processing time: {} ms\n", d2.count() / 1000.0);
  fmt::print("Writing time: {} ms\n", d3.count() / 1000.0);
  fmt::print("Total time: {} ms\n", d4.count() / 1000.0);

  return 0;
}
catch (invalid_io & e) {
  fmt::print(std::cerr,"{}\n", e.what());
}
catch (...) {
  fmt::print(std::cerr,"Unexpected exception\n");
}
