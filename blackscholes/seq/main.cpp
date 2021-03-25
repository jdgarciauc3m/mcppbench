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
#include <cstring>

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

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
template <typename fptype>
fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
    fptype strike, fptype rate, fptype volatility,
    fptype time, int option_type)
{
  fptype OptionPrice;

  // local private working variables for the calculation
  fptype xStockPrice;
  fptype xStrikePrice;
  fptype xRiskFreeRate;
  fptype xVolatility;
  fptype xTime;
  fptype xSqrtTime;

  fptype logValues;
  fptype xLogTerm;
  fptype xD1;
  fptype xD2;
  fptype xPowerTerm;
  fptype xDen;
  fptype d1;
  fptype d2;
  fptype FutureValueX;
  fptype NofXd1;
  fptype NofXd2;
  fptype NegNofXd1;
  fptype NegNofXd2;

  xStockPrice = sptprice;
  xStrikePrice = strike;
  xRiskFreeRate = rate;
  xVolatility = volatility;

  xTime = time;
  xSqrtTime = sqrt(xTime);

  logValues = std::log( sptprice / strike );

  xLogTerm = logValues;


  xPowerTerm = xVolatility * xVolatility;
  xPowerTerm = xPowerTerm * 0.5;

  xD1 = xRiskFreeRate + xPowerTerm;
  xD1 = xD1 * xTime;
  xD1 = xD1 + xLogTerm;

  xDen = xVolatility * xSqrtTime;
  xD1 = xD1 / xDen;
  xD2 = xD1 -  xDen;

  d1 = xD1;
  d2 = xD2;

  NofXd1 = cdf_normal(d1);
  NofXd2 = cdf_normal(d2);

  FutureValueX = strike * ( std::exp(-(rate)*(time) ) );
  if (option_type == 0) {
    OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
  } else {
    NegNofXd1 = (1.0 - NofXd1);
    NegNofXd2 = (1.0 - NofXd2);
    OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
  }

  return OptionPrice;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

template <typename fptype>
int bs_thread(void *tid_ptr) {
  int i, j;
  fptype price;
  fptype priceDelta;
  int tid = *(int *)tid_ptr;
  int start = tid * (numOptions / nThreads);
  int end = start + (numOptions / nThreads);

  for (j=0; j<NUM_RUNS; j++) {
    for (i=start; i<end; i++) {
      /* Calling main function to calculate option value based on
       * Black & Scholes's equation.
       */
      price = BlkSchlsEqEuroNoDiv( sptprice<fptype>[i], strike<fptype>[i],
          rate<fptype>[i], volatility<fptype>[i], otime<fptype>[i],
          otype[i]);
      prices<fptype>[i] = price;

#ifdef ERR_CHK
      priceDelta = data[i].DGrefval - price;
            if( fabs(priceDelta) >= 1e-4 ){
                printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
                       i, price, data[i].DGrefval, priceDelta);
                numError ++;
            }
#endif
    }
  }

  return 0;
}

int main (int argc, char **argv)
{
  FILE *file;
  int i;
  int loopnum;
  using fptype = float;
  fptype * buffer;
  int * buffer2;
  int rv;

  printf("PARSEC Benchmark Suite\n");
  fflush(NULL);

  if (argc != 4)
  {
    printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
    exit(1);
  }
  nThreads = atoi(argv[1]);
  char *inputFile = argv[2];
  char *outputFile = argv[3];

  //Read input data from file
  file = fopen(inputFile, "r");
  if(file == NULL) {
    printf("ERROR: Unable to open file `%s'.\n", inputFile);
    exit(1);
  }
  rv = fscanf(file, "%i", &numOptions);
  if(rv != 1) {
    printf("ERROR: Unable to read from file `%s'.\n", inputFile);
    fclose(file);
    exit(1);
  }
  if(nThreads > numOptions) {
    printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");
    nThreads = numOptions;
  }

  if(nThreads != 1) {
    printf("Error: <nthreads> must be 1 (serial version)\n");
    exit(1);
  }

  // alloc spaces for the option data
  data<fptype> = (OptionData<fptype>*)malloc(numOptions*sizeof(OptionData<fptype>));
  prices<fptype> = (fptype*)malloc(numOptions*sizeof(fptype));
  for ( loopnum = 0; loopnum < numOptions; ++ loopnum )
  {
    rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data<fptype>[loopnum].s, &data<fptype>[loopnum].strike, &data<fptype>[loopnum].r, &data<fptype>[loopnum].divq, &data<fptype>[loopnum].v, &data<fptype>[loopnum].t, &data<fptype>[loopnum].OptionType, &data<fptype>[loopnum].divs, &data<fptype>[loopnum].DGrefval);
    if(rv != 9) {
      printf("ERROR: Unable to read from file `%s'.\n", inputFile);
      fclose(file);
      exit(1);
    }
  }
  rv = fclose(file);
  if(rv != 0) {
    printf("ERROR: Unable to close file `%s'.\n", inputFile);
    exit(1);
  }

  printf("Num of Options: %d\n", numOptions);
  printf("Num of Runs: %d\n", NUM_RUNS);

  constexpr int PAD = 256;
  constexpr int LINESIZE = 64;

  buffer = (fptype *) malloc(5 * numOptions * sizeof(fptype) + PAD);
  sptprice<fptype> = (fptype *) (((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
  strike<fptype> = sptprice<fptype> + numOptions;
  rate<fptype> = strike<fptype> + numOptions;
  volatility<fptype> = rate<fptype> + numOptions;
  otime<fptype> = volatility<fptype> + numOptions;

  buffer2 = (int *) malloc(numOptions * sizeof(fptype) + PAD);
  otype = (int *) (((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

  for (i=0; i<numOptions; i++) {
    otype[i]      = (data<fptype>[i].OptionType == 'P') ? 1 : 0;
    sptprice<fptype>[i]   = data<fptype>[i].s;
    strike<fptype>[i]     = data<fptype>[i].strike;
    rate<fptype>[i]       = data<fptype>[i].r;
    volatility<fptype>[i] = data<fptype>[i].v;
    otime<fptype>[i]      = data<fptype>[i].t;
  }

  printf("Size of data: %d\n", numOptions * (sizeof(OptionData<fptype>) + sizeof(int)));

  //serial version
  int tid=0;
  bs_thread<fptype>(&tid);

  //Write prices to output file
  file = fopen(outputFile, "w");
  if(file == NULL) {
    printf("ERROR: Unable to open file `%s'.\n", outputFile);
    exit(1);
  }
  rv = fprintf(file, "%i\n", numOptions);
  if(rv < 0) {
    printf("ERROR: Unable to write to file `%s'.\n", outputFile);
    fclose(file);
    exit(1);
  }
  for(i=0; i<numOptions; i++) {
    rv = fprintf(file, "%.18f\n", prices<fptype>[i]);
    if(rv < 0) {
      printf("ERROR: Unable to write to file `%s'.\n", outputFile);
      fclose(file);
      exit(1);
    }
  }
  rv = fclose(file);
  if(rv != 0) {
    printf("ERROR: Unable to close file `%s'.\n", outputFile);
    exit(1);
  }

#ifdef ERR_CHK
  printf("Num Errors: %d\n", numError);
#endif
  free(data<fptype>);
  free(prices<fptype>);

  return 0;
}

