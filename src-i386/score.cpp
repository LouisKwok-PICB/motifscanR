#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List sliding_motif_score(Rcpp::NumericMatrix matrix, unsigned int matlen,
                               const double scoremax,
                               const std::vector<std::string> sequences){
  std::string seq;
  size_t nseqs = sequences.size();
  Rcpp::List sldscores(nseqs);
  for (size_t seqidx = 0; seqidx < nseqs; seqidx++)
  {
    seq = sequences[seqidx];
    size_t seqlen = seq.size();
    if (seqlen < matlen){
      throw std::range_error("sequence length should >= matrix length");
    }
    std::vector<size_t> rowidx(seqlen);
    for (size_t pos = 0; pos < seqlen; pos++)
    {
      switch (seq[pos])
      {
      case 'A':
      case 'a':
        rowidx[pos] = 0;
        break;
      case 'C':
      case 'c':
        rowidx[pos] = 1;
        break;
      case 'G':
      case 'g':
        rowidx[pos] = 2;
        break;
      case 'T':
      case 't':
        rowidx[pos] = 3;
        break;
      default:
        rowidx[pos] = -1;
      break;
      }
    }
    size_t scorelen = seqlen - matlen + 1;
    std::vector<double> sldscore(scorelen);
    for (size_t i = 0; i < scorelen; i++)
    {
      double score = 0;
      for (size_t j = 0; j < matlen; j++)
      {
        size_t pos = i + j;
        int rowidx_in = rowidx[pos];
        if (rowidx_in > -1)
        {
          score = score + matrix(rowidx_in, j);
        }
      }
      sldscore[i] = score / scoremax;
    }
    sldscores[seqidx] = sldscore;
  }
  return sldscores;
}

// [[Rcpp::export]]
std::vector<double> motif_score(Rcpp::NumericMatrix matrix, unsigned int matlen,
                                const double scoremax,
                                const std::vector<std::string> sequences){
  std::string seq;
  size_t nseqs = sequences.size();
  std::vector<double> scores(nseqs);

  for (size_t seqidx = 0; seqidx < nseqs; seqidx++)
  {
    seq = sequences[seqidx];
    size_t seqlen = seq.size();
    if (seqlen < matlen){
      throw std::range_error("sequence length should >= matrix length");
    }
    double score = 0;
    for (size_t pos = 0; pos < matlen; pos++)
    {
      switch (seq[pos])
      {
      case 'A':
      case 'a':
        score = score + matrix(0, pos);
        break;
      case 'C':
      case 'c':
        score = score + matrix(1, pos);
        break;
      case 'G':
      case 'g':
        score = score + matrix(2, pos);
        break;
      case 'T':
      case 't':
        score = score + matrix(3, pos);
        break;
      default:
        break;
      }
    }
    scores[seqidx] = score / scoremax;
  }
  return scores;
}
