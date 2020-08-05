
#include <iostream>
#include <armadillo>

template <typename T>
arma::uvec vec2vec(const std::vector<T>& input_vec)
{
  arma::uvec output_vec(input_vec.size());
  for (size_t i = 0; i < input_vec.size(); i++)
  {
    output_vec[i] = input_vec[i];
  }

  return output_vec;
}

class SampleLoader
{
private:
  std::vector<std::string> sample_names;
  std::vector<std::string> parsed_sample_names;
  std::vector<std::string> pair_names;
  std::map<std::string, std::vector<int>> sample_pairs;
  std::map<std::string, int> pair_sizes;

public:

  void init(std::string filename)
  {
    LoadFile(filename);
    ParseNames();
  }

  void LoadFile(std::string filename)
  {
    std::ifstream file(filename);
    std::string line;

    while(getline(file, line))
    {
      sample_names.push_back(line);
    }
  }

  std::vector<std::string> SplitString(const std::string& str, char delim)
  {
    std::stringstream ss(str);

    std::vector<std::string> token_list;
    std::string token;

    while(getline(ss, token, delim))
    {
      token_list.push_back(token);
    }

    return token_list;
  }

  void ParseNames()
  {
    sample_pairs.clear();
    pair_names.clear();
    parsed_sample_names.clear();

    int total_samples = sample_names.size();
    std::vector<std::string> token_list;

    int total_pairs=0;
    for (int i = 0; i < total_samples; i++)
    {
      token_list.clear();
      token_list = SplitString(sample_names[i], '_');
      std::string sample_name = token_list[1];

      if (sample_pairs.find(sample_name) == sample_pairs.end())
      {
        sample_pairs[sample_name];
        pair_names.push_back(sample_name);
        total_pairs++;
      }

      sample_pairs[sample_name].push_back(i);
      parsed_sample_names.push_back(sample_name);
    }

  }

  std::vector<std::string> GetPairNames()
  {
    return pair_names;
  }

  std::vector<std::string> GetSampleNames()
  {
    return sample_names;
  }

  std::vector<std::string> GetParsedSampleNames()
  {
    return parsed_sample_names;
  }

  std::map<std::string, int> GetPairSizes()
  {
    for (auto & x : sample_pairs)
    {
      int size = x.second.size();
      pair_sizes.insert(std::pair<std::string, int>(x.first, size));
    }

    return pair_sizes;
  }

  std::map<std::string, std::vector<int>> GetPairMap()
  {
    return sample_pairs;
  }

  int GetTotalPairs()
  {
    return pair_names.size();
  }

  int GetTotalSamples()
  {
    return sample_names.size();
  }

  void RemoveIndices(arma::uvec indices)
  {
    sort(indices.begin(), indices.end(), std::greater<int>());
    std::string sample_name;


    for (int && x : indices)
    {
      sample_name = parsed_sample_names[x];

      int iter = 0;
      for (auto & y : sample_pairs[sample_name])
      {
        if (y == x)
        {
          sample_pairs[sample_name].erase(sample_pairs[sample_name].begin() + iter);
        }
        iter++;
      }

      // clears empty pair name lists
      if (sample_pairs[sample_name].size() < 1)
      {
        sample_pairs.erase(sample_name);
      }

      sample_names.erase(sample_names.begin() + x);
      parsed_sample_names.erase(parsed_sample_names.begin() + x);
    }

    // recalculate indices
    ParseNames();
  }

  arma::uvec GetHealthyIndices()
  {
    std::vector<int> healthy_indices;

    for (auto & x : sample_pairs)
    {
      if (x.first.find("HC") != std::string::npos)
      {
        for (auto & y : x.second)
        {
          healthy_indices.push_back(y);
        }
      }
    }

    return vec2vec(healthy_indices);

  }

  arma::uvec GetPatientIndices()
  {
    std::vector<int> patient_indices;

    for (size_t i = 0; i < pair_names.size(); i++)
    {
      // kanungu
      if (pair_names[i].find("CK") != std::string::npos)
      {
        patient_indices.push_back(i);
      }

      // tororo
      if (pair_names[i].find("CT") != std::string::npos)
      {
        patient_indices.push_back(i);
      }

    }

    return vec2vec(patient_indices);
  }



};


class MatrixLoader
{
private:
  SampleLoader sl;
  arma::mat count_matrix;
  arma::mat zscore_matrix;

  arma::mat paired_zscore_matrix;
  arma::mat bool_matrix;

  arma::mat filtered_bool;
  arma::mat filtered_zscore;

  arma::rowvec sample_sums;
  arma::uvec failing_sample_indices;

  arma::colvec healthy_means;
  arma::colvec healthy_sds;

  arma::uvec passing_peptides;

public:
  MatrixLoader(arma::mat m, std::string sample_fn) :
    count_matrix(m)
    {
      sl.init(sample_fn);
      ConfirmDimensions();
    };

  void ConfirmDimensions()
  {
    int num_samples = sl.GetTotalSamples();
    int num_cols = count_matrix.n_cols;

    if (num_samples != num_cols)
    {
      fprintf(stderr, "Error : Expecting (%i) Samples but Sample List has (%lli)\n", num_samples, count_matrix.n_cols);
      exit(-1);
    }

  }

  void DiscardSamples(arma::uvec indices)
  {
    fprintf(stderr, "Removing (%lli) samples\n", indices.n_elem);

    count_matrix.shed_cols(indices);
    sample_sums.shed_cols(indices);
    sl.RemoveIndices(indices);
  }

  void FilterSamples(int minimum_readcount)
  {
    sample_sums = arma::sum(count_matrix, 0);
    failing_sample_indices = arma::find(sample_sums < minimum_readcount);
    DiscardSamples(failing_sample_indices);
  }

  // Confirms at least one replicate for each sample pair or drops sample
  void ValidateSamplePairs()
  {

    std::map<std::string, int> pair_sizes = sl.GetPairSizes();
    std::vector<std::string> parsed_sample_names = sl.GetParsedSampleNames();


    std::vector<std::string> pairs_to_remove;
    std::vector<int> indices_to_remove;

    for (auto & x : pair_sizes)
    {
      if (x.second < 2)
      {
        pairs_to_remove.push_back(x.first);
      }
    }

    for (std::string & x : pairs_to_remove)
    {

      for (size_t i = 0; i < parsed_sample_names.size(); i++)
      {
        if (parsed_sample_names[i] == x)
        {
          indices_to_remove.push_back(i);
        }
      }
    }

    arma::uvec indices(indices_to_remove.size());
    for (size_t i = 0; i < indices_to_remove.size(); i++)
    {
      indices[i] = indices_to_remove[i];
    }

    DiscardSamples(indices);
  }

  // Normalize read counts to given scalar
  void NormalizeCounts(int scalar)
  {
    sample_sums.clear();
    sample_sums = arma::sum(count_matrix, 0);

    count_matrix.each_row(
      [this, scalar] (arma::rowvec& x)
      {
        x = (x / sample_sums) * scalar;
      }
    );
  }

  // calculate mean and std on healthys for each target
  void HealthyParams()
  {
    arma::uvec healthy_indices = sl.GetHealthyIndices();

    // subset healthy indices in matrix
    arma::mat healthy_matrix = count_matrix.cols(healthy_indices);
    healthy_means = arma::mean(healthy_matrix, 1);
    healthy_sds = arma::stddev(healthy_matrix, 0, 1);

    healthy_sds.elem( find(healthy_sds < 1) ).ones();
  }

  // calculate zscore on entire set given healthy params
  void CalculateZscore()
  {
    zscore_matrix = count_matrix;

    zscore_matrix.each_col(
      [this] (arma::colvec& x)
      {
        x = ((x - healthy_means) / healthy_sds);
      }
    );

  }

  // cross validate zscore between pairs
  void CrossValidateZscore(int threshold)
  {
    std::map<std::string, std::vector<int>> sample_pairs = sl.GetPairMap();
    std::vector<std::string> pair_names = sl.GetPairNames();


    bool_matrix.resize(zscore_matrix.n_rows, pair_names.size());
    bool_matrix.zeros();

    paired_zscore_matrix.resize(zscore_matrix.n_rows, pair_names.size());
    paired_zscore_matrix.zeros();


    for (size_t i = 0; i < pair_names.size(); i++)
    {
      std::string p = pair_names[i];
      arma::uvec replicate_index = vec2vec(sample_pairs[p]);

      arma::mat submatrix = zscore_matrix.cols ( replicate_index );

      int row_counter=0;
      submatrix.each_row (
        [this, threshold, i, &row_counter] (arma::rowvec& x)
        {
          arma::uvec passing_indices = arma::find( x >= threshold);
          if (passing_indices.size() > 1)
          {
            bool_matrix(row_counter, i) = 1;
            paired_zscore_matrix(row_counter, i) = arma::mean( x );
          }

          row_counter++;
        }
      );

    }

  }

  // only accepts peptides passing a minimum count threshold
  void MinimumHitThreshold(int threshold, bool only_patients)
  {
    arma::colvec peptide_sums;
    if (only_patients)
    {
      arma::uvec patient_indices = sl.GetPatientIndices();
      peptide_sums = arma::sum(bool_matrix.cols(patient_indices), 1);
    }
    else
    {
      peptide_sums = arma::sum(bool_matrix, 1);
    }

    passing_peptides = arma::find(peptide_sums >= threshold);
  }

  void Write_BMAT()
  {

    filtered_bool = bool_matrix.rows(passing_peptides);

    // build sample name headers
    std::vector<std::string> pair_names = sl.GetPairNames();
    arma::field<std::string> header(pair_names.size());
    for (size_t i = 0; i < pair_names.size(); i++)
    {
      header[i] = pair_names[i];
    }


    // convert to integer matrix and save to output
    arma::conv_to<arma::imat>::from(filtered_bool).save(arma::csv_name("output_bool.csv", header));

  }

  void Write_ZMAT()
  {

    filtered_zscore = paired_zscore_matrix.rows(passing_peptides);

    // build sample name headers
    std::vector<std::string> pair_names = sl.GetPairNames();
    arma::field<std::string> header(pair_names.size());
    for (size_t i = 0; i < pair_names.size(); i++)
    {
      header[i] = pair_names[i];
    }

    // convert to integer matrix and save to output
    filtered_zscore.save(arma::csv_name("output_zscore.csv", header));

  }

};


int main()
{
  int minimum_readcount = 250000;
  int scalar = 500000;
  int z_min = 10;
  int c_min = 8;
  bool only_patients = true;
  arma::mat m;
  std::string sample_fn = "../data/sample_names.txt";


  m.load("../data/count_matrix.arma", arma::arma_binary);

  MatrixLoader ml(m, sample_fn);

  // filter samples with a minimum read count
  ml.FilterSamples(minimum_readcount);

  // confirm at least 1 technical replicate per sample
  ml.ValidateSamplePairs();

  // normalize counts to a given scalar
  ml.NormalizeCounts(scalar);

  // calculate healthy mean and std over each peptide
  ml.HealthyParams();

  // calculate zscore enrichment over all peptides
  ml.CalculateZscore();

  // cross validate hits (zscore above threshold) between replicates
  ml.CrossValidateZscore(z_min);

  // Apply a minimum hit threshold over the entire sample set or only over the
  // patient set
  ml.MinimumHitThreshold(c_min, only_patients);

  // write enriched peptide set as a bool matrix and as a zscore matrix
  ml.Write_BMAT();
  ml.Write_ZMAT();

  return 0;
}


// int main()
// {
//   arma::mat m;
//   m.load("../data/count_matrix.csv", arma::csv_ascii);
//   m.save("../data/count_matrix.arma", arma::arma_binary);
// }
