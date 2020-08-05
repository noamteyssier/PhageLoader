
#include <iostream>
#include <armadillo>
#include <boost/program_options.hpp>

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
  void MinimumHitThreshold(int threshold, bool all_samples)
  {
    arma::colvec peptide_sums;
    if (all_samples)
    {
      peptide_sums = arma::sum(bool_matrix, 1);
    }
    else
    {
      arma::uvec patient_indices = sl.GetPatientIndices();
      peptide_sums = arma::sum(bool_matrix.cols(patient_indices), 1);
    }

    passing_peptides = arma::find(peptide_sums >= threshold);
  }

  void Write_BMAT(std::string output_prefix)
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
    arma::conv_to<arma::imat>::from(filtered_bool).save(arma::csv_name(output_prefix + "_bool.csv", header));

  }

  void Write_ZMAT(std::string output_prefix)
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
    filtered_zscore.save(arma::csv_name(output_prefix + "_zscore.csv", header));

  }

};


namespace po = boost::program_options;
po::variables_map process_program_options(const int argc, const char *const argv[])
{
  po::variables_map vm;
  bool flag = false;
  try
  {
    po::options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("input_arma,i", po::value<std::string>()->required(), "Input Count Matrix to Preprocess in ARMA binary format")
      ("sample_names,n", po::value<std::string>()->required(), "Sample Names of columns in count matrix")
      ("output_prefix,o", po::value<std::string>()->required(), "Output prefix to save to")
      ("minimum_readcount,r", po::value<int>()->default_value(250000), "Minimum read count to consider a sample (default = 250000)")
      ("scalar,s", po::value<int>()->default_value(500000), "Scalar to normalize read counts to post filtering (default = 500000)")
      ("z_min,z", po::value<int>()->default_value(10), "Minimum Z-Score to accept a hit (default = 10)")
      ("c_min,c", po::value<int>()->default_value(8), "Minimum number of hits in sample set to accept a peptide as enriched (default = 8)")
      ("all_samples,p", po::bool_switch(&flag), "Include all samples in hit count threshold (default = Only Patient Samples)");

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout << desc << '\n';
      exit(-1);
    }
    else if (!vm.count("input_arma"))
    {
      std::cout
        << std::endl
        << "Error : Input Matrix in ARMA binary format required. (--input_arma,-i)"
        << std::endl << std::endl;

      std::cout << desc << '\n';
      exit(-1);
    }
    else if (!vm.count("sample_names"))
    {
      std::cout
        << std::endl
        << "Error : Sample names list required. (--sample_names,-n)"
        << std::endl << std::endl;

      std::cout << desc << '\n';
      exit(-1);
    }
    else if (!vm.count("output_prefix"))
    {
      std::cout
        << std::endl
        << "Error : Output Filename Prefix Required. (--output_prefix,-o)"
        << std::endl << std::endl;

      std::cout << desc << '\n';
      exit(-1);
    }
    po::notify(vm);
  }
  catch (const po::error &ex)
  {
    std::cerr << ex.what() << '\n';
  }

  return vm;
}

int main(const int argc, const char *const argv[])
{

  po::variables_map vm = process_program_options(argc, argv);

  std::string arma_fn = vm["input_arma"].as<std::string>();
  std::string sample_fn = vm["sample_names"].as<std::string>();
  std::string output_prefix = vm["output_prefix"].as<std::string>();

  int minimum_readcount = vm["minimum_readcount"].as<int>();
  int scalar = vm["scalar"].as<int>();
  int z_min = vm["z_min"].as<int>();
  int c_min = vm["c_min"].as<int>();
  bool all_samples = vm["all_samples"].as<bool>();


  // read in arma
  arma::mat m;
  m.load(arma_fn, arma::arma_binary);

  // Load MatrixLoader Object
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
  ml.MinimumHitThreshold(c_min, all_samples);

  // write enriched peptide set as a bool matrix and as a zscore matrix
  ml.Write_BMAT(output_prefix);
  ml.Write_ZMAT(output_prefix);

  return 0;
}
