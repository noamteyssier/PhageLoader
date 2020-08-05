#include <iostream>
#include <armadillo>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
po::variables_map process_program_options(const int argc, const char *const argv[])
{
  po::variables_map vm;

  try
  {
    po::options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help screen")
      ("input_csv,i", po::value<std::string>()->required(), "Input CSV of Count Matrix to Preprocess")
      ("output_arma,o", po::value<std::string>()->required(), "Output arma to save to");

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
      std::cout << desc << '\n';
      exit(-1);
    }
    else if (!vm.count("input_csv"))
    {
      std::cout
        << std::endl
        << "Error : Input Matrix Required. (--input_csv,-i)"
        << std::endl << std::endl;

      std::cout << desc << '\n';
      exit(-1);
    }
    else if (!vm.count("output_arma"))
    {
      std::cout
        << std::endl
        << "Error : Output Filename Required. (--output_arma,-o)"
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

  // assign arguments
  std::string raw_csv = vm["input_csv"].as<std::string>();
  std::string output_arma = vm["output_arma"].as<std::string>();

  arma::mat m;
  m.load(raw_csv, arma::csv_ascii);
  m.save(output_arma, arma::arma_binary);
  return 0;
}
