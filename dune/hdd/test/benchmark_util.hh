//
// Created by r_milk01 on 21.02.17.
//

#ifndef DUNE_HDD_BNECHMARK_UTIL_HH
#define DUNE_HDD_BNECHMARK_UTIL_HH

#include <dune/common/parallel/mpihelper.hh>

#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/string.hh>
#include <boost/filesystem.hpp>

namespace Dune {
namespace Stuff {

template <typename T>
class PrefixOutputIterator
{
  std::ostream& ostream;
  std::string prefix;
  bool first;

public:
  typedef std::size_t difference_type;
  typedef T value_type;
  typedef T* pointer;
  typedef T reference;
  typedef std::output_iterator_tag iterator_category;

  PrefixOutputIterator(std::ostream& o, std::string const& p = "")
      : ostream(o)
      , prefix(p)
      , first(true)
  {
  }

  PrefixOutputIterator& operator*()
  {
    return *this;
  }
  PrefixOutputIterator& operator++()
  {
    return *this;
  }
  PrefixOutputIterator& operator++(int)
  {
    return *this;
  }

  void operator=(T const& value)
  {
    if (first) {
      ostream << value;
      first = false;
    } else {
      ostream << prefix << value;
    }
  }
};


void dump_environment(boost::filesystem::ofstream& file, std::string csv_sep = ",")
{
  using namespace std;
  vector<string> header, values;
  for (char** current = environ; *current; current++) {
    string line(*current);
    const auto tokens = DSC::tokenize(line, "=");
    if (tokens.size() == 2) {
      header.push_back(tokens[0]);
      values.push_back(tokens[1]);
    }
  }
  copy(header.begin(), header.end(), PrefixOutputIterator<string>(file, csv_sep));
  file << '\n';
  copy(values.begin(), values.end(), PrefixOutputIterator<string>(file, csv_sep));
}

void mem_usage() {
  auto comm = Dune::MPIHelper::getCollectiveCommunication();
  // Compute the peak memory consumption of each processes
  int who = RUSAGE_SELF;
  struct rusage usage;
  getrusage(who, &usage);
  long peakMemConsumption = usage.ru_maxrss;
  // compute the maximum and mean peak memory consumption over all processes
  long maxPeakMemConsumption = comm.max(peakMemConsumption);
  long totalPeakMemConsumption = comm.sum(peakMemConsumption);
  long meanPeakMemConsumption = totalPeakMemConsumption / comm.size();
  // write output on rank zero
  if (comm.rank() == 0) {
    std::unique_ptr<boost::filesystem::ofstream> memoryConsFile(DSC::make_ofstream(
        std::string(DSC_CONFIG_GET("global.datadir", "data/")) + std::string("/memory.csv")));
    *memoryConsFile << "global.maxPeakMemoryConsumption,global.meanPeakMemoryConsumption\n" << maxPeakMemConsumption
                    << "," << meanPeakMemConsumption << std::endl;
  }
}

void dump_environment() {
  std::unique_ptr<boost::filesystem::ofstream> of(DSC::make_ofstream(
      std::string(DSC_CONFIG_GET("global.datadir", "data/")) + std::string("/env.txt")));
  dump_environment(*of);
}
int handle_exception(const Dune::Exception& exp) {
  std::cerr << "Failed with Dune::Exception: " << exp.what();
  DSC_PROFILER.output_per_rank("profiler");
  mem_usage();
  return Dune::Stuff::abort_all_mpi_processes();
}

int handle_exception(const std::exception& exp) {
  std::cerr << "Failed with std::exception: " << exp.what();
  DSC_PROFILER.output_per_rank("profiler");
  mem_usage();
  return Dune::Stuff::abort_all_mpi_processes();
}

int handle_exception(const tbb::tbb_exception& exp) {
  std::cerr << "Failed with tbb::exception" << exp.name() << ": " << exp.what();
  DSC_PROFILER.output_per_rank("profiler");
  mem_usage();
  return Dune::Stuff::abort_all_mpi_processes();
}

}
}
#endif //DUNE_HDD_BNECHMARK_UTIL_HH
