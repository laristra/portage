#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>

const double epsilon = 1.e-12;

// from https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
std::vector<std::string> getLineAndSplit(std::istream& stream)
{
    std::vector<std::string>   result;
    std::string                line;
    std::getline(stream, line);
    if (line.empty()) return result;

    std::stringstream          lineStream(line, std::ios_base::in);
    std::string                field;

    while(std::getline(lineStream,field, ','))
    {
        result.push_back(field);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && field.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}

int main(int narg, char** argv) {
  if (narg < 3) {
    std::cout << "not enough arguments\n";
    exit(5);
  }

  std::string goldfile, testfile;
  goldfile = std::string(argv[1]);
  testfile = std::string(argv[2]);

  std::vector<std::vector<double>> goldcoords(0), testcoords(0);
  std::vector<double> goldvalues(0), testvalues(0);
  
  std::ifstream goldin, testin;
  goldin.open(goldfile);
  testin.open(testfile);

  int ncol = 0;
  bool goon=true;
  getLineAndSplit(goldin); // skip gold header
  while (goon) {
    std::vector<std::string> values = getLineAndSplit(goldin);
    if (values.empty()) { // EOF
      goon = false;
    } else if (values.size() == 3) { // NORMAL
      if (ncol == 0) ncol = 3;
      std::vector<double> newcoord = {stod(values[0]), stod(values[1])};
      goldcoords.push_back(newcoord);
      goldvalues.push_back(stod(values[2]));
    } else if (values.size() == 4) { // NORMAL
      if (ncol == 0) ncol = 4;
      std::vector<double> newcoord = {stod(values[0]), stod(values[1]), stod(values[2])};
      goldcoords.push_back(newcoord);
      goldvalues.push_back(stod(values[3]));
    }
    if (values.size() != ncol and ncol > 0 and goon) { // wrong number of columns
      std::cout << "uneven columns in gold file\n";
      std::exit(1);
    }
  }
  int ncolgold = ncol;

  ncol = 0;
  goon = true;
  getLineAndSplit(testin); // skip test header
  while (goon) {
    std::vector<std::string> values = getLineAndSplit(testin);
    if (values.empty()) { // EOF
      goon = false;
    } else if (values.size() == 3) { // NORMAL
      if (ncol == 0) ncol = 3;
      std::vector<double> newcoord = {stod(values[0]), stod(values[1])};
      testcoords.push_back(newcoord);
      testvalues.push_back(stod(values[2]));
    } else if (values.size() == 4) { // NORMAL
      if (ncol == 0) ncol = 4;
      std::vector<double> newcoord = {stod(values[0]), stod(values[1]), stod(values[2])};
      testcoords.push_back(newcoord);
      testvalues.push_back(stod(values[3]));
    }
    if (values.size() != ncol and ncol > 0 and goon) { // wrong number of columns
      std::cout << "uneven columns in test file\n";
      std::exit(1);
    }
  }

  if (goldcoords.size() != testcoords.size() or 
      goldcoords.size() != goldvalues.size() or 
      testcoords.size() != testvalues.size() or 
      ncol != ncolgold) 
  {
    std::cout << "length of gold and test data not equal\n";
    std::exit(3);
  }

  bool isgood = true;
  for (size_t i=0; i<goldcoords.size(); i++) {
    bool isgood0, isgood1, isgood2=true, isgood3;
    isgood0 = fabs(goldcoords[i][0] - testcoords[i][0]) < epsilon;
    isgood1 = fabs(goldcoords[i][1] - testcoords[i][1]) < epsilon;
    if (ncol == 4) {
      isgood2 = fabs(goldcoords[i][2] - testcoords[i][2]) < epsilon;
    }
    isgood3 = fabs(goldvalues[i] - testvalues[i]) < epsilon;
    isgood = isgood0 and isgood1 and isgood2 and isgood3;
    if (not isgood0 or not isgood1 or not isgood2) {
      std::cout << "coordinates disagree: at line " << i << std::endl;
    }
    if (not isgood3) {
      std::cout << "values disagree at line " << i << std::endl;
    }
    isgood = isgood and isgood1 and isgood2;
  }

  if (isgood) {
    std::cout << "files agree" << std::endl;
  }
}
