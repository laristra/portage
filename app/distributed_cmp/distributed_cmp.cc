/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <sstream>
#include <map>
/*
#include <sys/time.h>

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <memory>
#include <cmath>
*/

void print_usage() {
  std::printf("Usage: distributed_cmp file_serial\n");
}

int main(int argc, char** argv) {

  if (argc !=2 ) {
    print_usage();
    return 1;
  }
  
  // vars for unpacking the file  
  int gid,matid;
  double value;
  
  std::pair<int,int> key;
  std::map<std::pair<int,int>,double> serial_map, distributed_map;
  
  // base filename
  std::string base_filename=argv[1];
  
  // error string
  std::stringstream error_string;
  
  /////////////////////////////// 
  // process the serial file
  /////////////////////////////// 
  
  std::ifstream f;
  
  f.open(base_filename);
  
  // make sure the serial file is good
  if (!f) throw std::runtime_error("Serial file " + base_filename +  " cannot be opened.");
  
  // read the file and unpack the tokens
  while (f >> gid >> matid >> value){
  
    std::cout << gid << " " << matid << " " << value << "\n";
    
    // construct the pair key
    key = std::make_pair(gid,matid);
    
    // does the key already exist
    if (serial_map.find(key)!=serial_map.end())
    {
      error_string << "Serial file " << base_filename <<  
        " had a duplicated key (gid, matid)=(" << gid << ", " << matid << ")";      
      throw std::runtime_error(error_string.str());
    }
    
    // add the key to the map
    serial_map[key]=value;
    
  }

  f.close();
  
  // loop through partitions
  int rank=0;
  
  while (true) {
  
    std::cout <<"\n";
    
    // open the serial file
    f.open(base_filename + "." + std::to_string(rank));
  
    // make sure the serial file is good
    if (!f) break;
  
    // read the file and unpack the tokens
    while (f >> gid >> matid >> value){
    
      std::cout<<gid<<" "<<matid<<" "<<value<<"\n";
      
      // construct the pair key
      key = std::make_pair(gid,matid);
      
      // check the distributed data is correct
      if (serial_map.find(key)==serial_map.end())
      {      
        // the key is not found in the serial map  
        error_string << "Distributed file: " << base_filename << "." << rank <<
          " had a key (gid, matid): (" << gid << ", " << matid << 
          ") not in the serial file";      
        throw std::runtime_error(error_string.str());
      } 
      else if (serial_map.find(key)!=serial_map.end() &&
               value != serial_map[key])
      {      
        // the key is already registered but with a different value  
        error_string << "Distributed file: " << base_filename << "." << rank <<
          " had a conflicting key: (" << gid << ", " << matid << 
          "). Serial value: " << serial_map[key] << " Distributed value: " << value;      
        throw std::runtime_error(error_string.str());
      }
      else
      {   
        // key is good, add the key to the map
        distributed_map[key]=value;
      }
      
    }
  
    // clean up
    f.close();
    
    // bump the partition
    rank++;
    
  }

  // check that there are distributed data files
  if (rank==0) 
    throw std::runtime_error("No partitions were found for  " + base_filename);
    
  // make sure there are no missing keys in the distributed files
  if (distributed_map.size() < serial_map.size()) 
    throw std::runtime_error("The distributed files missed keys from the serial run");

  return 0;
}
