/* Copyright 2017-2021 Fabian Koller
 *
 * This file is part of openPMD-api.
 *
 * openPMD-api is free software: you can redistribute it and/or modify
 * it under the terms of of either the GNU General Public License or
 * the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * openPMD-api is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License and the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with openPMD-api.
 * If not, see <http://www.gnu.org/licenses/>.
 */
#include <openPMD/openPMD.hpp>

#include <stdio.h> 
#include <unistd.h> 

#include <iostream>
#include <memory>
#include <cstddef>
#include <string>


using std::cout;
using namespace openPMD;


void getfullInfo(const std::string& fname){

  Series o = Series(fname, Access::READ_ONLY);

    std::cout << "Read iterations ";
    for( auto const& val : o.iterations )
        std::cout << '\t' << val.first;
    
    std::cout << "\n Read attributes in the root:\n";
    for( auto const& val : o.attributes() )
        std::cout << '\t' << val << '\n';
    std::cout << '\n';

    std::cout << "basePath - " << o.basePath() << '\n'
              << "iterationEncoding - " << o.iterationEncoding() << '\n'
              << "iterationFormat - " << o.iterationFormat() << '\n'
              << "meshesPath - " << o.meshesPath() << '\n'
              << "openPMD - " << o.openPMD() << '\n'
              << "openPMDextension - " << o.openPMDextension() << '\n'
      //  << "particlesPath - " << o.particlesPath() << '\n'
              << '\n';

    std::cout << "Read attributes in basePath:\n";
    for( auto const& a : o.iterations.attributes() )
        std::cout << '\t' << a << '\n';
    std::cout << '\n';

    std::cout << "Read iterations in basePath:\n";
    for( auto const& i : o.iterations )
        std::cout << '\t' << i.first << '\n';
    std::cout << '\n';

    for( auto const& i : o.iterations )
    {
        std::cout << "Read attributes in iteration " << i.first << ":\n";
        for( auto const& val : i.second.attributes() )
            std::cout << '\t' << val << '\n';
        std::cout << '\n';

        std::cout << i.first << ".time - " << i.second.time< float >() << '\n'
                  << i.first << ".dt - " << i.second.dt< float >() << '\n'
                  << i.first << ".timeUnitSI - " << i.second.timeUnitSI() << '\n'
                  << '\n';

        std::cout << "Read attributes in meshesPath in iteration " << i.first << ":\n";
        for( auto const& a : i.second.meshes.attributes() )
            std::cout << '\t' << a << '\n';
        std::cout << '\n';

        std::cout << "Read meshes in iteration " << i.first << ":\n";
        for( auto const& m : i.second.meshes )
            std::cout << '\t' << m.first << '\n';
        std::cout << '\n';

    }
}

int main(int argc, char *argv[])
{

    int opt;
    char* fname=NULL;
      
    while((opt = getopt(argc, argv, ":f:x")) != -1) 
    { 
        switch(opt) 
        { 
            case 'x': 
            case 'f': 
                printf("filename: %s\n", optarg);
		fname=optarg;
                break; 
        } 
    }   


    // Initiate the reading procedure

    getfullInfo(fname);
     

    std::cout << " Now create a serie " << std::endl;
    
    // Create a series and open the file
    Series series = Series( fname, Access::READ_ONLY);
    
    cout << "Read a Series with openPMD standard version "
         << series.openPMD() << '\n';
    cout << "The Series contains " << series.iterations.size() << " iterations:"; 

    // Loop over all iterations in the file
    int iter=0;
    for( auto const& i : series.iterations ){
      // Meshes  
      cout << "Iteration " << iter << " contains " << i.second.meshes.size() << " meshes:";
      for( auto const& m : i.second.meshes )
        cout << "\n\t" << m.first;
            
      cout << '\n';
      cout << "Iteration: "  <<  iter << "contains " << i.second.particles.size() << " particle species:";

      // Loop over species
      for( auto const& ps : i.second.particles ) {
        cout << "\n\t" << ps.first;
      }
      cout << '\n';	
      
      iter++;
      Iteration j = i.second;
    
      // Particles Species      
      openPMD::ParticleSpecies deuterons = j.particles["electron_l"];
      //std::shared_ptr<double> charge = electrons["charge"][openPMD::RecordComponent::SCALAR].loadChunk<double>();
      series.flush();
      
      // Access attribute within sub-group electron_gridx
      for( auto const& a : deuterons["momentum"].attributes() ){
	std::cout << '\t' << a << '\n';	
      }
      std::cout << '\n';
      
      auto p_x = deuterons["momentum"]["x"];
      Extent g_extent = p_x.getExtent();
      auto all_x  = p_x.loadChunk<double>();
      series.flush();

      cout << "Full Electron gridx starts with:\n\t{";
      for( size_t col = 0;  col < 5; ++col )
	cout << all_x.get()[col] << ", ";

      cout << "...}\n";
           
    }//!iteration++ 
    
    return 0;
}
