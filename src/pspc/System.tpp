#ifndef PSPC_SYSTEM_TPP
#define PSPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#if 0
#include <pspc/sweep/Sweep.h>
#include <pspc/sweep/SweepFactory.h>
#endif

#include <pspc/iterator/AmIterator.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/inter/Interaction.h>
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   System<D>::System()
    : mixture_(),
      unitCell_(),
      mesh_(),
      fft_(),
      groupName_(),
      basis_(),
      fileMaster_(),
      fieldIo_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      // sweepPtr_(0),
      // sweepFactoryPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      hasUnitCell_(false),
      isAllocated_(false),
      hasWFields_(false),
      hasCFields_(false)
      //hasSweep_(false)
   {  
      setClassName("System"); 

      fieldIo_.associate(unitCell_, mesh_, fft_, groupName_,
                         basis_, fileMaster_);

      interactionPtr_ = new ChiInteraction(); 
      iteratorPtr_ = new AmIterator<D>(this); 

      // sweepFactoryPtr_ = new SweepFactory(*this);
   }

   /*
   * Destructor.
   */
   template <int D>
   System<D>::~System()
   {
      if (interactionPtr_) {
         delete interactionPtr_;
      }
      if (iteratorPtr_) {
         delete iteratorPtr_;
      }
   }

   /*
   * Process command line options.
   */
   template <int D>
   void System<D>::setOptions(int argc, char **argv)
   {
      bool eflag = false;  // echo
      bool pFlag = false;  // param file 
      bool cFlag = false;  // command file 
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:f")) != -1) {
         switch (c) {
         case 'e':
            eflag = true;
            break;
         case 'p': // parameter file
            pFlag = true;
            pArg  = optarg;
            break;
         case 'c': // command file
            cFlag = true;
            cArg  = optarg;
            break;
         case 'i': // input prefix
            iFlag = true;
            iArg  = optarg;
            break;
         case 'o': // output prefix
            iFlag = true;
            oArg  = optarg;
            break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }
   
      // Set flag to echo parameters as they are read.
      if (eflag) {
         Util::ParamComponent::setEcho(true);
      }

      // If option -p, set parameter file name
      if (pFlag) {
         fileMaster().setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster().setCommandFileName(std::string(cArg));
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }

   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {
      readParamComposite(in, mixture());
      hasMixture_ = true;

      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      //int ns = mixture().nSolvent(); 
      int ns = 0;

      // Initialize homogeneous object
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      read(in, "unitCell", unitCell_);
      hasUnitCell_ = true;
      
      read(in, "mesh", mesh_);
      hasMesh_ = true;

      read(in, "groupName", groupName_);

      mixture().setMesh(mesh());
      mixture().setupUnitCell(unitCell());
      basis().makeBasis(mesh(), unitCell(), groupName_);

      allocate();
      isAllocated_ = true;

      // Initialize iterator
      readParamComposite(in, iterator());
      iterator().allocate();

      #if 0
      // Optionally instantiate a Sweep object
      readOptional<bool>(in, "hasSweep", hasSweep_);
      if (hasSweep_) {
         std::string className;
         bool isEnd;
         sweepPtr_ = 
            sweepFactoryPtr_->readObject(in, *this, className, isEnd);
         if (!sweepPtr_) {
            UTIL_THROW("Unrecognized Sweep subclass name");
         }
         sweepPtr_->setSystem(*this);
      }
      #endif
   }

   /*
   * Read default parameter file.
   */
   template <int D>
   void System<D>::readParam(std::istream& in)
   {
      readBegin(in, className().c_str());
      readParameters(in);  
      readEnd(in);
   }

   /*
   * Read default parameter file.
   */
   template <int D>
   void System<D>::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::allocate()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(hasMesh_);

      // Allocate wFields and cFields
      int nMonomer = mixture().nMonomer();
      wFields_.allocate(nMonomer);
      wFieldsRGrid_.allocate(nMonomer);
      wFieldsKGrid_.allocate(nMonomer);

      cFields_.allocate(nMonomer);
      cFieldsRGrid_.allocate(nMonomer);
      c2FieldsRGrid_.allocate(nMonomer*2);
      cFieldsKGrid_.allocate(nMonomer);
      
      for (int i = 0; i < nMonomer; ++i) {
         wField(i).allocate(basis().nStar());
         wFieldRGrid(i).allocate(mesh().dimensions());
         wFieldKGrid(i).allocate(mesh().dimensions());

         cField(i).allocate(basis().nStar());
         cFieldRGrid(i).allocate(mesh().dimensions());
         cFieldKGrid(i).allocate(mesh().dimensions());
      }
      
      for (int i = 0; i < nMonomer*2; ++i) {
         c2FieldRGrid(i).allocate(mesh().dimensions());
      }            
      
      isAllocated_ = true;
   }

   
   /*
   * Read and execute commands from a specified command file.
   */
   template <int D>
   void System<D>::readCommands(std::istream &in) 
   {
      UTIL_CHECK(isAllocated_);

      std::string command;
      std::string filename;

      bool readNext = true;
      while (readNext) {

         in >> command;
         Log::file() << command <<std::endl;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_W_BASIS") {
            in >> filename;
            Log::file() << " " << Str(filename, 20) <<std::endl;
            fieldIo().readFieldsBasis(filename, wFields());
            fieldIo().convertBasisToRGrid(wFields(), wFieldsRGrid());
            hasWFields_ = true;
            hasCFields_ = false;
         } else
         if (command == "READ_W_RGRID") {
            in >> filename;
            Log::file() << " " << Str(filename, 20) <<std::endl;
            fieldIo().readFieldsRGrid(filename, wFieldsRGrid());
            fieldIo().convertRGridToBasis(wFieldsRGrid(), wFields());
            hasWFields_ = true;
            hasCFields_ = false;
         } else
         if (command == "ITERATE") {

            Log::file() << std::endl;
            Log::file() << std::endl;

            // Read w (chemical potential fields) if not done previously 
            if (!hasWFields_) {
               in >> filename;
               Log::file() << " " << Str(filename, 20) <<std::endl;
               fieldIo().readFieldsBasis(filename, wFields());
               fieldIo().convertBasisToRGrid(wFields(), wFieldsRGrid());
               hasWFields_ = true;
            }

            // Iterative solution of SCF equations
            int fail = iterator().solve();
            hasCFields_ = true;

            if (fail) {
               Log::file() << "Iterator failed to converge\n";
            } else {
               computeFreeEnergy();
               outputThermo(Log::file());
            }

         } else
         if (command == "SOLVE_MDE") {
            Log::file() << std::endl;
            Log::file() << std::endl;
           
            // Read w (chemical potential fields) if not done previously 
            if (!hasWFields_) {
               in >> filename;
               Log::file() << " " << Str(filename, 20) <<std::endl;
               fieldIo().readFieldsBasis(filename, wFields());
               fieldIo().convertBasisToRGrid(wFields(), wFieldsRGrid());
               hasWFields_ = true;
            }

            // Solve the modified diffusion equation (without iteration)
            mixture().compute(wFieldsRGrid(), cFieldsRGrid(), c2FieldsRGrid());

            // Convert c fields from r-grid to basis
            fieldIo().convertRGridToBasis(cFieldsRGrid(), cFields());
            hasCFields_ = true;

         } else
         if (command == "WRITE_W_BASIS") {
            UTIL_CHECK(hasWFields_);
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsBasis(filename, wFields());
         } else 
         if (command == "WRITE_W_RGRID") {
            UTIL_CHECK(hasWFields_);
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsRGrid(filename, wFieldsRGrid());
         } else 
         if (command == "WRITE_C_BASIS") {
            UTIL_CHECK(hasCFields_);
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsBasis(filename, cFields());
         } else
         if (command == "WRITE_C_RGRID") {
            UTIL_CHECK(hasCFields_);
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsRGrid(filename, cFieldsRGrid());
            fieldIo().writeFieldsRGrid("rho_segment", c2FieldsRGrid());
         } else
         if (command == "BASIS_TO_RGRID") {

            // Note: This and other conversions use the c-field arrays
            // for storage, and thus corrupt previously stored values.
            hasCFields_ = false;

            // Read in basis format
            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;
            fieldIo().readFieldsBasis(inFileName, cFields());

            // Convert
            fieldIo().convertBasisToRGrid(cFields(), cFieldsRGrid());

            // Write in r-grid format
            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;
            fieldIo().writeFieldsRGrid(outFileName, cFieldsRGrid());

         } else 
         if (command == "RGRID_TO_BASIS") {
            hasCFields_ = false;

            // Read in r-grid format
            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;
            fieldIo().readFieldsRGrid(inFileName, cFieldsRGrid());

            // Convert from r-grid to basis
            fieldIo().convertRGridToBasis(cFieldsRGrid(), cFields());

            // Write in basis format
            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;
            fieldIo().writeFieldsBasis(outFileName, cFields());

         } else
         if (command == "KGRID_TO_RGRID") {
            hasCFields_ = false;

            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;


            fieldIo().readFieldsKGrid(inFileName, cFieldsKGrid());
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               fft().inverseTransform(cFieldKGrid(i), cFieldRGrid(i));
            }

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;
            fieldIo().writeFieldsRGrid(outFileName, cFieldsRGrid());

         } else
         if (command == "RGRID_TO_KGRID") {
            hasCFields_ = false;

            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;

            fieldIo().readFieldsRGrid(inFileName, cFieldsRGrid());
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               fft().forwardTransform(cFieldRGrid(i), cFieldKGrid(i));
            }

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;
            fieldIo().writeFieldsKGrid(outFileName, cFieldsKGrid());

         } else
         if (command == "RHO_TO_OMEGA") {
            hasCFields_ = false;
            hasWFields_ = false;

            // Open input file
            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;
            std::ifstream inFile;
            fileMaster().openInputFile(inFileName, inFile);

            // Read concentration field 
            std::string label;
            int nStar,nM;
            inFile >> label;
            inFile >> nStar;
            inFile >> label;
            inFile >> nM;
            int idum;
            for (int i = 0; i < nStar; ++i) {
               inFile >> idum;
               for (int j = 0; j < nM; ++j) {
                  inFile >> cField(j)[i];
               }
            }
            inFile.close();

            // Compute w fields from c fields
            for (int i = 0; i < basis().nStar(); ++i) {
               for (int j = 0; j < mixture().nMonomer(); ++j) {
                  wField(j)[i] = 0;
                  for (int k = 0; k < mixture().nMonomer(); ++k) {
                     wField(j)[i] += interaction().chi(j,k) * cField(k)[i];
                  }
               }
            }

            // Convert to r-grid format
            fieldIo().convertBasisToRGrid(wFields(), wFieldsRGrid());
            hasWFields_ = true;
            hasCFields_ = false;

            // Write w field
            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) << std::endl;
            fieldIo().writeFieldsBasis(outFileName, wFields());
            hasWFields_ = true;

         } else
         if (command == "OUTPUT_STARS") {

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) << std::endl;

            std::ofstream outFile;
            fileMaster().openOutputFile(outFileName, outFile);
            fieldIo().writeFieldHeader(outFile, mixture().nMonomer());
            basis().outputStars(outFile);

         } else
         if (command == "OUTPUT_WAVES") {

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) << std::endl;

            std::ofstream outFile;
            fileMaster().openOutputFile(outFileName, outFile);
            fieldIo().writeFieldHeader(outFile, mixture().nMonomer());
            basis().outputWaves(outFile);

         } else {
            Log::file() << "Error: Unknown command  " 
                        << command << std::endl;
            readNext = false;
         }
      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   template <int D>
   void System<D>::readCommands()
   {  
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile()); 
   }

  
   /*
   * Compute Helmoltz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      fHelmholtz_ = 0.0;
 
      // Compute ideal gas contributions to fHelhmoltz_
      Polymer<D>* polymerPtr;
      double phi, mu, length;
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         polymerPtr = &mixture().polymer(i);
         phi = polymerPtr->phi();
         mu = polymerPtr->mu();
         // Recall: mu = ln(phi/q)
         length = polymerPtr->length();
         fHelmholtz_ += phi*( mu - 1.0 )/length;
      }

      int nm  = mixture().nMonomer();
      int nStar = basis().nStar();
      double temp = 0;

      for (int i = 0; i < nm; ++i) {
         
         for (int j = i + 1; j < nm; ++j) {
            for (int k = 0; k < nStar; ++k) {
               fHelmholtz_+=
                  cFields_[i][k] * interaction().chi(i,j) * cFields_[j][k];
            }
         }

         for (int j = 0; j < nStar; ++j) {
            temp += wFields_[i][j] * cFields_[i][j];
         }

      }
      fHelmholtz_ -= temp;

      // Compute pressure
      pressure_ = -fHelmholtz_;
      for (int i = 0; i < np; ++i) {
         polymerPtr = &mixture().polymer(i);
         phi = polymerPtr->phi();
         mu = polymerPtr->mu();
         length = polymerPtr->length();

         pressure_ += mu * phi /length;
      }

   }

   template <int D>
   void System<D>::outputThermo(std::ostream& out)
   {
      out << std::endl;
      out << "fHelmholtz = " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure   = " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;

      out << "Polymers:" << std::endl;
      out << "    i"
          << "        phi[i]      "
          << "        mu[i]       " 
          << std::endl;
      for (int i = 0; i < mixture().nPolymer(); ++i) {
         out << Int(i, 5) 
             << "  " << Dbl(mixture().polymer(i).phi(),18, 11)
             << "  " << Dbl(mixture().polymer(i).mu(), 18, 11)  
             << std::endl;
      }
      out << std::endl;
   }

   template <int D>
   void System<D>::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      //int ns = mixture().nSolvent(); 
      int ns = 0;
      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      // Allocate c_ work array, if necessary
      if (c_.isAllocated()) {
         UTIL_CHECK(c_.capacity() == nm);
      } else {
         c_.allocate(nm);
      }

      int i;   // molecule index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      for (i = 0; i < np; ++i) {

         // Initial array of clump sizes 
         for (j = 0; j < nm; ++j) {
            c_[j] = 0.0;
         }

         // Compute clump sizes for all monomer types.
         nb = mixture().polymer(i).nBlock(); 
         for (k = 0; k < nb; ++k) {
            Block<D>& block = mixture().polymer(i).block(k);
            j = block.monomerId();
            c_[j] += block.length();
         }
 
         // Count the number of clumps of nonzero size
         nc = 0;
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 1.0E-8) {
               ++nc;
            }
         }
         homogeneous_.molecule(i).setNClump(nc);
 
         // Set clump properties for this Homogeneous::Molecule
         k = 0; // Clump index
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 1.0E-8) {
               homogeneous_.molecule(i).clump(k).setMonomerId(j);
               homogeneous_.molecule(i).clump(k).setSize(c_[j]);
               ++k;
            }
         }
         homogeneous_.molecule(i).computeSize();

      }

   }

} // namespace Pspc
} // namespace Pscf
#endif
