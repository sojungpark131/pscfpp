/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"

namespace Pscf
{ 

   using namespace Util;

   template class UnitCell<1>;
   template class UnitCell<2>;
   template class UnitCell<3>;

   /*
   * Extract a UnitCell<1>::LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, 
                              UnitCell<1>::LatticeSystem& lattice)
   {

      std::string buffer;
      in >> buffer;
      if (buffer == "Lamellar" || buffer == "lamellar") {
         lattice = UnitCell<1>::Lamellar;
      } else {
         UTIL_THROW("Invalid UnitCell<1>::LatticeSystem value input");
      }
      return in;
   }
   
   /* 
   * Insert a UnitCell<1>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, 
                            UnitCell<1>::LatticeSystem lattice) 
   {
      if (lattice == UnitCell<1>::Lamellar) {
         out << "lamellar";
      } else {
         UTIL_THROW("Invalid value of UnitCell<1>::Lamellar");
      } 
      return out;
   }

   /* 
   * Read the lattice system and set nParameter.
   */
   void UnitCell<1>::setNParameter()
   {
      if (lattice_ == UnitCell<1>::Lamellar) {
         nParameter_ = 1;
      } else {
         UTIL_THROW("Invalid lattice system value");
      } 
   }

   // Two-Dimensional Systems

   /*
   * Extract a UnitCell<2>::LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, 
                              UnitCell<2>::LatticeSystem& lattice)
   {

      std::string buffer;
      in >> buffer;
      if (buffer == "Square" || buffer == "square") {
         lattice = UnitCell<2>::Square;
      } else 
      if (buffer == "Rectangular" || buffer == "rectangular") {
         lattice = UnitCell<2>::Rectangular;
      } else 
      if (buffer == "Rhombic" || buffer == "rhombic") {
         lattice = UnitCell<2>::Rhombic;
      } else 
      if (buffer == "Hexagonal" || buffer == "hexagonal") {
         lattice = UnitCell<2>::Hexagonal; 
      } else 
      if (buffer == "Oblique" || buffer == "oblique") {
         lattice = UnitCell<2>::Oblique; 
      } else {
         #if 0
         Log::file() << "Unknown UnitCell<2>::LatticeSystem: " 
                     << buffer << std::endl;
         #endif
         UTIL_THROW("Invalid UnitCell<2>::LatticeSystem value input");
      }
      return in;
   }
   
   /* 
   * Insert a UnitCell<2>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator << (std::ostream& out, 
                              UnitCell<2>::LatticeSystem lattice) 
   {
      if (lattice == UnitCell<2>::Square) {
         out << "square";
      } else 
      if (lattice == UnitCell<2>::Rectangular) {
         out << "rectangular";
      } else
      if (lattice == UnitCell<2>::Rhombic) {
         out << "rhombic";
      } else
      if (lattice == UnitCell<2>::Hexagonal) {
         out << "hexagonal";
      } else
      if (lattice == UnitCell<2>::Oblique) {
         out << "oblique";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out;
   }

   /* 
   * Read the lattice system and set nParameter.
   */
   void UnitCell<2>::setNParameter()
   {
      if (lattice_ == UnitCell<2>::Square) {
         nParameter_ = 1;
      } else 
      if (lattice_ == UnitCell<2>::Rhombic) {
         nParameter_ = 1;
      } else
      if (lattice_ == UnitCell<2>::Hexagonal) {
         nParameter_ = 1;
      } else
      if (lattice_ == UnitCell<2>::Rectangular) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<2>::Oblique) {
         nParameter_ = 3;
      } else {
         UTIL_THROW("Invalid lattice system value");
      } 
   }

   // Three-Dimensional Systems

   /* 
   * Extract a UnitCell<3>::LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, 
                              UnitCell<3>::LatticeSystem& lattice)
   {

      std::string buffer;
      in >> buffer;
      if (buffer == "Cubic" || buffer == "cubic") {
         lattice = UnitCell<3>::Cubic;
      } else 
      if (buffer == "Tetragonal" || buffer == "tetragonal") {
         lattice = UnitCell<3>::Tetragonal;
      } else 
      if (buffer == "Orthorhombic" || buffer == "orthorhombic") {
         lattice = UnitCell<3>::Orthorhombic;
      } else 
      if (buffer == "Monoclinic" || buffer == "monoclinic") {
         lattice = UnitCell<3>::Monoclinic; 
      } else 
      if (buffer == "Triclinic" || buffer == "triclinic") {
         lattice = UnitCell<3>::Triclinic; 
      } else 
      if (buffer == "Rhombohedral" || buffer == "rhombohedral") {
         lattice = UnitCell<3>::Rhombohedral;
      } else 
      if (buffer == "Hexagonal" || buffer == "hexagonal") {
         lattice = UnitCell<3>::Hexagonal; 
      } else {
         #if 0
         Log::file() << "Unknown UnitCell<3>::LatticeSystem: " 
                     << buffer << std::endl;
         #endif
         UTIL_THROW("Invalid UnitCell<3>::LatticeSystem value input");
      }
      return in;
   }
   
   /* 
   * Insert a UnitCell<3>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, 
                            UnitCell<3>::LatticeSystem lattice) 
   {
      if (lattice == UnitCell<3>::Cubic) {
         out << "cubic";
      } else 
      if (lattice == UnitCell<3>::Tetragonal) {
         out << "tetragonal";
      } else
      if (lattice == UnitCell<3>::Orthorhombic) {
         out << "orthorhombic";
      } else
      if (lattice == UnitCell<3>::Monoclinic) {
         out << "monoclinic";
      } else
      if (lattice == UnitCell<3>::Triclinic) {
         out << "triclinic";
      } else
      if (lattice == UnitCell<3>::Rhombohedral) {
         out << "rhombohedral";
      } else
      if (lattice == UnitCell<3>::Hexagonal) {
         out << "hexagonal";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out;
   }

   /* 
   * Read the lattice system and set nParameter.
   */
   void UnitCell<3>::setNParameter()
   {
      if (lattice_ == UnitCell<3>::Cubic) {
         nParameter_ = 1;
      } else 
      if (lattice_ == UnitCell<3>::Tetragonal) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<3>::Orthorhombic) {
         nParameter_ = 3;
      } else
      if (lattice_ == UnitCell<3>::Monoclinic) {
         nParameter_ = 4;
      } else
      if (lattice_ == UnitCell<3>::Triclinic) {
         nParameter_ = 6;
      } else
      if (lattice_ == UnitCell<3>::Rhombohedral) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<3>::Hexagonal) {
         nParameter_ = 2;
      } else {
         UTIL_THROW("Invalid value");
      } 
   }

} 