#ifndef PSPC_POLYMER_TPP
#define PSPC_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pscf {
namespace Pspc { 

   template <int D>
   Polymer<D>::Polymer()
   {  setClassName("Polymer");}

   template <int D>
   Polymer<D>::~Polymer()
   {}

   template <int D>
   void Polymer<D>::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Species::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi;
   }

   template <int D>
   void Polymer<D>::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Species::Open);  
      mu_ = mu; 
   }

   /*
   * Set unit cell dimensions in all solvers.
   */ 
   template <int D>
   void Polymer<D>::setupUnitCell(UnitCell<D> const & unitCell)
   {
      // Set association to unitCell
      unitCellPtr_ = &unitCell;

      for (int j = 0; j < nBlock(); ++j) {
         block(j).setupUnitCell(unitCell);
      }
      std::cout<<"setup unit cell should only happens once, if you see this tell gk"<<std::endl;
   }

   /*
   * Compute solution to MDE and block concentrations.
   */ 
   template <int D>
   void Polymer<D>::compute(DArray<WField> const & wFields)
   {

      now = Timer::now();
      setupTimer.start(now);
      // Setup solvers for all blocks
      //not all blocks need solver
      int monomerId;
      int blockId;
      for(int i = 0; i < Base::firstMonomerSeenCount(); ++i) {
         if( i == 0) { //do all backbones
            for(int j = 0; j < Base::firstMonomerSeen(1); ++j) {
               monomerId = block(j).monomerId();
               block(j).setupSolver(wFields[monomerId]);
            }
         } else {
            blockId = Base::firstMonomerSeen(i);
            monomerId = block(blockId).monomerId();
            block(blockId).setupSolver(wFields[monomerId]);
         }
      }

      //for (int j = 0; j < nBlock(); ++j) {
      //   monomerId = block(j).monomerId();
      //   block(j).setupSolver(wFields[monomerId]);
      //}

      // Call base class PolymerTmpl solver() function
      //we implement this directly, a lot of assumption otherise
      //solve();

      // Clear all propagators
      for (int j = 0; j < minPropgCount(); ++j) {
         Base::propagator(j).setIsSolved(false);
      }

      now = Timer::now();
      setupTimer.stop(now);

      // Solve modified diffusion equation for all propagators
      // The brush polymer follows a specific order
      // 1) solve the forward propg for the side chain
      // 2-N-1) solve the forward and backward for the backbone
      // N) solve the summation of the backward propg of the side chain
      

      now = Timer::now();
      propagatorTimer.start(now);

      now = Timer::now();
      propTrueTimer.start(now);

      //determine how many forward propagators to solve
      for (int j = 0; j < minPropgCount() - (Base::firstMonomerSeen_.size() - 1); ++j) {
         UTIL_CHECK(Base::propagator(j).isReady());
         Base::propagator(j).solve();
         //std::cout<<j<< " is completed "<<std::endl;
      }
      
      now = Timer::now();
      propTrueTimer.stop(now);

      //get reference to head of propagator
      int nx = Base::propagator(minPropgCount() -1).meshSize();

      for (int j = 0; j < Base::firstMonomerSeen_.size() - 1; ++j) {

         int index = minPropgCount() - 1 - (Base::firstMonomerSeen_.size() - 2)+ j;
         //std::cout<<"index "<<index<<std::endl;
         WField& qh = Base::propagator(index).headFree();
         
         //for(int ix = 0;ix < nx; ++ix) {
         //   qh[ix] = 0;
         //}

         //for all blocks with the same type get the reverse direction
         int indexEnd;
         if(j + 2 >= Base::firstMonomerSeen_.size() ) {
            indexEnd = nBlock();
         } else {
            indexEnd = Base::firstMonomerSeen_[j + 2];
         }
         for(int blockId = Base::firstMonomerSeen_[j + 1]; 
             blockId < indexEnd; ++blockId) {

            now = Timer::now();
            propTrueTimer.start(now);

            block(blockId).propagator(1).computeHead();
            
            now = Timer::now();
            propTrueTimer.stop(now);

            //avoid adding the first propagator itself
            //qh is the first propagator
            if(blockId != Base::firstMonomerSeen_[j+1]) {
               qh += block(blockId).propagator(1).headFree();
            }
         }

         now = Timer::now();
         propTrueTimer.start(now);

         Base::propagator(index).solveFree();
         //Base::propagator(index).solve(qh);

         now = Timer::now();
         propTrueTimer.stop(now);


      }
      now = Timer::now();
      propagatorTimer.stop(now);


      now = Timer::now();
      rhoTimer.start(now);

      // Compute molecular partition function
      // Careful here -> ordering matters since not all propagator is solved
      // Still correct under the assumption of backbone specification
      double q = block(0).propagator(0).computeQ();

      //std::cout<<"big Q : "<<q<<std::endl;
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q);
      } 
      else if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q;
      }

      // Compute block concentration fields
      // Change backbone computation -> only once
      //minPropgcount is not exactly right here
      // should be backbone length + side chain amount
      double prefactor = phi_ / (q *length() );
      for (int i = 0; i < Base::firstMonomerSeen_[1]; ++i) {
         block(i).computeConcentration(prefactor);
      }
      for(int i = 1; i < Base::firstMonomerSeen_.size(); ++i) {
         block(Base::firstMonomerSeen_[i]).computeConcentration(prefactor);
      }

      now = Timer::now();
      rhoTimer.stop(now);


      double propTime = propagatorTimer.time();
      double setupTime = setupTimer.time();
      double rhoTime = rhoTimer.time();
      double propTrueTime = propTrueTimer.time();
      Log::file() << "setup time  = " << setupTime  << " s  \n";
      Log::file() << "prop time  = " << propTime  << " s  \n";
      Log::file() << "proptrue time  = " << propTrueTime  << " s  \n";
      Log::file() << "rho time  = " << rhoTime  << " s  \n";

   }

   /*
   * Compute stress from a polymer chain.
   */
   template <int D>
   void Polymer<D>::computeStress()
   {
     
      // Initialize all stress_ elements zero
      for (int i = 0; i < 6; ++i) {
        stress_[i] = 0.0;
      }

      // Compute and accumulate stress contributions from all blocks
      double prefactor = exp(mu_)/length();
      for (int i = 0; i < Base::firstMonomerSeen_[1]; ++i) {
         block(i).computeStress(prefactor);
         for (int j=0; j < unitCellPtr_->nParameter() ; ++j){
            stress_[j] += block(i).stress(j);
         }
      }
      for(int i = 1; i < Base::firstMonomerSeen_.size(); ++i) {
         block(Base::firstMonomerSeen_[i]).computeStress(prefactor);
         for (int j=0; j < unitCellPtr_->nParameter() ; ++j){
            stress_[j] += block(Base::firstMonomerSeen_[i]).stress(j);
         }

      }


   }

}
}
#endif
