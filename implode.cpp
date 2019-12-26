//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file   implode.cpp
//  \author Christopher Sharp 
//  \brief  Problem generator for spherical imploding cold CGM cloud problem.  Works 
//          as intended in Cartesian but will run in cylindrical and spherical
//          coordinates.

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <math.h>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../bvals/bvals.hpp"
#include "../utils/utils.hpp"

static Real grav;

Real rad;
Real pa;
Real da;
Real prat;
Real drat;
Real press_conv;

void CoolingFxn(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {

  std::cout<<"I at least made it to the cooling fxn!!!\n\n";

  Real g          = pmb->peos->GetGamma();
  Real gm1        = g - 1.0;

  Real temp6	  = 8.6136e-3;                                  // in k_B * T_6 = P_6 / rho_6 units converted to code units
  Real temp5	  = 0.1 * temp6;                                // T = 10^5 K as fraction of T = 10^6 K in code units
  Real temp5_5    = 0.316228 * temp6;                           // T = 10^5.5 K as fraction of T = 10^6 K in code units

  /**
  Real x_0        = std::abs(pmb->ie - pmb->is) / 2;
  Real y_0        = std::abs(pmb->je - pmb->js) / 2;
  Real z_0        = std::abs(pmb->ke - pmb->ks) / 2;
  **/

  Real x_0     	  = (pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min) / 2;
  Real y_0        = (pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min) / 2;
  Real z_0        = (pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min) / 2;

  Real x;
  Real y;
  Real z;  
  Real r;
  Real rad = std::pow((x_0 * x_0 + y_0 * y_0 + z_0 * z_0),0.5) / 5;  // def. a bit arbitrary but will help resolve spatial scale in here;

  // assuming some serious symmetry here for practical purposes
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    //z = -(k - z_0);
    z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      //y = -(j - y_0);
      y = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        //x   = -(i - x_0);
        x   = pmb->pcoord->x3v(i);
        r   = std::pow((x * x + y * y + z * z),0.5);

        std::cout << "x y z r rad " << x << " " << y << " " << z << " " << r << " " << rad << " " << "\n\n";

        Real pres = prim(IEN,k,j,i);
        Real dens = prim(IDN,k,j,i);
        Real temp = gm1 * pres / dens;

        Real t_cool = 250 * (pa / pres) * pow(temp / temp5_5,2.7);  // in Myr code units assuming metallicity ~ 0.3 (so correct to order-1)
        Real cooling = ((dt / t_cool) * cons(IEN,k,j,i) 
                        * std::exp(-(temp6 / temp) - (r / rad)));  // cooling that should die off exponentially with lower temp and higher distance

        std::cout << "coolratio " << cooling / cons(IEN,k,j,i) << " " << "\n\n";
        cons(IEN,k,j,i) -= cooling;

        // std::cout << "presureratio tempratio tcool " <<  pa / pres << " " << temp / temp5_5 << " " << t_cool << " " << "\n\n";

      }
    }
  }

  return;
}

//==========================================================================
//! \fn void Mesh::InitUserMeshData()
//  \brief Enroll cooling function
//==========================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {

  EnrollUserExplicitSourceFunction(CoolingFxn);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Cold underpressurized CGM cloud implosion problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real dens_conv  = 1.478e28;                                // in solar masses per cubic kpc  
  Real vel_conv   = 1.022e-6;                                // in kpc per Myr
  Real press_conv = 1.543e16;                                // in solar masses per kpc*squared Myr

  rad  = pin->GetReal("problem","radius");
  pa   = (pin->GetOrAddReal("problem","pamb",1.0)
              * press_conv);
  da   = (pin->GetOrAddReal("problem","damb",1.0)
              * dens_conv);
  prat = pin->GetReal("problem","prat");                // may don't need if using cooling
  Real drat = pin->GetOrAddReal("problem","drat",1.0);       // probably don't need altogether

  Real b0,angle;

  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = std::sqrt(pin->GetReal("problem","b0") * 2 * pa);
    angle = (PI/180.0)*pin->GetReal("problem","angle");
  }

  Real gamma = peos->GetGamma();
  Real gm1   = gamma - 1.0;
  
  grav = phydro->hsrc.GetG2();

  // get coordinates of center of blast and converts to Cartesian if needed
  Real x1_0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Real x2_0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Real x3_0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Real x0,y0,z0;
  if (COORDINATE_SYSTEM == "cartesian") {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  
  // setup uniform ambient medium with spherical cold region
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real r;
    if (COORDINATE_SYSTEM == "cartesian") {
      Real x = pcoord->x1v(i);
      Real y = pcoord->x2v(j);
      Real z = pcoord->x3v(k);
      r = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else if (COORDINATE_SYSTEM == "cylindrical") {
      Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
      Real z = pcoord->x3v(k);
      r = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else { // if (COORDINATE_SYSTEM == "spherical_polar")
      Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
      Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
      Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      r = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    }

    Real den = da;

    phydro->u(IDN,k,j,i) = den;    
    phydro->u(IM1,k,j,i) = 0.0 * den;
    phydro->u(IM2,k,j,i) = 0.0 * den;
    phydro->u(IM3,k,j,i) = 0.0 * den;
        
    // This sets gravitational energy for 2D and 3D cases seperately.
    // Gravitational acceleration should be passed as a user parameter
    // with units in accordance with kpc Myr^-2 (for us: 10^-5 kpc Myr^-2).
    // If desiring a gravitational source located "below" the simulated space:
    
    if (block_size.nx3 == 1) {
	  phydro->u(IEN,k,j,i) += grav * den * std::abs(pcoord->x2v(j) - 0.5 * (pcoord->x2v(js) - pcoord->x2v(je)));
	} else {
		phydro->u(IEN,k,j,i) += grav*den*(pcoord->x3v(k));
	}

    if (r < rad) {
      std::cout<<"Inside the appropriate area\n\n";
      CoolingFxn;
    }

    // shouldn't need below if using a cooling function
    
    if (NON_BAROTROPIC_EOS) {
      Real pres = pa;
      //if (rad < rin) {                     // originally less than rout
        //if (rad < rout) {                  // originally less than rin
          //pres = pa / prat;
        //} else {                           // add smooth ramp in pressure
          //Real f = (rad-rin) / (rout-rin);
         
        // original - as we enter the boundary, over-pressure dominates;
        // as we traverse through it, ambient pressure begins to dominate
         
        
        //Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
        
         
        // User - as we enter the boundary, ambient pressure dominates;
        // as we traverse through it, over-pressure begins to dominate
        // Note: doesn't seem to matter for now which is pa or prat*pa,
        // since there is no transition region at present, but this log
        // ramp does need to be included for good physical results
          //Real log_pres = f * std::log(prat*pa) + (1.0-f) * std::log(pa);
         
          //pres = std::exp(log_pres);
        }
      }

      // phydro->u(IEN,k,j,i) += pres/gm1;
      
//      if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
//        phydro->u(IEN,k,j,i) += den;
     
//    }

  }}

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {

    Real kpcPerCell = (pcoord->x1v(ie) - pcoord->x1v(is)) / (ie - is);
    Real k_vec      = 2 * PI / (20 * kpcPerCell);

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie+1; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {

            // long-term: consider adding function/object to handle
            // various magnetic field topologies independently of this
            // particular problem generator, according to user input
            // then this particular pgen could simply check divergence
            // req, but then is that already handled by Athena itself?
            
            // currently: spatially periodic field
            pfield->b.x1f(k,j,i) = b0 * std::cos(k_vec * pcoord->x1v(i)) * std::sin(k_vec * pcoord->x2v(j));
            
            // only Cartesian supported for now; some way to handle general
            // conversion in the future?
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            Real phi = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) =
                b0 * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0 * std::abs(std::sin(theta))
                * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {

            // again, may go with imported function/object for generic magnetic fields; here, spatially periodic
            pfield->b.x2f(k,j,i) = -1 * b0 * std::sin(k_vec * pcoord->x1v(i)) * std::cos(k_vec * pcoord->x2v(j)); 

            // again, only Cartesian supported at this time            
          } else if (COORDINATE_SYSTEM == "cylindrical") {
            Real phi = pcoord->x2v(j);
            pfield->b.x2f(k,j,i) =
              b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x2f(k,j,i) = b0 * std::cos(theta)
              * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
            if (std::sin(theta) < 0.0)
              pfield->b.x2f(k,j,i) *= -1.0;
          }
        }
      }
    }
    // this is interesting -- if we pursue a systematic approach where we iterate through
    // topologies of *vector potentials* then take the curl of the magnetic field (in which
    // case also a divergenceless check is not needed) then we face the possibility of 
    // having non-zero z components organically: but sometimes in the 2D plane. How to handle
    // correctly in the general case?
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          if (COORDINATE_SYSTEM == "cartesian" || COORDINATE_SYSTEM == "cylindrical") {
            pfield->b.x3f(k,j,i) = 0.0;
          } else { //if (COORDINATE_SYSTEM == "spherical_polar") {
            Real phi = pcoord->x3v(k);
            pfield->b.x3f(k,j,i) =
              b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0; // set magnetic pressure
        }
      }
    }
  }
}

/**

//==========================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Use cooling function to update cells
//==========================================================================

void MeshBlock::UserWorkInLoop() {

  int x,y,z;
  Real r;

  for (int k = ks; k <= ke; ++k) {
    z = pcoord->x3v(k);
    for (int j = js; j <= je; ++j) {
      y = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        x = pcoord->x1v(i);
        r = std::sqrt(SQR(x) + SQR(y));
        if (r <= rad) {
          //CoolingFxn(pmb,time,dt,phydro->w,pfield->b,phydro->u);
          CoolingFxn;
          std::cout<<"I am cooling!\n\n";
        } else {
          std::cout<<"I'm not cooling\n\n";
        }
      }
    }
  }

  return;
}

**/

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;
  AthenaArray<Real> pr;
  pr.DeleteAthenaArray();
  return;
}
