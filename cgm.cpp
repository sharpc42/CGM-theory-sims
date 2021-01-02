//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file   implode.cpp
//  \author Christopher Sharp
//  \brief  Problem generator for a cooling CGM environment in the limit with clouds.  Works
//          as intended in Cartesian but will run in cylindrical and spherical coordinates:
//          without error or good results.

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <math.h>
#include <vector>
#include <typeinfo>
#include <queue>

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

Real rad;
Real pa;
Real da;
Real prat;
Real drat;
Real press_conv;
Real trat;
Real H;

//Real drt   = 1.0;
Real temp6 = 0.012;                 // k_B * T_6 = P_6 / rho_6 in code units
Real m_H   = 8.42e-58;					    // hydrogen mass in solar masses

std::vector<string> celldE;

/**
 * 
 * This is jointly the cooling and heating function for the evolution of the mesh.
 * Inspired by McCourt, et al. (2012), it takes the form of first calculating the
 * cooling. This is as a function of the cooling time in proportion to dynamical time.
 * Then the cooling is spatially averaged at each height, and this magnitude becomes
 * the heating. The positive energy change is the difference between the heating and
 * cooling, while the negative energy change divides the energy by 1 + | dE / E |
 * in order to prevent negative energies.
 * 
**/

void CoolHeat(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {

  Real g          = pmb->peos->GetGamma();
  Real gm1        = g - 1.0;

  Real temp5	    = 0.1 * temp6;                       // T = 10^5 K in code units
  Real temp7      = 10 * temp6;                        // T = 10^7 K in code units

  Real x_0     	  = (pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min) / 2;
  Real y_0        = (pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min) / 2;
  Real z_0        = (pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min) / 2;
  Real x;
  Real y;
  Real z;

  Real pres;
  Real dens;
  Real temp;
  Real t_dyn;
  Real cooling;
  Real heating;

  Real cool_sum;  // for the spatial average of cooling values
  int cell_sum;   // to divide for the spatial average

  queue<Real> cool;
  queue<Real> cell_cool;

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    z = pmb->pcoord->x3v(k);

    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
      cell_sum = 0;
      cool_sum = 0;

      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);

	      t_dyn    = std::pow(2 * (pmb->pcoord->x2v(j) + H) / g,0.5);
	      cooling = cons(IEN,k,j,i) / trat / t_dyn;

        cell_cool.push(cooling);
        cool_sum += cooling;
        cell_sum++;
      }

      heating = cool_sum / cell_sum;    // spatial average of cooling for this height

      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);  

        pres = cons(IEN,k,j,i);
        dens = cons(IDN,k,j,i);
        temp = gm1 * pres / dens;

        cooling = cell_cool.front();
        cell_cool.pop();
        Real dE = heating - cooling;   // needs to be absolute value
        
        // update the cell energy with dE

          if ((cons(IEN,k,j,i) + (dt * dE ) > 0) && (temp > temp5) && (temp < temp7)) {
            cons(IEN,k,j,i) += dt * dE;
          }
        //}

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
  
  trat = pin->GetReal("problem", "trat");
  H    = pin->GetReal("problem", "H");

  EnrollUserExplicitSourceFunction(CoolHeat);

  return;
}

//==========================================================================
//! \fn Real perturbModeSum()
//  \brief Real superposition of modes for perturbation
//==========================================================================

Real perturbModeSum(Real sum,int ne,int n,int x,int y,Real L) {
  Real mode = (1 / ((rand() % 10000) + 1)) * std::cos(n * 2 * M_PI * x / L);
  if (n == ne) {
    return mode;
  } else {
    return mode -= perturbModeSum(mode,ne,++n,x,y,L);
  }
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
    b0 = std::sqrt(pin->GetReal("problem","pbrat") * 2 * pa);
    angle = (PI/180.0)*pin->GetReal("problem","angle");
  }

  Real gamma = peos->GetGamma();
  Real gm1   = gamma - 1.0;
  Real t_0   = gm1 * pa / da;

  // get coordinates of center of blast and convert to Cartesian if needed
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
    x0 = x1_0 * std::sin(x2_0) * std::cos(x3_0);
    y0 = x1_0 * std::sin(x2_0) * std::sin(x3_0);
    z0 = x1_0 * std::cos(x2_0);
  } else {
    // check for an unrecognized coordinate system at top of code block
    std::stringstream msg;
    msg << "### FATAL ERROR in cgm.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // setup uniform ambient medium with spherical cold region
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real r;
    Real x;
    Real y;
    Real z;
    r = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    if (COORDINATE_SYSTEM == "cylindrical") {
      x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
      z = pcoord->x3v(k);
      r = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else if (COORDINATE_SYSTEM == "spherical_polar") {
      x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
      y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
      z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
      r = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
    } else {
      x = pcoord->x1v(i);
      y = pcoord->x2v(j);
      z = pcoord->x3v(k);
    }

    Real L   = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;

    Real Y      = y + H;   // hard-coded; change tp user parameter-based
    Real den    = da * std::exp(-0.1 * (std::pow(1 + (10 * Y / H) * (10 * Y / H),0.5) - 1));      
                          // hard-coded a = H / 10 as a check;
												  // eventually have have user pass as
												  // a parameter the ratio H / a (= 10 here)
    Real k_ptrb = 4 * M_PI / L;  // originally 8, now 2
    //Real dro    = 0.1 * den * std::cos(k_ptrb * x);
    //Real dro    = 0.1 * den * modeSum(0.0,40,2,x,L);
    Real dro    = den * perturbModeSum(0.0,40,2,x,y,L);

    Real pres = den * temp6 / gm1;
    Real temp = pres * gm1 / den;
    Real temp5_5    = 0.316228 * temp6;
    Real temp_0 = pa * gm1 / da;

    phydro->u(IDN,k,j,i) = den + dro;
    //phydro->u(IDN,k,j,i) = den;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;

    phydro->u(IEN,k,j,i) = pres * gm1;

    }

  }}

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {

    Real kpcPerCell = (pcoord->x1v(ie) - pcoord->x1v(is)) / (ie - is);
    Real k_vec      = 2 * PI / (2 * kpcPerCell);

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie+1; ++i) {
          if (COORDINATE_SYSTEM == "cartesian") {
            
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real L   = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;

            // currently: linear field instead of spatially periodic field;
            // switch commented-out to change
            
            //pfield->b.x1f(k,j,i) = b0 * std::cos(k_vec * pcoord->x1v(i)) * std::sin(k_vec * pcoord->x2v(j));
            pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
            
            // Only Cartesian supported for now for spatially periodic field
            
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

            Real y = pcoord->x2v(j);
            Real L   = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
            
            // currently: linear field instead of spatially periodic field,
            // switch commented-out to change
            
            //pfield->b.x2f(k,j,i) = -1 * b0 * std::sin(k_vec * pcoord->x1v(i)) * std::cos(k_vec * pcoord->x2v(j)); 
            pfield->b.x2f(k,j,i) = b0 * std::sin(angle);

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
    
    // z-component of field, and sets magnetic pressure/energy density

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
          phydro->u(IEN,k,j,i) += 0.5*b0*b0; // magnetic pressure
        }
      }
    }
  }
}
