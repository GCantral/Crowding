// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
//fix 
#include "fix_wall_depletion.h"
#include <cmath>
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWalldepletion::FixWalldepletion(LAMMPS *lmp, int narg, char **arg) :
  FixWall(lmp, narg, arg)
{
  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

void FixWalldepletion::precompute(int m)
{
  coeff1[m] = 48.0 * 1 * pow(sigma[m],12.0);
  //coeff2[m] = 3.0 * epsilon[m] * pow(sigma[m],3.0);
  coeff2[m] = 4.0*1 * pow(sigma[m],12.0);
  //coeff4[m] = epsilon[m] * pow(sigma[m],3.0);
  coeff3[m]=3.14159265358979323846*epsilon[m]/24.0;

  double rinv = 1.0/2.0;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset[m] = coeff2[m]*r4inv*r4inv*r4inv + coeff3[m]*16.0*sigma[m]*
              sigma[m] * sigma[m] ;
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWalldepletion::wall_particle(int m, int which, double coord)
{
  double delta,fwall,rinv,r2inv,r4inv,r13inv;
  double vn;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      if (delta<(sigma[m]/2.0)){
      
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r13inv = r4inv*r4inv*r4inv*rinv;
      fwall = side * (coeff1[m]*r13inv);
      f[i][dim] -= fwall;
      ewall[0] += coeff2[m]*r4inv*r4inv*r4inv - offset[m];
      ewall[m+1] += fwall;
      continue;
      }
      if (delta>(sigma[m]/2.0)){
      
      fwall = side * coeff2[m]*(18*sigma[m]*sigma[m] +
                                24*delta*delta-24*sigma[m]*delta);
      f[i][dim] -= fwall;
      ewall[0] += (-coeff2[m]*(3*sigma[m]-2*delta)*(3*sigma[m]-2*delta)*
                  (3*sigma[m]+2*delta));
      ewall[m+1] += fwall;
      }
    }


  if (onflag) error->one(FLERR,"Particle on or inside fix wall surface");
}
