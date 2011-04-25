/*
 *  Copyright (C) 2006  Hin-Tak Leung
 *
 *  Some addition by David Clayton 22/01/2009 to add covariance and odds ratios
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 51 Franklin Street
 *  Fifth Floor, Boston, MA 02110-1301  USA.
 */

/*
   see accompanying notes which describes the methodology
*/

/*
  the geno_count is a mixture of inputs, intermediate
  outputs and final outputs; in a future tidy-up
  one might want to separate the three
*/
typedef struct {
  int count[9];
  double *expt;
  int total;
  int a;
  int b;
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int a0;
  int a1;
  int a2;
  int a3;
  double p;
  double u;
  double v;
  double w;
  double x;
  double bigD;
  double rsq2;
  double dprime;
  double lod;
  double odds_ratio;
  double covar;
  /* */
  double post_p;
  signed char sign_of_r; /* sign_of_r is really just sign of bigD, but eventually
			 we might want to separate/hide all the internal values
		      */
} geno_count;

typedef geno_count *geno_cptr;

/*
  given two snp.matrix position and length, read along the vector
  and populates the initial inputs in a new geno_count object
  and returns it.

  caller is responsible for free()
 */
geno_cptr get_geno_count(unsigned char *pos1, unsigned char *pos2, int length);
