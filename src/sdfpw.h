/*
 *  Copyright (C) 2006  Hin-Tak Leung
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
  This is the external R<->C interface of the pairwise LD calculation
  code. See sdfpw.c for implementations.
*/

/*
  calculates the pair-wise LD of two specific snps;
  return nothing. All output to console
*/
SEXP snp_pairwise(SEXP v, SEXP i, SEXP j);

/*
  draw graphics directly without generating an
  intermediate snp.dprime object
*/
SEXP snp_pair_graphics(SEXP v, SEXP fileoutput, SEXP i, SEXP j, SEXP depth, SEXP do_notes);

/*
  given a snp.dprime object, draw the graphics

  scheme: colour scheme - currently only two are supported,
  0 for the haploview dfault, 1 for r^2 grayscale  
*/
SEXP snp_dprime_draw(SEXP list_in, SEXP fileoutput, SEXP scheme, SEXP do_notes, SEXP metric);

/*
  calculates the pairwise-LD and returns as a snp.dprime object
 */
SEXP snp_pair_range(SEXP v, SEXP i, SEXP j, SEXP depth, SEXP signed_r);
