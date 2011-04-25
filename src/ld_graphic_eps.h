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
  This file contains the interfaces to the drawing code. Actual
  implementation in ld_graphic_eps.c
*/

/*
   ld_graphics_output_struct is the graphic state object
*/
typedef struct {
  FILE *fp;
  FILE *debug;
  int range;
  int offset;
  int headroom;
  int min_metric;
  int metric_length;
  int (*get_color)(geno_cptr res);
} ld_graphics_output_struct;

typedef ld_graphics_output_struct *ld_graphics_ptr;

/* initialize the graphic object, open outputfile, etc */
ld_graphics_ptr graphic_init(const char *filename, int ii, int jj, int i_depth, int scheme, int extra_headroom);

/* add_metric expands to accomodate the metric bar */
void graphic_add_metric(ld_graphics_ptr gh, int min_metric, int metric_length);

/*
   These two routines are no-OPs in EPS output.

   They are a historical relic from PPM output
   days which do need some routines for padding
   the beginning and ending of a scan line
*/
void graphic_scan_line_begin(ld_graphics_ptr gh, int idx_j);
void graphic_scan_line_end(ld_graphics_ptr gh, int idx_j);

/*
   Given a snp name and a position, do what's needed to be done
*/
void graphic_do_name(ld_graphics_ptr gh, int x, const char *name);

/* Given a snp listed order, and its chromosome position, draw the metric */
void graphic_do_metric(ld_graphics_ptr gh, int x, int trans_x);

/*
   Given a pair of positions, do what's needed to be done.
   calls get_geno_counts and draw pairs, etc.
   It is a composite routine which for drawing a paired-LD
   without generating an intermediate snp.dprime object
*/
void graphic_do_pair(ld_graphics_ptr gh, unsigned char *pos1, unsigned char *pos2, int i, int j, int rows, int do_notes);

/*
  Given a pair-wise stats, draw
 */
void graphic_draw_pair(ld_graphics_ptr gh, geno_cptr res, int i, int j, int do_notes);

/*
  Finalizing the graphic state and de-allocating resources, etc
*/
void graphic_close(ld_graphics_ptr gh);
