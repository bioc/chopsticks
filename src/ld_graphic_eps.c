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
   See ld_graphic_eps.h for description for external interfaces
   to routines here. Only internal routines are documented here.
*/

#include <stdio.h>  /* file operations */
#include <stdlib.h> /* calloc */
#include <string.h> /* memset */
#include "pairwise_linkage.h"
#include "ld_graphic_eps.h"
#include "ld_graphic_eps_priv.h"

/* internal routines */
/* 
   These are respondible for the color-schemes.
   For each new color-scheme, we add one more
   here of the exact same prototype, and hook up the "scheme" switch
   in graphic_init() to select between them.
*/
static int get_color_std(geno_cptr res);
static int get_color_rsq(geno_cptr res);

ld_graphics_ptr graphic_init(const char *filename, int ii, int jj, int i_depth, int scheme, int extra_headroom)
{
  int junk; /* useless to silent the compiler */
  ld_graphics_ptr gh = (ld_graphics_ptr)calloc(1, sizeof(ld_graphics_output_struct));
  double scale = 1.0;

  if (!gh)
    return NULL;

  gh->fp = fopen(filename, "w");
  gh->debug = NULL;
  gh->headroom = 50;
  gh->min_metric = 0;
  gh->metric_length = 0;
  /* gh->debug = fopen("/tmp/debug", "w"); */
  gh->range = jj - ii + 1;
  gh->offset = ii - 1;
  switch(scheme)
    {
    case 0:
      gh->get_color = get_color_std;
      break;
    case 1:
      gh->get_color = get_color_rsq;
      break;
    default:
      gh->get_color = get_color_std;
      break;
    }

  if(!gh->fp) {
    return NULL;
  }
  if (extra_headroom) {
    gh->headroom = 150;
  }
  
  if (gh->range) 
    {
      scale = 70.0 / gh->range; /* 70 ~ 842/12 */
    }

  junk = fwrite(EPS_HEADER, strlen(EPS_HEADER), 1, gh->fp);
  fprintf(gh->fp, EPS_BOUNDBOX(gh->range , i_depth, scale, gh->headroom));
  junk = fwrite(EPS_PROLOG_1, strlen(EPS_PROLOG_1), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_2, strlen(EPS_PROLOG_2), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_3, strlen(EPS_PROLOG_3), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_4, strlen(EPS_PROLOG_4), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_5, strlen(EPS_PROLOG_5), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_6, strlen(EPS_PROLOG_6), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_7, strlen(EPS_PROLOG_7), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_8, strlen(EPS_PROLOG_8), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_9, strlen(EPS_PROLOG_9), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_10, strlen(EPS_PROLOG_10), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_11, strlen(EPS_PROLOG_11), 1, gh->fp);
  junk = fwrite(EPS_PROLOG_12, strlen(EPS_PROLOG_12), 1, gh->fp);
  fprintf(gh->fp, EPS_SCALE_TO_FIT(scale));
  fprintf(gh->fp, EPS_CLIP_PATH(gh->range , i_depth, gh->headroom));
  junk = fwrite(EPS_BODY, strlen(EPS_BODY), 1, gh->fp);
  if(junk != strlen(EPS_BODY)) return NULL;
  fprintf(gh->fp, EPS_TRANSLATE(0, i_depth));
    return gh;
}

void graphic_add_metric(ld_graphics_ptr gh, int min_metric, int metric_length){
  gh->min_metric = min_metric;
  gh->metric_length = metric_length;
}

/* no-OPs */
void graphic_scan_line_begin(ld_graphics_ptr gh, int idx_j) 
{
  return;
}

/* no-OPs */
void graphic_scan_line_end(ld_graphics_ptr gh, int idx_j)
{
  return;
}

int graphic_close(ld_graphics_ptr gh)
{
  int junk;
  junk = fwrite(EPS_ENDING, strlen(EPS_ENDING), 1, gh->fp);
  if (junk != strlen(EPS_ENDING))
    return -1;
  fclose(gh->fp);
  if (gh->debug !=NULL) 
    fclose(gh->debug);
  free(gh);
  return 0;
}

void graphic_do_name(ld_graphics_ptr gh, int x, const char *name) 
{
  fprintf(gh->fp, EPS_TEXT(x-gh->offset, name));
}

void graphic_do_metric(ld_graphics_ptr gh, int x, int trans_x) 
{
  if(gh->metric_length) {
    fprintf(gh->fp, EPS_METRIC(x-gh->offset, ((float)(trans_x - gh->min_metric) * (gh->range -1))/gh->metric_length));
  }
}

void graphic_do_pair(ld_graphics_ptr gh, unsigned char *pos1, unsigned char *pos2, int i, int j, int rows, int do_notes)
{
  geno_cptr res = NULL;
  res = get_geno_count(pos1, pos2, rows);
  graphic_draw_pair(gh, res, i, j, do_notes);
  free(res->expt);
  free(res);
}

void graphic_draw_pair(ld_graphics_ptr gh, geno_cptr res, int i, int j, int do_notes)
{
  int color = 0;
  color = gh->get_color(res);
  fprintf(gh->fp, EPS_RHOMBUS(i-gh->offset, j, color));
  if (do_notes)
    fprintf(gh->fp, EPS_NOTE(i-gh->offset,j,res->dprime, res->rsq2, res->lod, color));
  if(gh->debug != NULL)
    fprintf(gh->debug, "%d\t%d\t%f\t%f\t%f\n", i, i + j + 1, res->dprime, res->lod, res->rsq2);
}

/* Haploview's color decision code in DPrimeDisplay.java */
/*
  double d = thisPair.getDPrime();
  double l = thisPair.getLOD();
  Color boxColor = null;
  if (l > 2) {
    if (d < 0.5) {
      //high LOD, low D'
      boxColor = new Color(255, 224, 224);
    } else {
      //high LOD, high D' shades of red
      double blgr = (255-32)*2*(1-d);
      boxColor = new Color(255, (int) blgr, (int) blgr);
    }
  } else if (d > 0.99) {
    //high D', low LOD blueish color
    boxColor = new Color(192, 192, 240);
  } else {
    //no LD
    boxColor = Color.white;
  }
*/

int get_color_std(geno_cptr res)
{
  if ((res->dprime < -0.01) || (res->rsq2 < -0.01))
    {
      return 500;
    }
  else if (res->lod > 2) {
    if (res->dprime < 0.5) {
      /* high LOD, low D' */
      return 224;
    } else {
      /* high LOD, high D' shades of red */
      double blgr = (255-32)*2*(1 - res->dprime);
      return (int) blgr;
    }
  } else if (res->dprime > 0.99) {
    /* high D', low LOD blueish color */
    return 300;
  } else {
    /* no LD */
    return 400;
  }
}

/* Haploview's ESQ_SCHEME */
/* 
   double rsq = thisPair.getRSquared();
   Color boxColor = null;
   
   int r, g, b;
   
   r = g = b = (int)(255.0 * (1.0 - rsq));
   
   boxColor = new Color(r, g, b);
*/


int get_color_rsq(geno_cptr res)
{
  if ((res->dprime < -0.01) || (res->rsq2 < -0.01))
    {
      return 500;
    }
  double blgr =  255 *(1 - res->rsq2);
  return (1000 + (int) blgr);
}
