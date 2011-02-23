/*
 * Copyright (C) 2007  Hin-Tak Leung
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA.
 */

/*
  Here is how the logic goes:
  We look at every sample of a particular snp.
  (1) if any sample contain an invalid input, the whole snp is invalid
  (2) if the whole snp, after considering all the sample, is multi-allelic,
      the whole snp is also invalid
  (3) if any allele is missing, the whole of that pair is missing/no call.
  (4) other than the cases above, we should have a definite call between
  AA, AB, and BB.
*/

/* we'll use the higher bit > 4 for INVALID bit mask */
#define INVALID_MASK      0x40
#define MULTIALLELIC_MASK 0x20
#define MISSING_MASK      0x10

#define BUILDIN_TABLE_LENGTH 11

static char buildin[BUILDIN_TABLE_LENGTH][2] = {
  {'A', 0x01},
  {'C', 0x02},
  {'G', 0x04},
  {'T', 0x08},

  {'N', MISSING_MASK},
  {'-', MISSING_MASK},

  {'1', 0x01},
  {'2', 0x02},
  {'3', 0x04},
  {'4', 0x08},

  {'0', MISSING_MASK},
};

static struct {
  char allele1 ;
  char allele2 ;
  char *message;
} resultmap[] = {
  {0x00, 0x00, "-/-"}, /* 0000 */ /* this is never used */
  {0x01, 0x00, "A/-"}, /* 0001 */
  {0x02, 0x00, "C/-"}, /* 0010 */
  {0x01, 0x02, "A/C"}, /* 0011 */
  {0x00, 0x04, "-/G"}, /* 0100 */ /* tested for the same as A/G */
  {0x01, 0x04, "A/G"}, /* 0101 */
  {0x02, 0x04, "C/G"}, /* 0110 */
  {MULTIALLELIC_MASK, MULTIALLELIC_MASK, ""   }, /* 0111 */
  {0x00, 0x08, "-/T"}, /* 1000 */ /* tested for the same as A/T */
  {0x01, 0x08, "A/T"}, /* 1001 */
  {0x02, 0x08, "C/T"}, /* 1010 */
  {MULTIALLELIC_MASK, MULTIALLELIC_MASK, ""   }, /* 1011 */
  {0x04, 0x08, "G/T"}, /* 1100 */
  {MULTIALLELIC_MASK, MULTIALLELIC_MASK, ""   }, /* 1101 */
  {MULTIALLELIC_MASK, MULTIALLELIC_MASK, ""   }, /* 1110 */
  {MULTIALLELIC_MASK, MULTIALLELIC_MASK, ""   }, /* 1111 */ 
};
