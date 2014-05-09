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


  /* careful with this macro: this goes to just after the 
     next null or space */
#define  goto_next_token(x) \
  while ( *(x) && (*(x) != ' ') && (*(x) != '\t')) { x++; }; x++

#define skip_3_tokens(x) \
  goto_next_token(x);\
  goto_next_token(x);\
  goto_next_token(x)

#define skip_5_tokens(x) \
  goto_next_token(x);\
  goto_next_token(x);\
  goto_next_token(x);\
  goto_next_token(x);\
  goto_next_token(x)
