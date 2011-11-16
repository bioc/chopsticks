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

/* raw postscript code from fig2dev 
   and hand edited
 */

#define EPS_HEADER \
"%!PS-Adobe-2.0 EPSF-2.0\n\
%%Title: graphics.eps\n\
%%Creator: CIMR snp.matrix \n\
%%CreationDate: Sat Apr  1 00:00:00 2006\n\
%%For: nobody@example.com ()\n\
%%BoundingBox: "

#define EPS_BOUNDBOX(x, y, scale, eps_headroom) "0 0 %d %d\n", (int) ((x) * 12 * scale + 0.5), (int) (((y) * 6 + eps_headroom) * scale + 0.5)

/* 
perl -e 'for ($i = 0; $i < 256; $i++) {printf "%i %3.3f%s", $i, $i/255, "\n";}'
perl -e 'for ($i = 0; $i < 256; $i++) {printf "/p1%03i \[%3.3f %3.3f %3.3f\] def\\n\\%s", $i, $i/255,  $i/255,  $i/255,"\n";
perl -e 'for ($i = 0; $i < 256; $i++) {printf "/c1%03i \{%3.3f setgray\} bind def\\n\\%s", $i, $i/255, "\n";}'
*/

#define EPS_PROLOG \
"%%Magnification: 1.0000\n\
%%EndComments\n\
%%BeginProlog\n\
/pdfmark where\n\
 { pop globaldict /?pdfmark /exec load put } % pdfmark is built-in\n\
 {\n\
 globaldict\n\
  begin\n\
  /?pdfmark /pop load def % pdfmark is absent\n\
  /pdfmark /cleartomark load def\n\
  end\n\
 }\n\
ifelse\n\
%%EndProlog\n\
/$F2psDict 1000 dict def\n\
$F2psDict begin\n\
$F2psDict /mtrx matrix put\n\
/c500 {0.60 setgray} bind def\n\
/c400 {1 setgray} bind def\n\
/c300 {0.753 0.753 0.941 srgb} bind def\n\
/cb {0.0 0.0 1.0 srgb} bind def\n\
/c-1 {1.0 0.000 0.000 srgb} bind def\n\
/c0 {1.0 0.000 0.000 srgb} bind def\n\
/c1 {1.0 0.004 0.004 srgb} bind def\n\
/c2 {1.0 0.008 0.008 srgb} bind def\n\
/c3 {1.0 0.012 0.012 srgb} bind def\n\
/c4 {1.0 0.016 0.016 srgb} bind def\n\
/c5 {1.0 0.020 0.020 srgb} bind def\n\
/c6 {1.0 0.024 0.024 srgb} bind def\n\
/c7 {1.0 0.027 0.027 srgb} bind def\n\
/c8 {1.0 0.031 0.031 srgb} bind def\n\
/c9 {1.0 0.035 0.035 srgb} bind def\n\
/c10 {1.0 0.039 0.039 srgb} bind def\n\
/c11 {1.0 0.043 0.043 srgb} bind def\n\
/c12 {1.0 0.047 0.047 srgb} bind def\n\
/c13 {1.0 0.051 0.051 srgb} bind def\n\
/c14 {1.0 0.055 0.055 srgb} bind def\n\
/c15 {1.0 0.059 0.059 srgb} bind def\n\
/c16 {1.0 0.063 0.063 srgb} bind def\n\
/c17 {1.0 0.067 0.067 srgb} bind def\n\
/c18 {1.0 0.071 0.071 srgb} bind def\n\
/c19 {1.0 0.075 0.075 srgb} bind def\n\
/c20 {1.0 0.078 0.078 srgb} bind def\n\
/c21 {1.0 0.082 0.082 srgb} bind def\n\
/c22 {1.0 0.086 0.086 srgb} bind def\n\
/c23 {1.0 0.090 0.090 srgb} bind def\n\
/c24 {1.0 0.094 0.094 srgb} bind def\n\
/c25 {1.0 0.098 0.098 srgb} bind def\n\
/c26 {1.0 0.102 0.102 srgb} bind def\n\
/c27 {1.0 0.106 0.106 srgb} bind def\n\
/c28 {1.0 0.110 0.110 srgb} bind def\n\
/c29 {1.0 0.114 0.114 srgb} bind def\n\
/c30 {1.0 0.118 0.118 srgb} bind def\n\
/c31 {1.0 0.122 0.122 srgb} bind def\n\
/c32 {1.0 0.125 0.125 srgb} bind def\n\
/c33 {1.0 0.129 0.129 srgb} bind def\n\
/c34 {1.0 0.133 0.133 srgb} bind def\n\
/c35 {1.0 0.137 0.137 srgb} bind def\n\
/c36 {1.0 0.141 0.141 srgb} bind def\n\
/c37 {1.0 0.145 0.145 srgb} bind def\n\
/c38 {1.0 0.149 0.149 srgb} bind def\n\
/c39 {1.0 0.153 0.153 srgb} bind def\n\
/c40 {1.0 0.157 0.157 srgb} bind def\n\
/c41 {1.0 0.161 0.161 srgb} bind def\n\
/c42 {1.0 0.165 0.165 srgb} bind def\n\
/c43 {1.0 0.169 0.169 srgb} bind def\n\
/c44 {1.0 0.173 0.173 srgb} bind def\n\
/c45 {1.0 0.176 0.176 srgb} bind def\n\
/c46 {1.0 0.180 0.180 srgb} bind def\n\
/c47 {1.0 0.184 0.184 srgb} bind def\n\
/c48 {1.0 0.188 0.188 srgb} bind def\n\
/c49 {1.0 0.192 0.192 srgb} bind def\n\
/c50 {1.0 0.196 0.196 srgb} bind def\n\
/c51 {1.0 0.200 0.200 srgb} bind def\n\
/c52 {1.0 0.204 0.204 srgb} bind def\n\
/c53 {1.0 0.208 0.208 srgb} bind def\n\
/c54 {1.0 0.212 0.212 srgb} bind def\n\
/c55 {1.0 0.216 0.216 srgb} bind def\n\
/c56 {1.0 0.220 0.220 srgb} bind def\n\
/c57 {1.0 0.224 0.224 srgb} bind def\n\
/c58 {1.0 0.227 0.227 srgb} bind def\n\
/c59 {1.0 0.231 0.231 srgb} bind def\n\
/c60 {1.0 0.235 0.235 srgb} bind def\n\
/c61 {1.0 0.239 0.239 srgb} bind def\n\
/c62 {1.0 0.243 0.243 srgb} bind def\n\
/c63 {1.0 0.247 0.247 srgb} bind def\n\
/c64 {1.0 0.251 0.251 srgb} bind def\n\
/c65 {1.0 0.255 0.255 srgb} bind def\n\
/c66 {1.0 0.259 0.259 srgb} bind def\n\
/c67 {1.0 0.263 0.263 srgb} bind def\n\
/c68 {1.0 0.267 0.267 srgb} bind def\n\
/c69 {1.0 0.271 0.271 srgb} bind def\n\
/c70 {1.0 0.275 0.275 srgb} bind def\n\
/c71 {1.0 0.278 0.278 srgb} bind def\n\
/c72 {1.0 0.282 0.282 srgb} bind def\n\
/c73 {1.0 0.286 0.286 srgb} bind def\n\
/c74 {1.0 0.290 0.290 srgb} bind def\n\
/c75 {1.0 0.294 0.294 srgb} bind def\n\
/c76 {1.0 0.298 0.298 srgb} bind def\n\
/c77 {1.0 0.302 0.302 srgb} bind def\n\
/c78 {1.0 0.306 0.306 srgb} bind def\n\
/c79 {1.0 0.310 0.310 srgb} bind def\n\
/c80 {1.0 0.314 0.314 srgb} bind def\n\
/c81 {1.0 0.318 0.318 srgb} bind def\n\
/c82 {1.0 0.322 0.322 srgb} bind def\n\
/c83 {1.0 0.325 0.325 srgb} bind def\n\
/c84 {1.0 0.329 0.329 srgb} bind def\n\
/c85 {1.0 0.333 0.333 srgb} bind def\n\
/c86 {1.0 0.337 0.337 srgb} bind def\n\
/c87 {1.0 0.341 0.341 srgb} bind def\n\
/c88 {1.0 0.345 0.345 srgb} bind def\n\
/c89 {1.0 0.349 0.349 srgb} bind def\n\
/c90 {1.0 0.353 0.353 srgb} bind def\n\
/c91 {1.0 0.357 0.357 srgb} bind def\n\
/c92 {1.0 0.361 0.361 srgb} bind def\n\
/c93 {1.0 0.365 0.365 srgb} bind def\n\
/c94 {1.0 0.369 0.369 srgb} bind def\n\
/c95 {1.0 0.373 0.373 srgb} bind def\n\
/c96 {1.0 0.376 0.376 srgb} bind def\n\
/c97 {1.0 0.380 0.380 srgb} bind def\n\
/c98 {1.0 0.384 0.384 srgb} bind def\n\
/c99 {1.0 0.388 0.388 srgb} bind def\n\
/c100 {1.0 0.392 0.392 srgb} bind def\n\
/c101 {1.0 0.396 0.396 srgb} bind def\n\
/c102 {1.0 0.400 0.400 srgb} bind def\n\
/c103 {1.0 0.404 0.404 srgb} bind def\n\
/c104 {1.0 0.408 0.408 srgb} bind def\n\
/c105 {1.0 0.412 0.412 srgb} bind def\n\
/c106 {1.0 0.416 0.416 srgb} bind def\n\
/c107 {1.0 0.420 0.420 srgb} bind def\n\
/c108 {1.0 0.424 0.424 srgb} bind def\n\
/c109 {1.0 0.427 0.427 srgb} bind def\n\
/c110 {1.0 0.431 0.431 srgb} bind def\n\
/c111 {1.0 0.435 0.435 srgb} bind def\n\
/c112 {1.0 0.439 0.439 srgb} bind def\n\
/c113 {1.0 0.443 0.443 srgb} bind def\n\
/c114 {1.0 0.447 0.447 srgb} bind def\n\
/c115 {1.0 0.451 0.451 srgb} bind def\n\
/c116 {1.0 0.455 0.455 srgb} bind def\n\
/c117 {1.0 0.459 0.459 srgb} bind def\n\
/c118 {1.0 0.463 0.463 srgb} bind def\n\
/c119 {1.0 0.467 0.467 srgb} bind def\n\
/c120 {1.0 0.471 0.471 srgb} bind def\n\
/c121 {1.0 0.475 0.475 srgb} bind def\n\
/c122 {1.0 0.478 0.478 srgb} bind def\n\
/c123 {1.0 0.482 0.482 srgb} bind def\n\
/c124 {1.0 0.486 0.486 srgb} bind def\n\
/c125 {1.0 0.490 0.490 srgb} bind def\n\
/c126 {1.0 0.494 0.494 srgb} bind def\n\
/c127 {1.0 0.498 0.498 srgb} bind def\n\
/c128 {1.0 0.502 0.502 srgb} bind def\n\
/c129 {1.0 0.506 0.506 srgb} bind def\n\
/c130 {1.0 0.510 0.510 srgb} bind def\n\
/c131 {1.0 0.514 0.514 srgb} bind def\n\
/c132 {1.0 0.518 0.518 srgb} bind def\n\
/c133 {1.0 0.522 0.522 srgb} bind def\n\
/c134 {1.0 0.525 0.525 srgb} bind def\n\
/c135 {1.0 0.529 0.529 srgb} bind def\n\
/c136 {1.0 0.533 0.533 srgb} bind def\n\
/c137 {1.0 0.537 0.537 srgb} bind def\n\
/c138 {1.0 0.541 0.541 srgb} bind def\n\
/c139 {1.0 0.545 0.545 srgb} bind def\n\
/c140 {1.0 0.549 0.549 srgb} bind def\n\
/c141 {1.0 0.553 0.553 srgb} bind def\n\
/c142 {1.0 0.557 0.557 srgb} bind def\n\
/c143 {1.0 0.561 0.561 srgb} bind def\n\
/c144 {1.0 0.565 0.565 srgb} bind def\n\
/c145 {1.0 0.569 0.569 srgb} bind def\n\
/c146 {1.0 0.573 0.573 srgb} bind def\n\
/c147 {1.0 0.576 0.576 srgb} bind def\n\
/c148 {1.0 0.580 0.580 srgb} bind def\n\
/c149 {1.0 0.584 0.584 srgb} bind def\n\
/c150 {1.0 0.588 0.588 srgb} bind def\n\
/c151 {1.0 0.592 0.592 srgb} bind def\n\
/c152 {1.0 0.596 0.596 srgb} bind def\n\
/c153 {1.0 0.600 0.600 srgb} bind def\n\
/c154 {1.0 0.604 0.604 srgb} bind def\n\
/c155 {1.0 0.608 0.608 srgb} bind def\n\
/c156 {1.0 0.612 0.612 srgb} bind def\n\
/c157 {1.0 0.616 0.616 srgb} bind def\n\
/c158 {1.0 0.620 0.620 srgb} bind def\n\
/c159 {1.0 0.624 0.624 srgb} bind def\n\
/c160 {1.0 0.627 0.627 srgb} bind def\n\
/c161 {1.0 0.631 0.631 srgb} bind def\n\
/c162 {1.0 0.635 0.635 srgb} bind def\n\
/c163 {1.0 0.639 0.639 srgb} bind def\n\
/c164 {1.0 0.643 0.643 srgb} bind def\n\
/c165 {1.0 0.647 0.647 srgb} bind def\n\
/c166 {1.0 0.651 0.651 srgb} bind def\n\
/c167 {1.0 0.655 0.655 srgb} bind def\n\
/c168 {1.0 0.659 0.659 srgb} bind def\n\
/c169 {1.0 0.663 0.663 srgb} bind def\n\
/c170 {1.0 0.667 0.667 srgb} bind def\n\
/c171 {1.0 0.671 0.671 srgb} bind def\n\
/c172 {1.0 0.675 0.675 srgb} bind def\n\
/c173 {1.0 0.678 0.678 srgb} bind def\n\
/c174 {1.0 0.682 0.682 srgb} bind def\n\
/c175 {1.0 0.686 0.686 srgb} bind def\n\
/c176 {1.0 0.690 0.690 srgb} bind def\n\
/c177 {1.0 0.694 0.694 srgb} bind def\n\
/c178 {1.0 0.698 0.698 srgb} bind def\n\
/c179 {1.0 0.702 0.702 srgb} bind def\n\
/c180 {1.0 0.706 0.706 srgb} bind def\n\
/c181 {1.0 0.710 0.710 srgb} bind def\n\
/c182 {1.0 0.714 0.714 srgb} bind def\n\
/c183 {1.0 0.718 0.718 srgb} bind def\n\
/c184 {1.0 0.722 0.722 srgb} bind def\n\
/c185 {1.0 0.725 0.725 srgb} bind def\n\
/c186 {1.0 0.729 0.729 srgb} bind def\n\
/c187 {1.0 0.733 0.733 srgb} bind def\n\
/c188 {1.0 0.737 0.737 srgb} bind def\n\
/c189 {1.0 0.741 0.741 srgb} bind def\n\
/c190 {1.0 0.745 0.745 srgb} bind def\n\
/c191 {1.0 0.749 0.749 srgb} bind def\n\
/c192 {1.0 0.753 0.753 srgb} bind def\n\
/c193 {1.0 0.757 0.757 srgb} bind def\n\
/c194 {1.0 0.761 0.761 srgb} bind def\n\
/c195 {1.0 0.765 0.765 srgb} bind def\n\
/c196 {1.0 0.769 0.769 srgb} bind def\n\
/c197 {1.0 0.773 0.773 srgb} bind def\n\
/c198 {1.0 0.776 0.776 srgb} bind def\n\
/c199 {1.0 0.780 0.780 srgb} bind def\n\
/c200 {1.0 0.784 0.784 srgb} bind def\n\
/c201 {1.0 0.788 0.788 srgb} bind def\n\
/c202 {1.0 0.792 0.792 srgb} bind def\n\
/c203 {1.0 0.796 0.796 srgb} bind def\n\
/c204 {1.0 0.800 0.800 srgb} bind def\n\
/c205 {1.0 0.804 0.804 srgb} bind def\n\
/c206 {1.0 0.808 0.808 srgb} bind def\n\
/c207 {1.0 0.812 0.812 srgb} bind def\n\
/c208 {1.0 0.816 0.816 srgb} bind def\n\
/c209 {1.0 0.820 0.820 srgb} bind def\n\
/c210 {1.0 0.824 0.824 srgb} bind def\n\
/c211 {1.0 0.827 0.827 srgb} bind def\n\
/c212 {1.0 0.831 0.831 srgb} bind def\n\
/c213 {1.0 0.835 0.835 srgb} bind def\n\
/c214 {1.0 0.839 0.839 srgb} bind def\n\
/c215 {1.0 0.843 0.843 srgb} bind def\n\
/c216 {1.0 0.847 0.847 srgb} bind def\n\
/c217 {1.0 0.851 0.851 srgb} bind def\n\
/c218 {1.0 0.855 0.855 srgb} bind def\n\
/c219 {1.0 0.859 0.859 srgb} bind def\n\
/c220 {1.0 0.863 0.863 srgb} bind def\n\
/c221 {1.0 0.867 0.867 srgb} bind def\n\
/c222 {1.0 0.871 0.871 srgb} bind def\n\
/c223 {1.0 0.875 0.875 srgb} bind def\n\
/c224 {1.0 0.878 0.878 srgb} bind def\n\
%pdfmark definitions\n\
/p500 [0.6 0.6 0.6] def\n\
/p400 [1 1 1] def\n\
/p300 [0.753 0.753 0.941] def\n\
/pb {0.0 0.0 1.0 srgb} bind def\n\
/p-1 [1.0 0.000 0.000] def\n\
/p0 [1.0 0.000 0.000] def\n\
/p1 [1.0 0.004 0.004] def\n\
/p2 [1.0 0.008 0.008] def\n\
/p3 [1.0 0.012 0.012] def\n\
/p4 [1.0 0.016 0.016] def\n\
/p5 [1.0 0.020 0.020] def\n\
/p6 [1.0 0.024 0.024] def\n\
/p7 [1.0 0.027 0.027] def\n\
/p8 [1.0 0.031 0.031] def\n\
/p9 [1.0 0.035 0.035] def\n\
/p10 [1.0 0.039 0.039] def\n\
/p11 [1.0 0.043 0.043] def\n\
/p12 [1.0 0.047 0.047] def\n\
/p13 [1.0 0.051 0.051] def\n\
/p14 [1.0 0.055 0.055] def\n\
/p15 [1.0 0.059 0.059] def\n\
/p16 [1.0 0.063 0.063] def\n\
/p17 [1.0 0.067 0.067] def\n\
/p18 [1.0 0.071 0.071] def\n\
/p19 [1.0 0.075 0.075] def\n\
/p20 [1.0 0.078 0.078] def\n\
/p21 [1.0 0.082 0.082] def\n\
/p22 [1.0 0.086 0.086] def\n\
/p23 [1.0 0.090 0.090] def\n\
/p24 [1.0 0.094 0.094] def\n\
/p25 [1.0 0.098 0.098] def\n\
/p26 [1.0 0.102 0.102] def\n\
/p27 [1.0 0.106 0.106] def\n\
/p28 [1.0 0.110 0.110] def\n\
/p29 [1.0 0.114 0.114] def\n\
/p30 [1.0 0.118 0.118] def\n\
/p31 [1.0 0.122 0.122] def\n\
/p32 [1.0 0.125 0.125] def\n\
/p33 [1.0 0.129 0.129] def\n\
/p34 [1.0 0.133 0.133] def\n\
/p35 [1.0 0.137 0.137] def\n\
/p36 [1.0 0.141 0.141] def\n\
/p37 [1.0 0.145 0.145] def\n\
/p38 [1.0 0.149 0.149] def\n\
/p39 [1.0 0.153 0.153] def\n\
/p40 [1.0 0.157 0.157] def\n\
/p41 [1.0 0.161 0.161] def\n\
/p42 [1.0 0.165 0.165] def\n\
/p43 [1.0 0.169 0.169] def\n\
/p44 [1.0 0.173 0.173] def\n\
/p45 [1.0 0.176 0.176] def\n\
/p46 [1.0 0.180 0.180] def\n\
/p47 [1.0 0.184 0.184] def\n\
/p48 [1.0 0.188 0.188] def\n\
/p49 [1.0 0.192 0.192] def\n\
/p50 [1.0 0.196 0.196] def\n\
/p51 [1.0 0.200 0.200] def\n\
/p52 [1.0 0.204 0.204] def\n\
/p53 [1.0 0.208 0.208] def\n\
/p54 [1.0 0.212 0.212] def\n\
/p55 [1.0 0.216 0.216] def\n\
/p56 [1.0 0.220 0.220] def\n\
/p57 [1.0 0.224 0.224] def\n\
/p58 [1.0 0.227 0.227] def\n\
/p59 [1.0 0.231 0.231] def\n\
/p60 [1.0 0.235 0.235] def\n\
/p61 [1.0 0.239 0.239] def\n\
/p62 [1.0 0.243 0.243] def\n\
/p63 [1.0 0.247 0.247] def\n\
/p64 [1.0 0.251 0.251] def\n\
/p65 [1.0 0.255 0.255] def\n\
/p66 [1.0 0.259 0.259] def\n\
/p67 [1.0 0.263 0.263] def\n\
/p68 [1.0 0.267 0.267] def\n\
/p69 [1.0 0.271 0.271] def\n\
/p70 [1.0 0.275 0.275] def\n\
/p71 [1.0 0.278 0.278] def\n\
/p72 [1.0 0.282 0.282] def\n\
/p73 [1.0 0.286 0.286] def\n\
/p74 [1.0 0.290 0.290] def\n\
/p75 [1.0 0.294 0.294] def\n\
/p76 [1.0 0.298 0.298] def\n\
/p77 [1.0 0.302 0.302] def\n\
/p78 [1.0 0.306 0.306] def\n\
/p79 [1.0 0.310 0.310] def\n\
/p80 [1.0 0.314 0.314] def\n\
/p81 [1.0 0.318 0.318] def\n\
/p82 [1.0 0.322 0.322] def\n\
/p83 [1.0 0.325 0.325] def\n\
/p84 [1.0 0.329 0.329] def\n\
/p85 [1.0 0.333 0.333] def\n\
/p86 [1.0 0.337 0.337] def\n\
/p87 [1.0 0.341 0.341] def\n\
/p88 [1.0 0.345 0.345] def\n\
/p89 [1.0 0.349 0.349] def\n\
/p90 [1.0 0.353 0.353] def\n\
/p91 [1.0 0.357 0.357] def\n\
/p92 [1.0 0.361 0.361] def\n\
/p93 [1.0 0.365 0.365] def\n\
/p94 [1.0 0.369 0.369] def\n\
/p95 [1.0 0.373 0.373] def\n\
/p96 [1.0 0.376 0.376] def\n\
/p97 [1.0 0.380 0.380] def\n\
/p98 [1.0 0.384 0.384] def\n\
/p99 [1.0 0.388 0.388] def\n\
/p100 [1.0 0.392 0.392] def\n\
/p101 [1.0 0.396 0.396] def\n\
/p102 [1.0 0.400 0.400] def\n\
/p103 [1.0 0.404 0.404] def\n\
/p104 [1.0 0.408 0.408] def\n\
/p105 [1.0 0.412 0.412] def\n\
/p106 [1.0 0.416 0.416] def\n\
/p107 [1.0 0.420 0.420] def\n\
/p108 [1.0 0.424 0.424] def\n\
/p109 [1.0 0.427 0.427] def\n\
/p110 [1.0 0.431 0.431] def\n\
/p111 [1.0 0.435 0.435] def\n\
/p112 [1.0 0.439 0.439] def\n\
/p113 [1.0 0.443 0.443] def\n\
/p114 [1.0 0.447 0.447] def\n\
/p115 [1.0 0.451 0.451] def\n\
/p116 [1.0 0.455 0.455] def\n\
/p117 [1.0 0.459 0.459] def\n\
/p118 [1.0 0.463 0.463] def\n\
/p119 [1.0 0.467 0.467] def\n\
/p120 [1.0 0.471 0.471] def\n\
/p121 [1.0 0.475 0.475] def\n\
/p122 [1.0 0.478 0.478] def\n\
/p123 [1.0 0.482 0.482] def\n\
/p124 [1.0 0.486 0.486] def\n\
/p125 [1.0 0.490 0.490] def\n\
/p126 [1.0 0.494 0.494] def\n\
/p127 [1.0 0.498 0.498] def\n\
/p128 [1.0 0.502 0.502] def\n\
/p129 [1.0 0.506 0.506] def\n\
/p130 [1.0 0.510 0.510] def\n\
/p131 [1.0 0.514 0.514] def\n\
/p132 [1.0 0.518 0.518] def\n\
/p133 [1.0 0.522 0.522] def\n\
/p134 [1.0 0.525 0.525] def\n\
/p135 [1.0 0.529 0.529] def\n\
/p136 [1.0 0.533 0.533] def\n\
/p137 [1.0 0.537 0.537] def\n\
/p138 [1.0 0.541 0.541] def\n\
/p139 [1.0 0.545 0.545] def\n\
/p140 [1.0 0.549 0.549] def\n\
/p141 [1.0 0.553 0.553] def\n\
/p142 [1.0 0.557 0.557] def\n\
/p143 [1.0 0.561 0.561] def\n\
/p144 [1.0 0.565 0.565] def\n\
/p145 [1.0 0.569 0.569] def\n\
/p146 [1.0 0.573 0.573] def\n\
/p147 [1.0 0.576 0.576] def\n\
/p148 [1.0 0.580 0.580] def\n\
/p149 [1.0 0.584 0.584] def\n\
/p150 [1.0 0.588 0.588] def\n\
/p151 [1.0 0.592 0.592] def\n\
/p152 [1.0 0.596 0.596] def\n\
/p153 [1.0 0.600 0.600] def\n\
/p154 [1.0 0.604 0.604] def\n\
/p155 [1.0 0.608 0.608] def\n\
/p156 [1.0 0.612 0.612] def\n\
/p157 [1.0 0.616 0.616] def\n\
/p158 [1.0 0.620 0.620] def\n\
/p159 [1.0 0.624 0.624] def\n\
/p160 [1.0 0.627 0.627] def\n\
/p161 [1.0 0.631 0.631] def\n\
/p162 [1.0 0.635 0.635] def\n\
/p163 [1.0 0.639 0.639] def\n\
/p164 [1.0 0.643 0.643] def\n\
/p165 [1.0 0.647 0.647] def\n\
/p166 [1.0 0.651 0.651] def\n\
/p167 [1.0 0.655 0.655] def\n\
/p168 [1.0 0.659 0.659] def\n\
/p169 [1.0 0.663 0.663] def\n\
/p170 [1.0 0.667 0.667] def\n\
/p171 [1.0 0.671 0.671] def\n\
/p172 [1.0 0.675 0.675] def\n\
/p173 [1.0 0.678 0.678] def\n\
/p174 [1.0 0.682 0.682] def\n\
/p175 [1.0 0.686 0.686] def\n\
/p176 [1.0 0.690 0.690] def\n\
/p177 [1.0 0.694 0.694] def\n\
/p178 [1.0 0.698 0.698] def\n\
/p179 [1.0 0.702 0.702] def\n\
/p180 [1.0 0.706 0.706] def\n\
/p181 [1.0 0.710 0.710] def\n\
/p182 [1.0 0.714 0.714] def\n\
/p183 [1.0 0.718 0.718] def\n\
/p184 [1.0 0.722 0.722] def\n\
/p185 [1.0 0.725 0.725] def\n\
/p186 [1.0 0.729 0.729] def\n\
/p187 [1.0 0.733 0.733] def\n\
/p188 [1.0 0.737 0.737] def\n\
/p189 [1.0 0.741 0.741] def\n\
/p190 [1.0 0.745 0.745] def\n\
/p191 [1.0 0.749 0.749] def\n\
/p192 [1.0 0.753 0.753] def\n\
/p193 [1.0 0.757 0.757] def\n\
/p194 [1.0 0.761 0.761] def\n\
/p195 [1.0 0.765 0.765] def\n\
/p196 [1.0 0.769 0.769] def\n\
/p197 [1.0 0.773 0.773] def\n\
/p198 [1.0 0.776 0.776] def\n\
/p199 [1.0 0.780 0.780] def\n\
/p200 [1.0 0.784 0.784] def\n\
/p201 [1.0 0.788 0.788] def\n\
/p202 [1.0 0.792 0.792] def\n\
/p203 [1.0 0.796 0.796] def\n\
/p204 [1.0 0.800 0.800] def\n\
/p205 [1.0 0.804 0.804] def\n\
/p206 [1.0 0.808 0.808] def\n\
/p207 [1.0 0.812 0.812] def\n\
/p208 [1.0 0.816 0.816] def\n\
/p209 [1.0 0.820 0.820] def\n\
/p210 [1.0 0.824 0.824] def\n\
/p211 [1.0 0.827 0.827] def\n\
/p212 [1.0 0.831 0.831] def\n\
/p213 [1.0 0.835 0.835] def\n\
/p214 [1.0 0.839 0.839] def\n\
/p215 [1.0 0.843 0.843] def\n\
/p216 [1.0 0.847 0.847] def\n\
/p217 [1.0 0.851 0.851] def\n\
/p218 [1.0 0.855 0.855] def\n\
/p219 [1.0 0.859 0.859] def\n\
/p220 [1.0 0.863 0.863] def\n\
/p221 [1.0 0.867 0.867] def\n\
/p222 [1.0 0.871 0.871] def\n\
/p223 [1.0 0.875 0.875] def\n\
/p224 [1.0 0.878 0.878] def\n\
/c1000 {0.000 setgray} bind def\n\
/c1001 {0.004 setgray} bind def\n\
/c1002 {0.008 setgray} bind def\n\
/c1003 {0.012 setgray} bind def\n\
/c1004 {0.016 setgray} bind def\n\
/c1005 {0.020 setgray} bind def\n\
/c1006 {0.024 setgray} bind def\n\
/c1007 {0.027 setgray} bind def\n\
/c1008 {0.031 setgray} bind def\n\
/c1009 {0.035 setgray} bind def\n\
/c1010 {0.039 setgray} bind def\n\
/c1011 {0.043 setgray} bind def\n\
/c1012 {0.047 setgray} bind def\n\
/c1013 {0.051 setgray} bind def\n\
/c1014 {0.055 setgray} bind def\n\
/c1015 {0.059 setgray} bind def\n\
/c1016 {0.063 setgray} bind def\n\
/c1017 {0.067 setgray} bind def\n\
/c1018 {0.071 setgray} bind def\n\
/c1019 {0.075 setgray} bind def\n\
/c1020 {0.078 setgray} bind def\n\
/c1021 {0.082 setgray} bind def\n\
/c1022 {0.086 setgray} bind def\n\
/c1023 {0.090 setgray} bind def\n\
/c1024 {0.094 setgray} bind def\n\
/c1025 {0.098 setgray} bind def\n\
/c1026 {0.102 setgray} bind def\n\
/c1027 {0.106 setgray} bind def\n\
/c1028 {0.110 setgray} bind def\n\
/c1029 {0.114 setgray} bind def\n\
/c1030 {0.118 setgray} bind def\n\
/c1031 {0.122 setgray} bind def\n\
/c1032 {0.125 setgray} bind def\n\
/c1033 {0.129 setgray} bind def\n\
/c1034 {0.133 setgray} bind def\n\
/c1035 {0.137 setgray} bind def\n\
/c1036 {0.141 setgray} bind def\n\
/c1037 {0.145 setgray} bind def\n\
/c1038 {0.149 setgray} bind def\n\
/c1039 {0.153 setgray} bind def\n\
/c1040 {0.157 setgray} bind def\n\
/c1041 {0.161 setgray} bind def\n\
/c1042 {0.165 setgray} bind def\n\
/c1043 {0.169 setgray} bind def\n\
/c1044 {0.173 setgray} bind def\n\
/c1045 {0.176 setgray} bind def\n\
/c1046 {0.180 setgray} bind def\n\
/c1047 {0.184 setgray} bind def\n\
/c1048 {0.188 setgray} bind def\n\
/c1049 {0.192 setgray} bind def\n\
/c1050 {0.196 setgray} bind def\n\
/c1051 {0.200 setgray} bind def\n\
/c1052 {0.204 setgray} bind def\n\
/c1053 {0.208 setgray} bind def\n\
/c1054 {0.212 setgray} bind def\n\
/c1055 {0.216 setgray} bind def\n\
/c1056 {0.220 setgray} bind def\n\
/c1057 {0.224 setgray} bind def\n\
/c1058 {0.227 setgray} bind def\n\
/c1059 {0.231 setgray} bind def\n\
/c1060 {0.235 setgray} bind def\n\
/c1061 {0.239 setgray} bind def\n\
/c1062 {0.243 setgray} bind def\n\
/c1063 {0.247 setgray} bind def\n\
/c1064 {0.251 setgray} bind def\n\
/c1065 {0.255 setgray} bind def\n\
/c1066 {0.259 setgray} bind def\n\
/c1067 {0.263 setgray} bind def\n\
/c1068 {0.267 setgray} bind def\n\
/c1069 {0.271 setgray} bind def\n\
/c1070 {0.275 setgray} bind def\n\
/c1071 {0.278 setgray} bind def\n\
/c1072 {0.282 setgray} bind def\n\
/c1073 {0.286 setgray} bind def\n\
/c1074 {0.290 setgray} bind def\n\
/c1075 {0.294 setgray} bind def\n\
/c1076 {0.298 setgray} bind def\n\
/c1077 {0.302 setgray} bind def\n\
/c1078 {0.306 setgray} bind def\n\
/c1079 {0.310 setgray} bind def\n\
/c1080 {0.314 setgray} bind def\n\
/c1081 {0.318 setgray} bind def\n\
/c1082 {0.322 setgray} bind def\n\
/c1083 {0.325 setgray} bind def\n\
/c1084 {0.329 setgray} bind def\n\
/c1085 {0.333 setgray} bind def\n\
/c1086 {0.337 setgray} bind def\n\
/c1087 {0.341 setgray} bind def\n\
/c1088 {0.345 setgray} bind def\n\
/c1089 {0.349 setgray} bind def\n\
/c1090 {0.353 setgray} bind def\n\
/c1091 {0.357 setgray} bind def\n\
/c1092 {0.361 setgray} bind def\n\
/c1093 {0.365 setgray} bind def\n\
/c1094 {0.369 setgray} bind def\n\
/c1095 {0.373 setgray} bind def\n\
/c1096 {0.376 setgray} bind def\n\
/c1097 {0.380 setgray} bind def\n\
/c1098 {0.384 setgray} bind def\n\
/c1099 {0.388 setgray} bind def\n\
/c1100 {0.392 setgray} bind def\n\
/c1101 {0.396 setgray} bind def\n\
/c1102 {0.400 setgray} bind def\n\
/c1103 {0.404 setgray} bind def\n\
/c1104 {0.408 setgray} bind def\n\
/c1105 {0.412 setgray} bind def\n\
/c1106 {0.416 setgray} bind def\n\
/c1107 {0.420 setgray} bind def\n\
/c1108 {0.424 setgray} bind def\n\
/c1109 {0.427 setgray} bind def\n\
/c1110 {0.431 setgray} bind def\n\
/c1111 {0.435 setgray} bind def\n\
/c1112 {0.439 setgray} bind def\n\
/c1113 {0.443 setgray} bind def\n\
/c1114 {0.447 setgray} bind def\n\
/c1115 {0.451 setgray} bind def\n\
/c1116 {0.455 setgray} bind def\n\
/c1117 {0.459 setgray} bind def\n\
/c1118 {0.463 setgray} bind def\n\
/c1119 {0.467 setgray} bind def\n\
/c1120 {0.471 setgray} bind def\n\
/c1121 {0.475 setgray} bind def\n\
/c1122 {0.478 setgray} bind def\n\
/c1123 {0.482 setgray} bind def\n\
/c1124 {0.486 setgray} bind def\n\
/c1125 {0.490 setgray} bind def\n\
/c1126 {0.494 setgray} bind def\n\
/c1127 {0.498 setgray} bind def\n\
/c1128 {0.502 setgray} bind def\n\
/c1129 {0.506 setgray} bind def\n\
/c1130 {0.510 setgray} bind def\n\
/c1131 {0.514 setgray} bind def\n\
/c1132 {0.518 setgray} bind def\n\
/c1133 {0.522 setgray} bind def\n\
/c1134 {0.525 setgray} bind def\n\
/c1135 {0.529 setgray} bind def\n\
/c1136 {0.533 setgray} bind def\n\
/c1137 {0.537 setgray} bind def\n\
/c1138 {0.541 setgray} bind def\n\
/c1139 {0.545 setgray} bind def\n\
/c1140 {0.549 setgray} bind def\n\
/c1141 {0.553 setgray} bind def\n\
/c1142 {0.557 setgray} bind def\n\
/c1143 {0.561 setgray} bind def\n\
/c1144 {0.565 setgray} bind def\n\
/c1145 {0.569 setgray} bind def\n\
/c1146 {0.573 setgray} bind def\n\
/c1147 {0.576 setgray} bind def\n\
/c1148 {0.580 setgray} bind def\n\
/c1149 {0.584 setgray} bind def\n\
/c1150 {0.588 setgray} bind def\n\
/c1151 {0.592 setgray} bind def\n\
/c1152 {0.596 setgray} bind def\n\
/c1153 {0.600 setgray} bind def\n\
/c1154 {0.604 setgray} bind def\n\
/c1155 {0.608 setgray} bind def\n\
/c1156 {0.612 setgray} bind def\n\
/c1157 {0.616 setgray} bind def\n\
/c1158 {0.620 setgray} bind def\n\
/c1159 {0.624 setgray} bind def\n\
/c1160 {0.627 setgray} bind def\n\
/c1161 {0.631 setgray} bind def\n\
/c1162 {0.635 setgray} bind def\n\
/c1163 {0.639 setgray} bind def\n\
/c1164 {0.643 setgray} bind def\n\
/c1165 {0.647 setgray} bind def\n\
/c1166 {0.651 setgray} bind def\n\
/c1167 {0.655 setgray} bind def\n\
/c1168 {0.659 setgray} bind def\n\
/c1169 {0.663 setgray} bind def\n\
/c1170 {0.667 setgray} bind def\n\
/c1171 {0.671 setgray} bind def\n\
/c1172 {0.675 setgray} bind def\n\
/c1173 {0.678 setgray} bind def\n\
/c1174 {0.682 setgray} bind def\n\
/c1175 {0.686 setgray} bind def\n\
/c1176 {0.690 setgray} bind def\n\
/c1177 {0.694 setgray} bind def\n\
/c1178 {0.698 setgray} bind def\n\
/c1179 {0.702 setgray} bind def\n\
/c1180 {0.706 setgray} bind def\n\
/c1181 {0.710 setgray} bind def\n\
/c1182 {0.714 setgray} bind def\n\
/c1183 {0.718 setgray} bind def\n\
/c1184 {0.722 setgray} bind def\n\
/c1185 {0.725 setgray} bind def\n\
/c1186 {0.729 setgray} bind def\n\
/c1187 {0.733 setgray} bind def\n\
/c1188 {0.737 setgray} bind def\n\
/c1189 {0.741 setgray} bind def\n\
/c1190 {0.745 setgray} bind def\n\
/c1191 {0.749 setgray} bind def\n\
/c1192 {0.753 setgray} bind def\n\
/c1193 {0.757 setgray} bind def\n\
/c1194 {0.761 setgray} bind def\n\
/c1195 {0.765 setgray} bind def\n\
/c1196 {0.769 setgray} bind def\n\
/c1197 {0.773 setgray} bind def\n\
/c1198 {0.776 setgray} bind def\n\
/c1199 {0.780 setgray} bind def\n\
/c1200 {0.784 setgray} bind def\n\
/c1201 {0.788 setgray} bind def\n\
/c1202 {0.792 setgray} bind def\n\
/c1203 {0.796 setgray} bind def\n\
/c1204 {0.800 setgray} bind def\n\
/c1205 {0.804 setgray} bind def\n\
/c1206 {0.808 setgray} bind def\n\
/c1207 {0.812 setgray} bind def\n\
/c1208 {0.816 setgray} bind def\n\
/c1209 {0.820 setgray} bind def\n\
/c1210 {0.824 setgray} bind def\n\
/c1211 {0.827 setgray} bind def\n\
/c1212 {0.831 setgray} bind def\n\
/c1213 {0.835 setgray} bind def\n\
/c1214 {0.839 setgray} bind def\n\
/c1215 {0.843 setgray} bind def\n\
/c1216 {0.847 setgray} bind def\n\
/c1217 {0.851 setgray} bind def\n\
/c1218 {0.855 setgray} bind def\n\
/c1219 {0.859 setgray} bind def\n\
/c1220 {0.863 setgray} bind def\n\
/c1221 {0.867 setgray} bind def\n\
/c1222 {0.871 setgray} bind def\n\
/c1223 {0.875 setgray} bind def\n\
/c1224 {0.878 setgray} bind def\n\
/c1225 {0.882 setgray} bind def\n\
/c1226 {0.886 setgray} bind def\n\
/c1227 {0.890 setgray} bind def\n\
/c1228 {0.894 setgray} bind def\n\
/c1229 {0.898 setgray} bind def\n\
/c1230 {0.902 setgray} bind def\n\
/c1231 {0.906 setgray} bind def\n\
/c1232 {0.910 setgray} bind def\n\
/c1233 {0.914 setgray} bind def\n\
/c1234 {0.918 setgray} bind def\n\
/c1235 {0.922 setgray} bind def\n\
/c1236 {0.925 setgray} bind def\n\
/c1237 {0.929 setgray} bind def\n\
/c1238 {0.933 setgray} bind def\n\
/c1239 {0.937 setgray} bind def\n\
/c1240 {0.941 setgray} bind def\n\
/c1241 {0.945 setgray} bind def\n\
/c1242 {0.949 setgray} bind def\n\
/c1243 {0.953 setgray} bind def\n\
/c1244 {0.957 setgray} bind def\n\
/c1245 {0.961 setgray} bind def\n\
/c1246 {0.965 setgray} bind def\n\
/c1247 {0.969 setgray} bind def\n\
/c1248 {0.973 setgray} bind def\n\
/c1249 {0.976 setgray} bind def\n\
/c1250 {0.980 setgray} bind def\n\
/c1251 {0.984 setgray} bind def\n\
/c1252 {0.988 setgray} bind def\n\
/c1253 {0.992 setgray} bind def\n\
/c1254 {0.996 setgray} bind def\n\
/c1255 {1.000 setgray} bind def\n\
/p1000 [0.000 0.000 0.000] def\n\
/p1001 [0.004 0.004 0.004] def\n\
/p1002 [0.008 0.008 0.008] def\n\
/p1003 [0.012 0.012 0.012] def\n\
/p1004 [0.016 0.016 0.016] def\n\
/p1005 [0.020 0.020 0.020] def\n\
/p1006 [0.024 0.024 0.024] def\n\
/p1007 [0.027 0.027 0.027] def\n\
/p1008 [0.031 0.031 0.031] def\n\
/p1009 [0.035 0.035 0.035] def\n\
/p1010 [0.039 0.039 0.039] def\n\
/p1011 [0.043 0.043 0.043] def\n\
/p1012 [0.047 0.047 0.047] def\n\
/p1013 [0.051 0.051 0.051] def\n\
/p1014 [0.055 0.055 0.055] def\n\
/p1015 [0.059 0.059 0.059] def\n\
/p1016 [0.063 0.063 0.063] def\n\
/p1017 [0.067 0.067 0.067] def\n\
/p1018 [0.071 0.071 0.071] def\n\
/p1019 [0.075 0.075 0.075] def\n\
/p1020 [0.078 0.078 0.078] def\n\
/p1021 [0.082 0.082 0.082] def\n\
/p1022 [0.086 0.086 0.086] def\n\
/p1023 [0.090 0.090 0.090] def\n\
/p1024 [0.094 0.094 0.094] def\n\
/p1025 [0.098 0.098 0.098] def\n\
/p1026 [0.102 0.102 0.102] def\n\
/p1027 [0.106 0.106 0.106] def\n\
/p1028 [0.110 0.110 0.110] def\n\
/p1029 [0.114 0.114 0.114] def\n\
/p1030 [0.118 0.118 0.118] def\n\
/p1031 [0.122 0.122 0.122] def\n\
/p1032 [0.125 0.125 0.125] def\n\
/p1033 [0.129 0.129 0.129] def\n\
/p1034 [0.133 0.133 0.133] def\n\
/p1035 [0.137 0.137 0.137] def\n\
/p1036 [0.141 0.141 0.141] def\n\
/p1037 [0.145 0.145 0.145] def\n\
/p1038 [0.149 0.149 0.149] def\n\
/p1039 [0.153 0.153 0.153] def\n\
/p1040 [0.157 0.157 0.157] def\n\
/p1041 [0.161 0.161 0.161] def\n\
/p1042 [0.165 0.165 0.165] def\n\
/p1043 [0.169 0.169 0.169] def\n\
/p1044 [0.173 0.173 0.173] def\n\
/p1045 [0.176 0.176 0.176] def\n\
/p1046 [0.180 0.180 0.180] def\n\
/p1047 [0.184 0.184 0.184] def\n\
/p1048 [0.188 0.188 0.188] def\n\
/p1049 [0.192 0.192 0.192] def\n\
/p1050 [0.196 0.196 0.196] def\n\
/p1051 [0.200 0.200 0.200] def\n\
/p1052 [0.204 0.204 0.204] def\n\
/p1053 [0.208 0.208 0.208] def\n\
/p1054 [0.212 0.212 0.212] def\n\
/p1055 [0.216 0.216 0.216] def\n\
/p1056 [0.220 0.220 0.220] def\n\
/p1057 [0.224 0.224 0.224] def\n\
/p1058 [0.227 0.227 0.227] def\n\
/p1059 [0.231 0.231 0.231] def\n\
/p1060 [0.235 0.235 0.235] def\n\
/p1061 [0.239 0.239 0.239] def\n\
/p1062 [0.243 0.243 0.243] def\n\
/p1063 [0.247 0.247 0.247] def\n\
/p1064 [0.251 0.251 0.251] def\n\
/p1065 [0.255 0.255 0.255] def\n\
/p1066 [0.259 0.259 0.259] def\n\
/p1067 [0.263 0.263 0.263] def\n\
/p1068 [0.267 0.267 0.267] def\n\
/p1069 [0.271 0.271 0.271] def\n\
/p1070 [0.275 0.275 0.275] def\n\
/p1071 [0.278 0.278 0.278] def\n\
/p1072 [0.282 0.282 0.282] def\n\
/p1073 [0.286 0.286 0.286] def\n\
/p1074 [0.290 0.290 0.290] def\n\
/p1075 [0.294 0.294 0.294] def\n\
/p1076 [0.298 0.298 0.298] def\n\
/p1077 [0.302 0.302 0.302] def\n\
/p1078 [0.306 0.306 0.306] def\n\
/p1079 [0.310 0.310 0.310] def\n\
/p1080 [0.314 0.314 0.314] def\n\
/p1081 [0.318 0.318 0.318] def\n\
/p1082 [0.322 0.322 0.322] def\n\
/p1083 [0.325 0.325 0.325] def\n\
/p1084 [0.329 0.329 0.329] def\n\
/p1085 [0.333 0.333 0.333] def\n\
/p1086 [0.337 0.337 0.337] def\n\
/p1087 [0.341 0.341 0.341] def\n\
/p1088 [0.345 0.345 0.345] def\n\
/p1089 [0.349 0.349 0.349] def\n\
/p1090 [0.353 0.353 0.353] def\n\
/p1091 [0.357 0.357 0.357] def\n\
/p1092 [0.361 0.361 0.361] def\n\
/p1093 [0.365 0.365 0.365] def\n\
/p1094 [0.369 0.369 0.369] def\n\
/p1095 [0.373 0.373 0.373] def\n\
/p1096 [0.376 0.376 0.376] def\n\
/p1097 [0.380 0.380 0.380] def\n\
/p1098 [0.384 0.384 0.384] def\n\
/p1099 [0.388 0.388 0.388] def\n\
/p1100 [0.392 0.392 0.392] def\n\
/p1101 [0.396 0.396 0.396] def\n\
/p1102 [0.400 0.400 0.400] def\n\
/p1103 [0.404 0.404 0.404] def\n\
/p1104 [0.408 0.408 0.408] def\n\
/p1105 [0.412 0.412 0.412] def\n\
/p1106 [0.416 0.416 0.416] def\n\
/p1107 [0.420 0.420 0.420] def\n\
/p1108 [0.424 0.424 0.424] def\n\
/p1109 [0.427 0.427 0.427] def\n\
/p1110 [0.431 0.431 0.431] def\n\
/p1111 [0.435 0.435 0.435] def\n\
/p1112 [0.439 0.439 0.439] def\n\
/p1113 [0.443 0.443 0.443] def\n\
/p1114 [0.447 0.447 0.447] def\n\
/p1115 [0.451 0.451 0.451] def\n\
/p1116 [0.455 0.455 0.455] def\n\
/p1117 [0.459 0.459 0.459] def\n\
/p1118 [0.463 0.463 0.463] def\n\
/p1119 [0.467 0.467 0.467] def\n\
/p1120 [0.471 0.471 0.471] def\n\
/p1121 [0.475 0.475 0.475] def\n\
/p1122 [0.478 0.478 0.478] def\n\
/p1123 [0.482 0.482 0.482] def\n\
/p1124 [0.486 0.486 0.486] def\n\
/p1125 [0.490 0.490 0.490] def\n\
/p1126 [0.494 0.494 0.494] def\n\
/p1127 [0.498 0.498 0.498] def\n\
/p1128 [0.502 0.502 0.502] def\n\
/p1129 [0.506 0.506 0.506] def\n\
/p1130 [0.510 0.510 0.510] def\n\
/p1131 [0.514 0.514 0.514] def\n\
/p1132 [0.518 0.518 0.518] def\n\
/p1133 [0.522 0.522 0.522] def\n\
/p1134 [0.525 0.525 0.525] def\n\
/p1135 [0.529 0.529 0.529] def\n\
/p1136 [0.533 0.533 0.533] def\n\
/p1137 [0.537 0.537 0.537] def\n\
/p1138 [0.541 0.541 0.541] def\n\
/p1139 [0.545 0.545 0.545] def\n\
/p1140 [0.549 0.549 0.549] def\n\
/p1141 [0.553 0.553 0.553] def\n\
/p1142 [0.557 0.557 0.557] def\n\
/p1143 [0.561 0.561 0.561] def\n\
/p1144 [0.565 0.565 0.565] def\n\
/p1145 [0.569 0.569 0.569] def\n\
/p1146 [0.573 0.573 0.573] def\n\
/p1147 [0.576 0.576 0.576] def\n\
/p1148 [0.580 0.580 0.580] def\n\
/p1149 [0.584 0.584 0.584] def\n\
/p1150 [0.588 0.588 0.588] def\n\
/p1151 [0.592 0.592 0.592] def\n\
/p1152 [0.596 0.596 0.596] def\n\
/p1153 [0.600 0.600 0.600] def\n\
/p1154 [0.604 0.604 0.604] def\n\
/p1155 [0.608 0.608 0.608] def\n\
/p1156 [0.612 0.612 0.612] def\n\
/p1157 [0.616 0.616 0.616] def\n\
/p1158 [0.620 0.620 0.620] def\n\
/p1159 [0.624 0.624 0.624] def\n\
/p1160 [0.627 0.627 0.627] def\n\
/p1161 [0.631 0.631 0.631] def\n\
/p1162 [0.635 0.635 0.635] def\n\
/p1163 [0.639 0.639 0.639] def\n\
/p1164 [0.643 0.643 0.643] def\n\
/p1165 [0.647 0.647 0.647] def\n\
/p1166 [0.651 0.651 0.651] def\n\
/p1167 [0.655 0.655 0.655] def\n\
/p1168 [0.659 0.659 0.659] def\n\
/p1169 [0.663 0.663 0.663] def\n\
/p1170 [0.667 0.667 0.667] def\n\
/p1171 [0.671 0.671 0.671] def\n\
/p1172 [0.675 0.675 0.675] def\n\
/p1173 [0.678 0.678 0.678] def\n\
/p1174 [0.682 0.682 0.682] def\n\
/p1175 [0.686 0.686 0.686] def\n\
/p1176 [0.690 0.690 0.690] def\n\
/p1177 [0.694 0.694 0.694] def\n\
/p1178 [0.698 0.698 0.698] def\n\
/p1179 [0.702 0.702 0.702] def\n\
/p1180 [0.706 0.706 0.706] def\n\
/p1181 [0.710 0.710 0.710] def\n\
/p1182 [0.714 0.714 0.714] def\n\
/p1183 [0.718 0.718 0.718] def\n\
/p1184 [0.722 0.722 0.722] def\n\
/p1185 [0.725 0.725 0.725] def\n\
/p1186 [0.729 0.729 0.729] def\n\
/p1187 [0.733 0.733 0.733] def\n\
/p1188 [0.737 0.737 0.737] def\n\
/p1189 [0.741 0.741 0.741] def\n\
/p1190 [0.745 0.745 0.745] def\n\
/p1191 [0.749 0.749 0.749] def\n\
/p1192 [0.753 0.753 0.753] def\n\
/p1193 [0.757 0.757 0.757] def\n\
/p1194 [0.761 0.761 0.761] def\n\
/p1195 [0.765 0.765 0.765] def\n\
/p1196 [0.769 0.769 0.769] def\n\
/p1197 [0.773 0.773 0.773] def\n\
/p1198 [0.776 0.776 0.776] def\n\
/p1199 [0.780 0.780 0.780] def\n\
/p1200 [0.784 0.784 0.784] def\n\
/p1201 [0.788 0.788 0.788] def\n\
/p1202 [0.792 0.792 0.792] def\n\
/p1203 [0.796 0.796 0.796] def\n\
/p1204 [0.800 0.800 0.800] def\n\
/p1205 [0.804 0.804 0.804] def\n\
/p1206 [0.808 0.808 0.808] def\n\
/p1207 [0.812 0.812 0.812] def\n\
/p1208 [0.816 0.816 0.816] def\n\
/p1209 [0.820 0.820 0.820] def\n\
/p1210 [0.824 0.824 0.824] def\n\
/p1211 [0.827 0.827 0.827] def\n\
/p1212 [0.831 0.831 0.831] def\n\
/p1213 [0.835 0.835 0.835] def\n\
/p1214 [0.839 0.839 0.839] def\n\
/p1215 [0.843 0.843 0.843] def\n\
/p1216 [0.847 0.847 0.847] def\n\
/p1217 [0.851 0.851 0.851] def\n\
/p1218 [0.855 0.855 0.855] def\n\
/p1219 [0.859 0.859 0.859] def\n\
/p1220 [0.863 0.863 0.863] def\n\
/p1221 [0.867 0.867 0.867] def\n\
/p1222 [0.871 0.871 0.871] def\n\
/p1223 [0.875 0.875 0.875] def\n\
/p1224 [0.878 0.878 0.878] def\n\
/p1225 [0.882 0.882 0.882] def\n\
/p1226 [0.886 0.886 0.886] def\n\
/p1227 [0.890 0.890 0.890] def\n\
/p1228 [0.894 0.894 0.894] def\n\
/p1229 [0.898 0.898 0.898] def\n\
/p1230 [0.902 0.902 0.902] def\n\
/p1231 [0.906 0.906 0.906] def\n\
/p1232 [0.910 0.910 0.910] def\n\
/p1233 [0.914 0.914 0.914] def\n\
/p1234 [0.918 0.918 0.918] def\n\
/p1235 [0.922 0.922 0.922] def\n\
/p1236 [0.925 0.925 0.925] def\n\
/p1237 [0.929 0.929 0.929] def\n\
/p1238 [0.933 0.933 0.933] def\n\
/p1239 [0.937 0.937 0.937] def\n\
/p1240 [0.941 0.941 0.941] def\n\
/p1241 [0.945 0.945 0.945] def\n\
/p1242 [0.949 0.949 0.949] def\n\
/p1243 [0.953 0.953 0.953] def\n\
/p1244 [0.957 0.957 0.957] def\n\
/p1245 [0.961 0.961 0.961] def\n\
/p1246 [0.965 0.965 0.965] def\n\
/p1247 [0.969 0.969 0.969] def\n\
/p1248 [0.973 0.973 0.973] def\n\
/p1249 [0.976 0.976 0.976] def\n\
/p1250 [0.980 0.980 0.980] def\n\
/p1251 [0.984 0.984 0.984] def\n\
/p1252 [0.988 0.988 0.988] def\n\
/p1253 [0.992 0.992 0.992] def\n\
/p1254 [0.996 0.996 0.996] def\n\
/p1255 [1.000 1.000 1.000] def\n\
end\n\
/cp {closepath} bind def\n\
/ef {eofill} bind def\n\
/gr {grestore} bind def\n\
/gs {gsave} bind def\n\
/sa {save} bind def\n\
/rs {restore} bind def\n\
/l {lineto} bind def\n\
/rl {rlineto} bind def\n\
/m {moveto} bind def\n\
/rm {rmoveto} bind def\n\
/n {newpath} bind def\n\
/s {stroke} bind def\n\
/sh {show} bind def\n\
/slc {setlinecap} bind def\n\
/slj {setlinejoin} bind def\n\
/slw {setlinewidth} bind def\n\
/srgb {setrgbcolor} bind def\n\
/rot {rotate} bind def\n\
/sc {scale} bind def\n\
/ff {findfont} bind def\n\
/sf {setfont} bind def\n\
/scf {scalefont} bind def\n\
/sw {stringwidth} bind def\n\
/tr {translate} bind def\n\
/tnt {dup dup currentrgbcolor\n\
  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n\
  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n\
  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}\n\
  bind def\n\
/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul\n\
  4 -2 roll mul srgb} bind def\n\
/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def\n\
/$F2psEnd {$F2psEnteredState restore end} def\n\
sa\n"

#define EPS_SCALE_TO_FIT(x) "%.6f %.6f sc\n", (x), (x)

#define EPS_CLIP_PATH(x,y, eps_headroom)  "n 0 %d m 0 0 l %d 0 l %d %d l cp clip n\n", \
((y)*6 + eps_headroom), ((x)*12), ((x)*12), ((y)*6 + eps_headroom)

#define EPS_BODY \
"1 -1 sc\n\
$F2psBegin\n\
10 setmiterlimit\n\
0 slj 0 slc\n\
 2 2 sc\n\
/Times-Roman ff 3.6 scf sf\n\
0.075 slw\n"

#define EPS_TRANSLATE(x,y) "%d -%d tr\n",  ((x) *3 + 6), ((y) *3 + 3)
#define EPS_TEXT(x, text) "%d 0.75 m gs 1 -1 sc  90.0 rot (%s) cb sh gr\n", ((x) *6 - 2), (text)

#define EPS_METRIC(x, trans_x) "n %d -20 m %.1f -40 l gs c1128 s gr n  %.1f -40 m %.1f -50 l gs 0.15 slw c1000 s gr\n", \
((x) *6 - 2), ((trans_x) *6 - 2), ((trans_x) *6 - 2), ((trans_x) *6 - 2)


#define EPS_RHOMBUS(x, y, color) \
  "n %d %d m -3 3 rl 3 3 rl 3 -3 rl cp gs c%d 1.00 shd ef gr gs c500 s gr\n" , ((x) *6 + (y) * 3), ((y) *3), (color)

#define EPS_NOTE(x, y, dprime, rsq2, lod, color) \
"[ /Rect [%d %d %d %d]\n\
/Subtype /Circle\n\
/Title (pair %d, %d)\n\
/Contents (d'=%.3f\n\
r^2=%.3f\n\
lod=%.3f)\n\
/C p%d\n\
/ANN pdfmark\n" , \
((x) *6 + (y) * 3 -1), ((y) *3 + 4), ((x) *6 + (y) * 3 +1), ((y) *3 + 2), (x +1), (y + x +2), (dprime), \
    (rsq2), (lod), (color)

#define EPS_ENDING \
"% here ends figure;\n\
$F2psEnd\n\
rs\n\
showpage\n"
