// This file was generated by autolevclean.c and contains C code generated by Autolev

#include <math.h>
#include "gait2dc.h"
void gait2dc_FK_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],
   double fk[NFK], double dfk_dq[NFK][NDOF],
   double fkdot[NFK], double dfkdot_dq[NFK][NDOF]) {
	double q1 = q[0];
	double q1p = qd[0];
	double q1pp = qdd[0];
	double q2 = q[1];
	double q2p = qd[1];
	double q2pp = qdd[1];
	double q3 = q[2];
	double q3p = qd[2];
	double q3pp = qdd[2];
	double q4 = q[3];
	double q4p = qd[3];
	double q4pp = qdd[3];
	double q5 = q[4];
	double q5p = qd[4];
	double q5pp = qdd[4];
	double q6 = q[5];
	double q6p = qd[5];
	double q6pp = qdd[5];
	double q7 = q[6];
	double q7p = qd[6];
	double q7pp = qdd[6];
	double q8 = q[7];
	double q8p = qd[7];
	double q8pp = qdd[7];
	double q9 = q[8];
	double q9p = qd[8];
	double q9pp = qdd[8];
	static double z[199];
  z[1] = cos(q3);
  z[2] = sin(q3);
  z[3] = cos(q4);
  z[4] = sin(q4);
  z[5] = cos(q5);
  z[6] = sin(q5);
  z[7] = z[1]*z[3] - z[2]*z[4];
  z[8] = -z[1]*z[4] - z[2]*z[3];
  z[9] = z[1]*z[4] + z[2]*z[3];
  z[12] = cos(q6);
  z[13] = sin(q6);
  z[14] = z[5]*z[9] + z[6]*z[7];
  z[15] = z[5]*z[7] + z[6]*z[8];
  z[17] = z[5]*z[7] - z[6]*z[9];
  z[18] = z[5]*z[8] - z[6]*z[7];
  z[25] = cos(q7);
  z[26] = sin(q7);
  z[27] = cos(q8);
  z[28] = sin(q8);
  z[29] = z[1]*z[25] - z[2]*z[26];
  z[30] = -z[1]*z[26] - z[2]*z[25];
  z[31] = z[1]*z[26] + z[2]*z[25];
  z[34] = cos(q9);
  z[35] = sin(q9);
  z[36] = z[27]*z[31] + z[28]*z[29];
  z[37] = z[27]*z[29] + z[28]*z[30];
  z[39] = z[27]*z[29] - z[28]*z[31];
  z[40] = z[27]*z[30] - z[28]*z[29];
  z[47] = z[12]*z[15] + z[13]*z[18];
  z[48] = z[12]*z[18] - z[13]*z[15];
  z[49] = z[12]*z[14] + z[13]*z[17];
  z[50] = z[12]*z[17] - z[13]*z[14];
  z[51] = z[34]*z[37] + z[35]*z[40];
  z[52] = z[34]*z[40] - z[35]*z[37];
  z[53] = z[34]*z[36] + z[35]*z[39];
  z[54] = z[34]*z[39] - z[35]*z[36];
  z[55] = z[2]*z[4] - z[1]*z[3];
  z[56] = z[5]*z[8] + z[6]*z[55];
  z[57] = z[5]*z[55] - z[6]*z[8];
  z[58] = -z[5]*z[7] - z[6]*z[8];
  z[59] = -z[5]*z[9] - z[6]*z[7];
  z[60] = z[12]*z[56] + z[13]*z[57];
  z[61] = z[12]*z[18] + z[13]*z[58];
  z[62] = z[12]*z[57] - z[13]*z[56];
  z[63] = z[12]*z[58] - z[13]*z[18];
  z[64] = -z[12]*z[15] - z[13]*z[18];
  z[65] = z[12]*z[17] + z[13]*z[59];
  z[66] = z[12]*z[59] - z[13]*z[17];
  z[67] = -z[12]*z[14] - z[13]*z[17];
  z[68] = z[2]*z[26] - z[1]*z[25];
  z[69] = z[27]*z[30] + z[28]*z[68];
  z[70] = z[27]*z[68] - z[28]*z[30];
  z[71] = -z[27]*z[29] - z[28]*z[30];
  z[72] = -z[27]*z[31] - z[28]*z[29];
  z[73] = z[34]*z[69] + z[35]*z[70];
  z[74] = z[34]*z[40] + z[35]*z[71];
  z[75] = z[34]*z[70] - z[35]*z[69];
  z[76] = z[34]*z[71] - z[35]*z[40];
  z[77] = -z[34]*z[37] - z[35]*z[40];
  z[78] = z[34]*z[39] + z[35]*z[72];
  z[79] = z[34]*z[72] - z[35]*z[39];
  z[80] = -z[34]*z[36] - z[35]*z[39];
  z[81] = par->ThighLen*z[55];
  z[82] = par->ThighLen*z[8];
  z[83] = -par->ShankLen*z[57] - par->ThighLen*z[55];
  z[84] = par->ShankLen*z[58];
  z[85] = -par->ShankLen*z[18] - par->ThighLen*z[8];
  z[86] = par->ShankLen*z[59];
  z[87] = par->ThighLen*z[68];
  z[88] = par->ThighLen*z[30];
  z[89] = -par->ShankLen*z[70] - par->ThighLen*z[68];
  z[90] = par->ShankLen*z[71];
  z[91] = -par->ShankLen*z[40] - par->ThighLen*z[30];
  z[92] = par->ShankLen*z[72];
  z[93] = z[2]*q3p;
  z[94] = z[1]*q3p;
  z[95] = z[8]*(q3p+q4p);
  z[96] = z[55]*(q3p+q4p);
  z[97] = z[7]*(q3p+q4p);
  z[98] = q1p - par->ThighLen*z[55]*(q3p+q4p);
  z[99] = q2p - par->ThighLen*z[8]*(q3p+q4p);
  z[100] = z[18]*q5p + z[56]*q3p + z[56]*q4p;
  z[101] = z[57]*q3p + z[57]*q4p + z[58]*q5p;
  z[102] = z[15]*q3p + z[15]*q4p + z[17]*q5p;
  z[103] = z[18]*q3p + z[18]*q4p + z[59]*q5p;
  z[104] = q1p - par->ThighLen*z[55]*(q3p+q4p) - par->ShankLen*(z[57]*q3p+
  z[57]*q4p+z[58]*q5p);
  z[105] = q2p - par->ThighLen*z[8]*(q3p+q4p) - par->ShankLen*(z[18]*q3p+
  z[18]*q4p+z[59]*q5p);
  z[106] = z[48]*q6p + z[60]*q3p + z[60]*q4p + z[61]*q5p;
  z[107] = z[62]*q3p + z[62]*q4p + z[63]*q5p + z[64]*q6p;
  z[108] = z[47]*q3p + z[47]*q4p + z[50]*q6p + z[65]*q5p;
  z[109] = z[48]*q3p + z[48]*q4p + z[66]*q5p + z[67]*q6p;
  z[110] = z[30]*(q3p+q7p);
  z[111] = z[68]*(q3p+q7p);
  z[112] = z[29]*(q3p+q7p);
  z[113] = q1p - par->ThighLen*z[68]*(q3p+q7p);
  z[114] = q2p - par->ThighLen*z[30]*(q3p+q7p);
  z[115] = z[40]*q8p + z[69]*q3p + z[69]*q7p;
  z[116] = z[70]*q3p + z[70]*q7p + z[71]*q8p;
  z[117] = z[37]*q3p + z[37]*q7p + z[39]*q8p;
  z[118] = z[40]*q3p + z[40]*q7p + z[72]*q8p;
  z[119] = q1p - par->ThighLen*z[68]*(q3p+q7p) - par->ShankLen*(z[70]*q3p+
  z[70]*q7p+z[71]*q8p);
  z[120] = q2p - par->ThighLen*z[30]*(q3p+q7p) - par->ShankLen*(z[40]*q3p+
  z[40]*q7p+z[72]*q8p);
  z[121] = z[52]*q9p + z[73]*q3p + z[73]*q7p + z[74]*q8p;
  z[122] = z[75]*q3p + z[75]*q7p + z[76]*q8p + z[77]*q9p;
  z[123] = z[51]*q3p + z[51]*q7p + z[54]*q9p + z[78]*q8p;
  z[124] = z[52]*q3p + z[52]*q7p + z[79]*q8p + z[80]*q9p;
  z[125] = par->ThighLen*z[9]*(q3p+q4p);
  z[126] = par->ThighLen*z[55]*(q3p+q4p);
  z[127] = z[5]*z[55] + z[6]*z[9];
  z[128] = z[57]*q5p + z[127]*q3p + z[127]*q4p;
  z[129] = z[5]*z[9] - z[6]*z[55];
  z[130] = -z[5]*z[8] - z[6]*z[55];
  z[131] = z[129]*q3p + z[129]*q4p + z[130]*q5p;
  z[132] = z[6]*z[7] - z[5]*z[8];
  z[133] = z[130]*q3p + z[130]*q4p + z[132]*q5p;
  z[134] = z[6]*z[9] - z[5]*z[7];
  z[135] = z[58]*q3p + z[58]*q4p + z[134]*q5p;
  z[136] = -par->ThighLen*z[9]*(q3p+q4p) - par->ShankLen*(z[129]*q3p+z[129]*
  q4p+z[130]*q5p);
  z[137] = par->ShankLen*(z[130]*q3p+z[130]*q4p+z[132]*q5p);
  z[138] = -par->ThighLen*z[55]*(q3p+q4p) - par->ShankLen*(z[57]*q3p+z[57]*
  q4p+z[58]*q5p);
  z[139] = par->ShankLen*(z[58]*q3p+z[58]*q4p+z[134]*q5p);
  z[140] = z[12]*z[127] + z[13]*z[129];
  z[141] = z[12]*z[57] + z[13]*z[130];
  z[142] = z[62]*q6p + z[140]*q3p + z[140]*q4p + z[141]*q5p;
  z[143] = z[12]*z[58] + z[13]*z[132];
  z[144] = z[63]*q6p + z[141]*q3p + z[141]*q4p + z[143]*q5p;
  z[145] = z[12]*z[129] - z[13]*z[127];
  z[146] = z[12]*z[130] - z[13]*z[57];
  z[147] = -z[12]*z[56] - z[13]*z[57];
  z[148] = z[145]*q3p + z[145]*q4p + z[146]*q5p + z[147]*q6p;
  z[149] = z[12]*z[132] - z[13]*z[58];
  z[150] = -z[12]*z[18] - z[13]*z[58];
  z[151] = z[146]*q3p + z[146]*q4p + z[149]*q5p + z[150]*q6p;
  z[152] = z[13]*z[15] - z[12]*z[18];
  z[153] = z[147]*q3p + z[147]*q4p + z[150]*q5p + z[152]*q6p;
  z[154] = z[12]*z[59] + z[13]*z[134];
  z[155] = z[61]*q3p + z[61]*q4p + z[66]*q6p + z[154]*q5p;
  z[156] = z[12]*z[134] - z[13]*z[59];
  z[157] = -z[12]*z[17] - z[13]*z[59];
  z[158] = z[63]*q3p + z[63]*q4p + z[156]*q5p + z[157]*q6p;
  z[159] = z[13]*z[14] - z[12]*z[17];
  z[160] = z[64]*q3p + z[64]*q4p + z[157]*q5p + z[159]*q6p;
  z[161] = par->ThighLen*z[31]*(q3p+q7p);
  z[162] = par->ThighLen*z[68]*(q3p+q7p);
  z[163] = z[27]*z[68] + z[28]*z[31];
  z[164] = z[70]*q8p + z[163]*q3p + z[163]*q7p;
  z[165] = z[27]*z[31] - z[28]*z[68];
  z[166] = -z[27]*z[30] - z[28]*z[68];
  z[167] = z[165]*q3p + z[165]*q7p + z[166]*q8p;
  z[168] = z[28]*z[29] - z[27]*z[30];
  z[169] = z[166]*q3p + z[166]*q7p + z[168]*q8p;
  z[170] = z[28]*z[31] - z[27]*z[29];
  z[171] = z[71]*q3p + z[71]*q7p + z[170]*q8p;
  z[172] = -par->ThighLen*z[31]*(q3p+q7p) - par->ShankLen*(z[165]*q3p+z[165]*
  q7p+z[166]*q8p);
  z[173] = par->ShankLen*(z[166]*q3p+z[166]*q7p+z[168]*q8p);
  z[174] = -par->ThighLen*z[68]*(q3p+q7p) - par->ShankLen*(z[70]*q3p+z[70]*
  q7p+z[71]*q8p);
  z[175] = par->ShankLen*(z[71]*q3p+z[71]*q7p+z[170]*q8p);
  z[176] = z[34]*z[163] + z[35]*z[165];
  z[177] = z[34]*z[70] + z[35]*z[166];
  z[178] = z[75]*q9p + z[176]*q3p + z[176]*q7p + z[177]*q8p;
  z[179] = z[34]*z[71] + z[35]*z[168];
  z[180] = z[76]*q9p + z[177]*q3p + z[177]*q7p + z[179]*q8p;
  z[181] = z[34]*z[165] - z[35]*z[163];
  z[182] = z[34]*z[166] - z[35]*z[70];
  z[183] = -z[34]*z[69] - z[35]*z[70];
  z[184] = z[181]*q3p + z[181]*q7p + z[182]*q8p + z[183]*q9p;
  z[185] = z[34]*z[168] - z[35]*z[71];
  z[186] = -z[34]*z[40] - z[35]*z[71];
  z[187] = z[182]*q3p + z[182]*q7p + z[185]*q8p + z[186]*q9p;
  z[188] = z[35]*z[37] - z[34]*z[40];
  z[189] = z[183]*q3p + z[183]*q7p + z[186]*q8p + z[188]*q9p;
  z[190] = z[34]*z[72] + z[35]*z[170];
  z[191] = z[74]*q3p + z[74]*q7p + z[79]*q9p + z[190]*q8p;
  z[192] = z[34]*z[170] - z[35]*z[72];
  z[193] = -z[34]*z[39] - z[35]*z[72];
  z[194] = z[76]*q3p + z[76]*q7p + z[192]*q8p + z[193]*q9p;
  z[195] = z[35]*z[36] - z[34]*z[39];
  z[196] = z[77]*q3p + z[77]*q7p + z[193]*q8p + z[195]*q9p;
  z[197] = z[9]*(q3p+q4p);
  z[198] = z[31]*(q3p+q7p);


  fk[0] = q1;
  fk[1] = q2;
  fk[2] = z[1];
  fk[3] = -z[2];
  fk[4] = z[2];
  fk[5] = z[1];
  fk[6] = q1;
  fk[7] = q2;
  fk[8] = z[7];
  fk[9] = z[8];
  fk[10] = z[9];
  fk[11] = z[7];
  fk[12] = q1 - par->ThighLen*z[8];
  fk[13] = q2 - par->ThighLen*z[7];
  fk[14] = z[15];
  fk[15] = z[18];
  fk[16] = z[14];
  fk[17] = z[17];
  fk[18] = q1 - par->ShankLen*z[18] - par->ThighLen*z[8];
  fk[19] = q2 - par->ShankLen*z[17] - par->ThighLen*z[7];
  fk[20] = z[47];
  fk[21] = z[48];
  fk[22] = z[49];
  fk[23] = z[50];
  fk[24] = q1;
  fk[25] = q2;
  fk[26] = z[29];
  fk[27] = z[30];
  fk[28] = z[31];
  fk[29] = z[29];
  fk[30] = q1 - par->ThighLen*z[30];
  fk[31] = q2 - par->ThighLen*z[29];
  fk[32] = z[37];
  fk[33] = z[40];
  fk[34] = z[36];
  fk[35] = z[39];
  fk[36] = q1 - par->ShankLen*z[40] - par->ThighLen*z[30];
  fk[37] = q2 - par->ShankLen*z[39] - par->ThighLen*z[29];
  fk[38] = z[51];
  fk[39] = z[52];
  fk[40] = z[53];
  fk[41] = z[54];
  dfk_dq[0][0] = 1;
  dfk_dq[0][1] = 0;
  dfk_dq[0][2] = 0;
  dfk_dq[0][3] = 0;
  dfk_dq[0][4] = 0;
  dfk_dq[0][5] = 0;
  dfk_dq[0][6] = 0;
  dfk_dq[0][7] = 0;
  dfk_dq[0][8] = 0;
  dfk_dq[1][0] = 0;
  dfk_dq[1][1] = 1;
  dfk_dq[1][2] = 0;
  dfk_dq[1][3] = 0;
  dfk_dq[1][4] = 0;
  dfk_dq[1][5] = 0;
  dfk_dq[1][6] = 0;
  dfk_dq[1][7] = 0;
  dfk_dq[1][8] = 0;
  dfk_dq[2][0] = 0;
  dfk_dq[2][1] = 0;
  dfk_dq[2][2] = -z[2];
  dfk_dq[2][3] = 0;
  dfk_dq[2][4] = 0;
  dfk_dq[2][5] = 0;
  dfk_dq[2][6] = 0;
  dfk_dq[2][7] = 0;
  dfk_dq[2][8] = 0;
  dfk_dq[3][0] = 0;
  dfk_dq[3][1] = 0;
  dfk_dq[3][2] = -z[1];
  dfk_dq[3][3] = 0;
  dfk_dq[3][4] = 0;
  dfk_dq[3][5] = 0;
  dfk_dq[3][6] = 0;
  dfk_dq[3][7] = 0;
  dfk_dq[3][8] = 0;
  dfk_dq[4][0] = 0;
  dfk_dq[4][1] = 0;
  dfk_dq[4][2] = z[1];
  dfk_dq[4][3] = 0;
  dfk_dq[4][4] = 0;
  dfk_dq[4][5] = 0;
  dfk_dq[4][6] = 0;
  dfk_dq[4][7] = 0;
  dfk_dq[4][8] = 0;
  dfk_dq[5][0] = 0;
  dfk_dq[5][1] = 0;
  dfk_dq[5][2] = -z[2];
  dfk_dq[5][3] = 0;
  dfk_dq[5][4] = 0;
  dfk_dq[5][5] = 0;
  dfk_dq[5][6] = 0;
  dfk_dq[5][7] = 0;
  dfk_dq[5][8] = 0;
  dfk_dq[6][0] = 1;
  dfk_dq[6][1] = 0;
  dfk_dq[6][2] = 0;
  dfk_dq[6][3] = 0;
  dfk_dq[6][4] = 0;
  dfk_dq[6][5] = 0;
  dfk_dq[6][6] = 0;
  dfk_dq[6][7] = 0;
  dfk_dq[6][8] = 0;
  dfk_dq[7][0] = 0;
  dfk_dq[7][1] = 1;
  dfk_dq[7][2] = 0;
  dfk_dq[7][3] = 0;
  dfk_dq[7][4] = 0;
  dfk_dq[7][5] = 0;
  dfk_dq[7][6] = 0;
  dfk_dq[7][7] = 0;
  dfk_dq[7][8] = 0;
  dfk_dq[8][0] = 0;
  dfk_dq[8][1] = 0;
  dfk_dq[8][2] = z[8];
  dfk_dq[8][3] = z[8];
  dfk_dq[8][4] = 0;
  dfk_dq[8][5] = 0;
  dfk_dq[8][6] = 0;
  dfk_dq[8][7] = 0;
  dfk_dq[8][8] = 0;
  dfk_dq[9][0] = 0;
  dfk_dq[9][1] = 0;
  dfk_dq[9][2] = z[55];
  dfk_dq[9][3] = z[55];
  dfk_dq[9][4] = 0;
  dfk_dq[9][5] = 0;
  dfk_dq[9][6] = 0;
  dfk_dq[9][7] = 0;
  dfk_dq[9][8] = 0;
  dfk_dq[10][0] = 0;
  dfk_dq[10][1] = 0;
  dfk_dq[10][2] = z[7];
  dfk_dq[10][3] = z[7];
  dfk_dq[10][4] = 0;
  dfk_dq[10][5] = 0;
  dfk_dq[10][6] = 0;
  dfk_dq[10][7] = 0;
  dfk_dq[10][8] = 0;
  dfk_dq[11][0] = 0;
  dfk_dq[11][1] = 0;
  dfk_dq[11][2] = z[8];
  dfk_dq[11][3] = z[8];
  dfk_dq[11][4] = 0;
  dfk_dq[11][5] = 0;
  dfk_dq[11][6] = 0;
  dfk_dq[11][7] = 0;
  dfk_dq[11][8] = 0;
  dfk_dq[12][0] = 1;
  dfk_dq[12][1] = 0;
  dfk_dq[12][2] = -z[81];
  dfk_dq[12][3] = -z[81];
  dfk_dq[12][4] = 0;
  dfk_dq[12][5] = 0;
  dfk_dq[12][6] = 0;
  dfk_dq[12][7] = 0;
  dfk_dq[12][8] = 0;
  dfk_dq[13][0] = 0;
  dfk_dq[13][1] = 1;
  dfk_dq[13][2] = -z[82];
  dfk_dq[13][3] = -z[82];
  dfk_dq[13][4] = 0;
  dfk_dq[13][5] = 0;
  dfk_dq[13][6] = 0;
  dfk_dq[13][7] = 0;
  dfk_dq[13][8] = 0;
  dfk_dq[14][0] = 0;
  dfk_dq[14][1] = 0;
  dfk_dq[14][2] = z[56];
  dfk_dq[14][3] = z[56];
  dfk_dq[14][4] = z[18];
  dfk_dq[14][5] = 0;
  dfk_dq[14][6] = 0;
  dfk_dq[14][7] = 0;
  dfk_dq[14][8] = 0;
  dfk_dq[15][0] = 0;
  dfk_dq[15][1] = 0;
  dfk_dq[15][2] = z[57];
  dfk_dq[15][3] = z[57];
  dfk_dq[15][4] = z[58];
  dfk_dq[15][5] = 0;
  dfk_dq[15][6] = 0;
  dfk_dq[15][7] = 0;
  dfk_dq[15][8] = 0;
  dfk_dq[16][0] = 0;
  dfk_dq[16][1] = 0;
  dfk_dq[16][2] = z[15];
  dfk_dq[16][3] = z[15];
  dfk_dq[16][4] = z[17];
  dfk_dq[16][5] = 0;
  dfk_dq[16][6] = 0;
  dfk_dq[16][7] = 0;
  dfk_dq[16][8] = 0;
  dfk_dq[17][0] = 0;
  dfk_dq[17][1] = 0;
  dfk_dq[17][2] = z[18];
  dfk_dq[17][3] = z[18];
  dfk_dq[17][4] = z[59];
  dfk_dq[17][5] = 0;
  dfk_dq[17][6] = 0;
  dfk_dq[17][7] = 0;
  dfk_dq[17][8] = 0;
  dfk_dq[18][0] = 1;
  dfk_dq[18][1] = 0;
  dfk_dq[18][2] = z[83];
  dfk_dq[18][3] = z[83];
  dfk_dq[18][4] = -z[84];
  dfk_dq[18][5] = 0;
  dfk_dq[18][6] = 0;
  dfk_dq[18][7] = 0;
  dfk_dq[18][8] = 0;
  dfk_dq[19][0] = 0;
  dfk_dq[19][1] = 1;
  dfk_dq[19][2] = z[85];
  dfk_dq[19][3] = z[85];
  dfk_dq[19][4] = -z[86];
  dfk_dq[19][5] = 0;
  dfk_dq[19][6] = 0;
  dfk_dq[19][7] = 0;
  dfk_dq[19][8] = 0;
  dfk_dq[20][0] = 0;
  dfk_dq[20][1] = 0;
  dfk_dq[20][2] = z[60];
  dfk_dq[20][3] = z[60];
  dfk_dq[20][4] = z[61];
  dfk_dq[20][5] = z[48];
  dfk_dq[20][6] = 0;
  dfk_dq[20][7] = 0;
  dfk_dq[20][8] = 0;
  dfk_dq[21][0] = 0;
  dfk_dq[21][1] = 0;
  dfk_dq[21][2] = z[62];
  dfk_dq[21][3] = z[62];
  dfk_dq[21][4] = z[63];
  dfk_dq[21][5] = z[64];
  dfk_dq[21][6] = 0;
  dfk_dq[21][7] = 0;
  dfk_dq[21][8] = 0;
  dfk_dq[22][0] = 0;
  dfk_dq[22][1] = 0;
  dfk_dq[22][2] = z[47];
  dfk_dq[22][3] = z[47];
  dfk_dq[22][4] = z[65];
  dfk_dq[22][5] = z[50];
  dfk_dq[22][6] = 0;
  dfk_dq[22][7] = 0;
  dfk_dq[22][8] = 0;
  dfk_dq[23][0] = 0;
  dfk_dq[23][1] = 0;
  dfk_dq[23][2] = z[48];
  dfk_dq[23][3] = z[48];
  dfk_dq[23][4] = z[66];
  dfk_dq[23][5] = z[67];
  dfk_dq[23][6] = 0;
  dfk_dq[23][7] = 0;
  dfk_dq[23][8] = 0;
  dfk_dq[24][0] = 1;
  dfk_dq[24][1] = 0;
  dfk_dq[24][2] = 0;
  dfk_dq[24][3] = 0;
  dfk_dq[24][4] = 0;
  dfk_dq[24][5] = 0;
  dfk_dq[24][6] = 0;
  dfk_dq[24][7] = 0;
  dfk_dq[24][8] = 0;
  dfk_dq[25][0] = 0;
  dfk_dq[25][1] = 1;
  dfk_dq[25][2] = 0;
  dfk_dq[25][3] = 0;
  dfk_dq[25][4] = 0;
  dfk_dq[25][5] = 0;
  dfk_dq[25][6] = 0;
  dfk_dq[25][7] = 0;
  dfk_dq[25][8] = 0;
  dfk_dq[26][0] = 0;
  dfk_dq[26][1] = 0;
  dfk_dq[26][2] = z[30];
  dfk_dq[26][3] = 0;
  dfk_dq[26][4] = 0;
  dfk_dq[26][5] = 0;
  dfk_dq[26][6] = z[30];
  dfk_dq[26][7] = 0;
  dfk_dq[26][8] = 0;
  dfk_dq[27][0] = 0;
  dfk_dq[27][1] = 0;
  dfk_dq[27][2] = z[68];
  dfk_dq[27][3] = 0;
  dfk_dq[27][4] = 0;
  dfk_dq[27][5] = 0;
  dfk_dq[27][6] = z[68];
  dfk_dq[27][7] = 0;
  dfk_dq[27][8] = 0;
  dfk_dq[28][0] = 0;
  dfk_dq[28][1] = 0;
  dfk_dq[28][2] = z[29];
  dfk_dq[28][3] = 0;
  dfk_dq[28][4] = 0;
  dfk_dq[28][5] = 0;
  dfk_dq[28][6] = z[29];
  dfk_dq[28][7] = 0;
  dfk_dq[28][8] = 0;
  dfk_dq[29][0] = 0;
  dfk_dq[29][1] = 0;
  dfk_dq[29][2] = z[30];
  dfk_dq[29][3] = 0;
  dfk_dq[29][4] = 0;
  dfk_dq[29][5] = 0;
  dfk_dq[29][6] = z[30];
  dfk_dq[29][7] = 0;
  dfk_dq[29][8] = 0;
  dfk_dq[30][0] = 1;
  dfk_dq[30][1] = 0;
  dfk_dq[30][2] = -z[87];
  dfk_dq[30][3] = 0;
  dfk_dq[30][4] = 0;
  dfk_dq[30][5] = 0;
  dfk_dq[30][6] = -z[87];
  dfk_dq[30][7] = 0;
  dfk_dq[30][8] = 0;
  dfk_dq[31][0] = 0;
  dfk_dq[31][1] = 1;
  dfk_dq[31][2] = -z[88];
  dfk_dq[31][3] = 0;
  dfk_dq[31][4] = 0;
  dfk_dq[31][5] = 0;
  dfk_dq[31][6] = -z[88];
  dfk_dq[31][7] = 0;
  dfk_dq[31][8] = 0;
  dfk_dq[32][0] = 0;
  dfk_dq[32][1] = 0;
  dfk_dq[32][2] = z[69];
  dfk_dq[32][3] = 0;
  dfk_dq[32][4] = 0;
  dfk_dq[32][5] = 0;
  dfk_dq[32][6] = z[69];
  dfk_dq[32][7] = z[40];
  dfk_dq[32][8] = 0;
  dfk_dq[33][0] = 0;
  dfk_dq[33][1] = 0;
  dfk_dq[33][2] = z[70];
  dfk_dq[33][3] = 0;
  dfk_dq[33][4] = 0;
  dfk_dq[33][5] = 0;
  dfk_dq[33][6] = z[70];
  dfk_dq[33][7] = z[71];
  dfk_dq[33][8] = 0;
  dfk_dq[34][0] = 0;
  dfk_dq[34][1] = 0;
  dfk_dq[34][2] = z[37];
  dfk_dq[34][3] = 0;
  dfk_dq[34][4] = 0;
  dfk_dq[34][5] = 0;
  dfk_dq[34][6] = z[37];
  dfk_dq[34][7] = z[39];
  dfk_dq[34][8] = 0;
  dfk_dq[35][0] = 0;
  dfk_dq[35][1] = 0;
  dfk_dq[35][2] = z[40];
  dfk_dq[35][3] = 0;
  dfk_dq[35][4] = 0;
  dfk_dq[35][5] = 0;
  dfk_dq[35][6] = z[40];
  dfk_dq[35][7] = z[72];
  dfk_dq[35][8] = 0;
  dfk_dq[36][0] = 1;
  dfk_dq[36][1] = 0;
  dfk_dq[36][2] = z[89];
  dfk_dq[36][3] = 0;
  dfk_dq[36][4] = 0;
  dfk_dq[36][5] = 0;
  dfk_dq[36][6] = z[89];
  dfk_dq[36][7] = -z[90];
  dfk_dq[36][8] = 0;
  dfk_dq[37][0] = 0;
  dfk_dq[37][1] = 1;
  dfk_dq[37][2] = z[91];
  dfk_dq[37][3] = 0;
  dfk_dq[37][4] = 0;
  dfk_dq[37][5] = 0;
  dfk_dq[37][6] = z[91];
  dfk_dq[37][7] = -z[92];
  dfk_dq[37][8] = 0;
  dfk_dq[38][0] = 0;
  dfk_dq[38][1] = 0;
  dfk_dq[38][2] = z[73];
  dfk_dq[38][3] = 0;
  dfk_dq[38][4] = 0;
  dfk_dq[38][5] = 0;
  dfk_dq[38][6] = z[73];
  dfk_dq[38][7] = z[74];
  dfk_dq[38][8] = z[52];
  dfk_dq[39][0] = 0;
  dfk_dq[39][1] = 0;
  dfk_dq[39][2] = z[75];
  dfk_dq[39][3] = 0;
  dfk_dq[39][4] = 0;
  dfk_dq[39][5] = 0;
  dfk_dq[39][6] = z[75];
  dfk_dq[39][7] = z[76];
  dfk_dq[39][8] = z[77];
  dfk_dq[40][0] = 0;
  dfk_dq[40][1] = 0;
  dfk_dq[40][2] = z[51];
  dfk_dq[40][3] = 0;
  dfk_dq[40][4] = 0;
  dfk_dq[40][5] = 0;
  dfk_dq[40][6] = z[51];
  dfk_dq[40][7] = z[78];
  dfk_dq[40][8] = z[54];
  dfk_dq[41][0] = 0;
  dfk_dq[41][1] = 0;
  dfk_dq[41][2] = z[52];
  dfk_dq[41][3] = 0;
  dfk_dq[41][4] = 0;
  dfk_dq[41][5] = 0;
  dfk_dq[41][6] = z[52];
  dfk_dq[41][7] = z[79];
  dfk_dq[41][8] = z[80];
  fkdot[0] = q1p;
  fkdot[1] = q2p;
  fkdot[2] = -z[93];
  fkdot[3] = -z[94];
  fkdot[4] = z[94];
  fkdot[5] = -z[93];
  fkdot[6] = q1p;
  fkdot[7] = q2p;
  fkdot[8] = z[95];
  fkdot[9] = z[96];
  fkdot[10] = z[97];
  fkdot[11] = z[95];
  fkdot[12] = z[98];
  fkdot[13] = z[99];
  fkdot[14] = z[100];
  fkdot[15] = z[101];
  fkdot[16] = z[102];
  fkdot[17] = z[103];
  fkdot[18] = z[104];
  fkdot[19] = z[105];
  fkdot[20] = z[106];
  fkdot[21] = z[107];
  fkdot[22] = z[108];
  fkdot[23] = z[109];
  fkdot[24] = q1p;
  fkdot[25] = q2p;
  fkdot[26] = z[110];
  fkdot[27] = z[111];
  fkdot[28] = z[112];
  fkdot[29] = z[110];
  fkdot[30] = z[113];
  fkdot[31] = z[114];
  fkdot[32] = z[115];
  fkdot[33] = z[116];
  fkdot[34] = z[117];
  fkdot[35] = z[118];
  fkdot[36] = z[119];
  fkdot[37] = z[120];
  fkdot[38] = z[121];
  fkdot[39] = z[122];
  fkdot[40] = z[123];
  fkdot[41] = z[124];
  dfkdot_dq[0][0] = 0;
  dfkdot_dq[0][1] = 0;
  dfkdot_dq[0][2] = 0;
  dfkdot_dq[0][3] = 0;
  dfkdot_dq[0][4] = 0;
  dfkdot_dq[0][5] = 0;
  dfkdot_dq[0][6] = 0;
  dfkdot_dq[0][7] = 0;
  dfkdot_dq[0][8] = 0;
  dfkdot_dq[1][0] = 0;
  dfkdot_dq[1][1] = 0;
  dfkdot_dq[1][2] = 0;
  dfkdot_dq[1][3] = 0;
  dfkdot_dq[1][4] = 0;
  dfkdot_dq[1][5] = 0;
  dfkdot_dq[1][6] = 0;
  dfkdot_dq[1][7] = 0;
  dfkdot_dq[1][8] = 0;
  dfkdot_dq[2][0] = 0;
  dfkdot_dq[2][1] = 0;
  dfkdot_dq[2][2] = -z[94];
  dfkdot_dq[2][3] = 0;
  dfkdot_dq[2][4] = 0;
  dfkdot_dq[2][5] = 0;
  dfkdot_dq[2][6] = 0;
  dfkdot_dq[2][7] = 0;
  dfkdot_dq[2][8] = 0;
  dfkdot_dq[3][0] = 0;
  dfkdot_dq[3][1] = 0;
  dfkdot_dq[3][2] = z[93];
  dfkdot_dq[3][3] = 0;
  dfkdot_dq[3][4] = 0;
  dfkdot_dq[3][5] = 0;
  dfkdot_dq[3][6] = 0;
  dfkdot_dq[3][7] = 0;
  dfkdot_dq[3][8] = 0;
  dfkdot_dq[4][0] = 0;
  dfkdot_dq[4][1] = 0;
  dfkdot_dq[4][2] = -z[93];
  dfkdot_dq[4][3] = 0;
  dfkdot_dq[4][4] = 0;
  dfkdot_dq[4][5] = 0;
  dfkdot_dq[4][6] = 0;
  dfkdot_dq[4][7] = 0;
  dfkdot_dq[4][8] = 0;
  dfkdot_dq[5][0] = 0;
  dfkdot_dq[5][1] = 0;
  dfkdot_dq[5][2] = -z[94];
  dfkdot_dq[5][3] = 0;
  dfkdot_dq[5][4] = 0;
  dfkdot_dq[5][5] = 0;
  dfkdot_dq[5][6] = 0;
  dfkdot_dq[5][7] = 0;
  dfkdot_dq[5][8] = 0;
  dfkdot_dq[6][0] = 0;
  dfkdot_dq[6][1] = 0;
  dfkdot_dq[6][2] = 0;
  dfkdot_dq[6][3] = 0;
  dfkdot_dq[6][4] = 0;
  dfkdot_dq[6][5] = 0;
  dfkdot_dq[6][6] = 0;
  dfkdot_dq[6][7] = 0;
  dfkdot_dq[6][8] = 0;
  dfkdot_dq[7][0] = 0;
  dfkdot_dq[7][1] = 0;
  dfkdot_dq[7][2] = 0;
  dfkdot_dq[7][3] = 0;
  dfkdot_dq[7][4] = 0;
  dfkdot_dq[7][5] = 0;
  dfkdot_dq[7][6] = 0;
  dfkdot_dq[7][7] = 0;
  dfkdot_dq[7][8] = 0;
  dfkdot_dq[8][0] = 0;
  dfkdot_dq[8][1] = 0;
  dfkdot_dq[8][2] = z[96];
  dfkdot_dq[8][3] = z[96];
  dfkdot_dq[8][4] = 0;
  dfkdot_dq[8][5] = 0;
  dfkdot_dq[8][6] = 0;
  dfkdot_dq[8][7] = 0;
  dfkdot_dq[8][8] = 0;
  dfkdot_dq[9][0] = 0;
  dfkdot_dq[9][1] = 0;
  dfkdot_dq[9][2] = z[197];
  dfkdot_dq[9][3] = z[197];
  dfkdot_dq[9][4] = 0;
  dfkdot_dq[9][5] = 0;
  dfkdot_dq[9][6] = 0;
  dfkdot_dq[9][7] = 0;
  dfkdot_dq[9][8] = 0;
  dfkdot_dq[10][0] = 0;
  dfkdot_dq[10][1] = 0;
  dfkdot_dq[10][2] = z[95];
  dfkdot_dq[10][3] = z[95];
  dfkdot_dq[10][4] = 0;
  dfkdot_dq[10][5] = 0;
  dfkdot_dq[10][6] = 0;
  dfkdot_dq[10][7] = 0;
  dfkdot_dq[10][8] = 0;
  dfkdot_dq[11][0] = 0;
  dfkdot_dq[11][1] = 0;
  dfkdot_dq[11][2] = z[96];
  dfkdot_dq[11][3] = z[96];
  dfkdot_dq[11][4] = 0;
  dfkdot_dq[11][5] = 0;
  dfkdot_dq[11][6] = 0;
  dfkdot_dq[11][7] = 0;
  dfkdot_dq[11][8] = 0;
  dfkdot_dq[12][0] = 0;
  dfkdot_dq[12][1] = 0;
  dfkdot_dq[12][2] = -z[125];
  dfkdot_dq[12][3] = -z[125];
  dfkdot_dq[12][4] = 0;
  dfkdot_dq[12][5] = 0;
  dfkdot_dq[12][6] = 0;
  dfkdot_dq[12][7] = 0;
  dfkdot_dq[12][8] = 0;
  dfkdot_dq[13][0] = 0;
  dfkdot_dq[13][1] = 0;
  dfkdot_dq[13][2] = -z[126];
  dfkdot_dq[13][3] = -z[126];
  dfkdot_dq[13][4] = 0;
  dfkdot_dq[13][5] = 0;
  dfkdot_dq[13][6] = 0;
  dfkdot_dq[13][7] = 0;
  dfkdot_dq[13][8] = 0;
  dfkdot_dq[14][0] = 0;
  dfkdot_dq[14][1] = 0;
  dfkdot_dq[14][2] = z[128];
  dfkdot_dq[14][3] = z[128];
  dfkdot_dq[14][4] = z[101];
  dfkdot_dq[14][5] = 0;
  dfkdot_dq[14][6] = 0;
  dfkdot_dq[14][7] = 0;
  dfkdot_dq[14][8] = 0;
  dfkdot_dq[15][0] = 0;
  dfkdot_dq[15][1] = 0;
  dfkdot_dq[15][2] = z[131];
  dfkdot_dq[15][3] = z[131];
  dfkdot_dq[15][4] = z[133];
  dfkdot_dq[15][5] = 0;
  dfkdot_dq[15][6] = 0;
  dfkdot_dq[15][7] = 0;
  dfkdot_dq[15][8] = 0;
  dfkdot_dq[16][0] = 0;
  dfkdot_dq[16][1] = 0;
  dfkdot_dq[16][2] = z[100];
  dfkdot_dq[16][3] = z[100];
  dfkdot_dq[16][4] = z[103];
  dfkdot_dq[16][5] = 0;
  dfkdot_dq[16][6] = 0;
  dfkdot_dq[16][7] = 0;
  dfkdot_dq[16][8] = 0;
  dfkdot_dq[17][0] = 0;
  dfkdot_dq[17][1] = 0;
  dfkdot_dq[17][2] = z[101];
  dfkdot_dq[17][3] = z[101];
  dfkdot_dq[17][4] = z[135];
  dfkdot_dq[17][5] = 0;
  dfkdot_dq[17][6] = 0;
  dfkdot_dq[17][7] = 0;
  dfkdot_dq[17][8] = 0;
  dfkdot_dq[18][0] = 0;
  dfkdot_dq[18][1] = 0;
  dfkdot_dq[18][2] = z[136];
  dfkdot_dq[18][3] = z[136];
  dfkdot_dq[18][4] = -z[137];
  dfkdot_dq[18][5] = 0;
  dfkdot_dq[18][6] = 0;
  dfkdot_dq[18][7] = 0;
  dfkdot_dq[18][8] = 0;
  dfkdot_dq[19][0] = 0;
  dfkdot_dq[19][1] = 0;
  dfkdot_dq[19][2] = z[138];
  dfkdot_dq[19][3] = z[138];
  dfkdot_dq[19][4] = -z[139];
  dfkdot_dq[19][5] = 0;
  dfkdot_dq[19][6] = 0;
  dfkdot_dq[19][7] = 0;
  dfkdot_dq[19][8] = 0;
  dfkdot_dq[20][0] = 0;
  dfkdot_dq[20][1] = 0;
  dfkdot_dq[20][2] = z[142];
  dfkdot_dq[20][3] = z[142];
  dfkdot_dq[20][4] = z[144];
  dfkdot_dq[20][5] = z[107];
  dfkdot_dq[20][6] = 0;
  dfkdot_dq[20][7] = 0;
  dfkdot_dq[20][8] = 0;
  dfkdot_dq[21][0] = 0;
  dfkdot_dq[21][1] = 0;
  dfkdot_dq[21][2] = z[148];
  dfkdot_dq[21][3] = z[148];
  dfkdot_dq[21][4] = z[151];
  dfkdot_dq[21][5] = z[153];
  dfkdot_dq[21][6] = 0;
  dfkdot_dq[21][7] = 0;
  dfkdot_dq[21][8] = 0;
  dfkdot_dq[22][0] = 0;
  dfkdot_dq[22][1] = 0;
  dfkdot_dq[22][2] = z[106];
  dfkdot_dq[22][3] = z[106];
  dfkdot_dq[22][4] = z[155];
  dfkdot_dq[22][5] = z[109];
  dfkdot_dq[22][6] = 0;
  dfkdot_dq[22][7] = 0;
  dfkdot_dq[22][8] = 0;
  dfkdot_dq[23][0] = 0;
  dfkdot_dq[23][1] = 0;
  dfkdot_dq[23][2] = z[107];
  dfkdot_dq[23][3] = z[107];
  dfkdot_dq[23][4] = z[158];
  dfkdot_dq[23][5] = z[160];
  dfkdot_dq[23][6] = 0;
  dfkdot_dq[23][7] = 0;
  dfkdot_dq[23][8] = 0;
  dfkdot_dq[24][0] = 0;
  dfkdot_dq[24][1] = 0;
  dfkdot_dq[24][2] = 0;
  dfkdot_dq[24][3] = 0;
  dfkdot_dq[24][4] = 0;
  dfkdot_dq[24][5] = 0;
  dfkdot_dq[24][6] = 0;
  dfkdot_dq[24][7] = 0;
  dfkdot_dq[24][8] = 0;
  dfkdot_dq[25][0] = 0;
  dfkdot_dq[25][1] = 0;
  dfkdot_dq[25][2] = 0;
  dfkdot_dq[25][3] = 0;
  dfkdot_dq[25][4] = 0;
  dfkdot_dq[25][5] = 0;
  dfkdot_dq[25][6] = 0;
  dfkdot_dq[25][7] = 0;
  dfkdot_dq[25][8] = 0;
  dfkdot_dq[26][0] = 0;
  dfkdot_dq[26][1] = 0;
  dfkdot_dq[26][2] = z[111];
  dfkdot_dq[26][3] = 0;
  dfkdot_dq[26][4] = 0;
  dfkdot_dq[26][5] = 0;
  dfkdot_dq[26][6] = z[111];
  dfkdot_dq[26][7] = 0;
  dfkdot_dq[26][8] = 0;
  dfkdot_dq[27][0] = 0;
  dfkdot_dq[27][1] = 0;
  dfkdot_dq[27][2] = z[198];
  dfkdot_dq[27][3] = 0;
  dfkdot_dq[27][4] = 0;
  dfkdot_dq[27][5] = 0;
  dfkdot_dq[27][6] = z[198];
  dfkdot_dq[27][7] = 0;
  dfkdot_dq[27][8] = 0;
  dfkdot_dq[28][0] = 0;
  dfkdot_dq[28][1] = 0;
  dfkdot_dq[28][2] = z[110];
  dfkdot_dq[28][3] = 0;
  dfkdot_dq[28][4] = 0;
  dfkdot_dq[28][5] = 0;
  dfkdot_dq[28][6] = z[110];
  dfkdot_dq[28][7] = 0;
  dfkdot_dq[28][8] = 0;
  dfkdot_dq[29][0] = 0;
  dfkdot_dq[29][1] = 0;
  dfkdot_dq[29][2] = z[111];
  dfkdot_dq[29][3] = 0;
  dfkdot_dq[29][4] = 0;
  dfkdot_dq[29][5] = 0;
  dfkdot_dq[29][6] = z[111];
  dfkdot_dq[29][7] = 0;
  dfkdot_dq[29][8] = 0;
  dfkdot_dq[30][0] = 0;
  dfkdot_dq[30][1] = 0;
  dfkdot_dq[30][2] = -z[161];
  dfkdot_dq[30][3] = 0;
  dfkdot_dq[30][4] = 0;
  dfkdot_dq[30][5] = 0;
  dfkdot_dq[30][6] = -z[161];
  dfkdot_dq[30][7] = 0;
  dfkdot_dq[30][8] = 0;
  dfkdot_dq[31][0] = 0;
  dfkdot_dq[31][1] = 0;
  dfkdot_dq[31][2] = -z[162];
  dfkdot_dq[31][3] = 0;
  dfkdot_dq[31][4] = 0;
  dfkdot_dq[31][5] = 0;
  dfkdot_dq[31][6] = -z[162];
  dfkdot_dq[31][7] = 0;
  dfkdot_dq[31][8] = 0;
  dfkdot_dq[32][0] = 0;
  dfkdot_dq[32][1] = 0;
  dfkdot_dq[32][2] = z[164];
  dfkdot_dq[32][3] = 0;
  dfkdot_dq[32][4] = 0;
  dfkdot_dq[32][5] = 0;
  dfkdot_dq[32][6] = z[164];
  dfkdot_dq[32][7] = z[116];
  dfkdot_dq[32][8] = 0;
  dfkdot_dq[33][0] = 0;
  dfkdot_dq[33][1] = 0;
  dfkdot_dq[33][2] = z[167];
  dfkdot_dq[33][3] = 0;
  dfkdot_dq[33][4] = 0;
  dfkdot_dq[33][5] = 0;
  dfkdot_dq[33][6] = z[167];
  dfkdot_dq[33][7] = z[169];
  dfkdot_dq[33][8] = 0;
  dfkdot_dq[34][0] = 0;
  dfkdot_dq[34][1] = 0;
  dfkdot_dq[34][2] = z[115];
  dfkdot_dq[34][3] = 0;
  dfkdot_dq[34][4] = 0;
  dfkdot_dq[34][5] = 0;
  dfkdot_dq[34][6] = z[115];
  dfkdot_dq[34][7] = z[118];
  dfkdot_dq[34][8] = 0;
  dfkdot_dq[35][0] = 0;
  dfkdot_dq[35][1] = 0;
  dfkdot_dq[35][2] = z[116];
  dfkdot_dq[35][3] = 0;
  dfkdot_dq[35][4] = 0;
  dfkdot_dq[35][5] = 0;
  dfkdot_dq[35][6] = z[116];
  dfkdot_dq[35][7] = z[171];
  dfkdot_dq[35][8] = 0;
  dfkdot_dq[36][0] = 0;
  dfkdot_dq[36][1] = 0;
  dfkdot_dq[36][2] = z[172];
  dfkdot_dq[36][3] = 0;
  dfkdot_dq[36][4] = 0;
  dfkdot_dq[36][5] = 0;
  dfkdot_dq[36][6] = z[172];
  dfkdot_dq[36][7] = -z[173];
  dfkdot_dq[36][8] = 0;
  dfkdot_dq[37][0] = 0;
  dfkdot_dq[37][1] = 0;
  dfkdot_dq[37][2] = z[174];
  dfkdot_dq[37][3] = 0;
  dfkdot_dq[37][4] = 0;
  dfkdot_dq[37][5] = 0;
  dfkdot_dq[37][6] = z[174];
  dfkdot_dq[37][7] = -z[175];
  dfkdot_dq[37][8] = 0;
  dfkdot_dq[38][0] = 0;
  dfkdot_dq[38][1] = 0;
  dfkdot_dq[38][2] = z[178];
  dfkdot_dq[38][3] = 0;
  dfkdot_dq[38][4] = 0;
  dfkdot_dq[38][5] = 0;
  dfkdot_dq[38][6] = z[178];
  dfkdot_dq[38][7] = z[180];
  dfkdot_dq[38][8] = z[122];
  dfkdot_dq[39][0] = 0;
  dfkdot_dq[39][1] = 0;
  dfkdot_dq[39][2] = z[184];
  dfkdot_dq[39][3] = 0;
  dfkdot_dq[39][4] = 0;
  dfkdot_dq[39][5] = 0;
  dfkdot_dq[39][6] = z[184];
  dfkdot_dq[39][7] = z[187];
  dfkdot_dq[39][8] = z[189];
  dfkdot_dq[40][0] = 0;
  dfkdot_dq[40][1] = 0;
  dfkdot_dq[40][2] = z[121];
  dfkdot_dq[40][3] = 0;
  dfkdot_dq[40][4] = 0;
  dfkdot_dq[40][5] = 0;
  dfkdot_dq[40][6] = z[121];
  dfkdot_dq[40][7] = z[191];
  dfkdot_dq[40][8] = z[124];
  dfkdot_dq[41][0] = 0;
  dfkdot_dq[41][1] = 0;
  dfkdot_dq[41][2] = z[122];
  dfkdot_dq[41][3] = 0;
  dfkdot_dq[41][4] = 0;
  dfkdot_dq[41][5] = 0;
  dfkdot_dq[41][6] = z[122];
  dfkdot_dq[41][7] = z[194];
  dfkdot_dq[41][8] = z[196];

}
