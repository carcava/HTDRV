#include<stdio.h>
#include<stdlib.h>
#include "pwinput.h"

const char CONTROL_NAMELIST[] =
"&control\n"
"    calculation = 'scf'\n"
"    restart_mode='from_scratch'\n"
"    prefix='H2O-60',\n"
"    tprnfor = .true.\n"
"    tstress = .true.\n"
"    pseudo_dir = '../'\n"
"    outdir = './'\n"
"/\n";

const char SYSTEM_NAMELIST[] =
" &system\n"
"      ibrav=0,\n"
"      nat=93,\n"
"      ntyp=2,\n"
"      ecutwfc=80.0\n"
"      ecfixed=70.0\n"
"      qcutz=70.0\n"
"      q2sigma=5.0\n"
"    occupations='smearing',\n"
"    smearing='gaussian',\n"
"    degauss=0.10,\n"
"/\n";

const char ELECTRONS_NAMELIST[] =
"&electrons\n"
"    electron_maxstep = 50,\n"
"    mixing_mode = 'plain'\n"
"    mixing_beta = 0.3\n"
"    conv_thr =  1.0d-06\n"
"/\n";

const char ATOMIC_SPECIES[] =
"ATOMIC_SPECIES\n"
" O  16.0  O.BLYP.UPF\n"
" H   1.0  H.fpmd.UPF\n";

const char K_POINTS[] =
"K_POINTS {gamma}\n"
"1 1 1   0 0 0\n";

const char CELL_PARAMETERS[] =
"CELL_PARAMETERS (bohr)\n";

const char ATOMIC_POSITIONS[] =
"ATOMIC_POSITIONS (bohr)\n";

Pos::Pos() {
	stepid = 0;
	no = NOXYGEN;
	nh = NHYDROGEN;
	for (int i = 0; i<no; i++)   O[i][0] = O[i][1] = O[i][2] = 0.0;
	for (int i = 0; i<nh; i++) H[i][0] = H[i][1] = H[i][2] = 0.0;
};
void Pos::Readin(FILE * fp) {
	char line[1024];
	fgets(line, 1023, fp);
	sscanf(line, "%d", &stepid);
	for (int i = 0; i<no; i++) {
		fgets(line, 1023, fp);
		sscanf(line, "%lf %lf %lf", &O[i][0], &O[i][1], &O[i][2]);
	}
	for (int i = 0; i<nh; i++) {
		fgets(line, 1023, fp);
		sscanf(line, "%lf %lf %lf", &H[i][0], &H[i][1], &H[i][2]);
	}
};

void Pos::PrintOut(FILE *fp) const {
	for (int i = 0; i<no; i++) {
		fprintf(fp, "O %lf %lf %lf\n", O[i][0], O[i][1], O[i][2]);
	}
	for (int i = 0; i<nh; i++) {
		fprintf(fp, "H %lf %lf %lf\n", H[i][0], H[i][1], H[i][2]);
	}
	return;
};

Cel::Cel() {
	stepid = 0;
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++) h[i][j] = 0.0;
};
void Cel::Readin(FILE * fp) {
	char line[1024];
	fgets(line, 1023, fp);
	sscanf(line, "%*s %d", &stepid);
	for (int i = 0; i<3; i++) {
		fgets(line, 1023, fp);
		sscanf(line, "%lf %lf %lf", &h[i][0], &h[i][1], &h[i][2]);
	}
};

void Cel::PrintOut(FILE *fp) const {
	for (int i = 0; i<3; i++) {
		fprintf(fp, "%lf %lf %lf\n", h[i][0], h[i][1], h[i][2]);
	}
	return;
};


void WritePwInput( const char * filename, const Pos & p, const Cel & c ) {
	FILE * fp = fopen(filename, "w");
	if (!fp) {
		fprintf(stderr, "ERROR opening file %s\n", filename);
		exit(1);
	}
	fprintf(fp, "%s", CONTROL_NAMELIST);
	fprintf(fp, "%s", SYSTEM_NAMELIST);
	fprintf(fp, "%s", ELECTRONS_NAMELIST);
	fprintf(fp, "%s", ATOMIC_SPECIES);
	fprintf(fp, "%s", K_POINTS);
	fprintf(fp, "%s", CELL_PARAMETERS);
	c.PrintOut(fp);
	fprintf(fp, "%s", ATOMIC_POSITIONS);
	p.PrintOut(fp);
	fclose(fp);
}