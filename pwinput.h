#ifndef _PWINPUT_H_
#define _PWINPUT_H_



const int NOXYGEN = 31;
const int NHYDROGEN = 62;




using namespace std;

class Pos {
	int stepid;
	int no;
	int nh;
	double O[NOXYGEN][3];
	double H[NHYDROGEN][3];
public:
	Pos();
	void Readin(FILE * fp);
	int StepID() const { return stepid; };
	void PrintOut(FILE *fp) const;
};

class Cel {
	int stepid;
	double h[3][3];
public:
	Cel();
	void Readin(FILE * fp);
	int StepID() const { return stepid; };
	void PrintOut(FILE *fp) const;
};


void WritePwInput(const char * filename, const Pos & p, const Cel & c);

class PwInput {

};

#endif
