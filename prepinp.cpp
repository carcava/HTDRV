#include<stdio.h>
#include<stdlib.h>
#include<list>

#if defined(_WIN32)
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <string.h>
#include "comm.h"

#ifdef __cplusplus
extern "C" {
#endif

/* interface to pw.x */
/* launch a pw.x-like calculation */
void ht_pw_drv(int lib_comm, int nimage, int npot, int npool, int ntaskgroup,
                int nband, int ndiag, int *exit_status, char *input_file);

int c_mkdir_safe( const char * dirname );

double cclock_();

#ifdef __cplusplus
}
#endif


const int NOXYGEN=31;
const int NHYDROGEN=62;
const int SHUTDOWN_MESSAGE=-1;
const int SHUTDOWN_MESSAGE_SIZE=1;

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

using namespace std;

class Pos {
	int stepid;
	double O[ NOXYGEN ][3];
	double H[ NHYDROGEN ][3];
public:
	Pos(){
		stepid=0;
		for(int i=0;i<NOXYGEN;i++)   O[i][0]=O[i][1]=O[i][2]=0.0;
		for(int i=0;i<NHYDROGEN;i++) H[i][0]=H[i][1]=H[i][2]=0.0;
	};
	void Readin( FILE * fp ) {
		char line[1024];
		fgets(line, 1023, fp);
		sscanf(line,"%d",&stepid);
		for(int i=0;i<NOXYGEN;i++) {
			fgets(line, 1023, fp);
			sscanf(line,"%lf %lf %lf",&O[i][0],&O[i][1],&O[i][2]);
		}
		for(int i=0;i<NHYDROGEN;i++) {
			fgets(line, 1023, fp);
			sscanf(line,"%lf %lf %lf",&H[i][0],&H[i][1],&H[i][2]);
		}
	};
	int StepID() const { return stepid; };
	void PrintOut( FILE *fp ) const { 
		for(int i=0;i<NOXYGEN;i++) {
			fprintf(fp,"O %lf %lf %lf\n",O[i][0],O[i][1],O[i][2]);
		}
		for(int i=0;i<NHYDROGEN;i++) {
			fprintf(fp,"H %lf %lf %lf\n",H[i][0],H[i][1],H[i][2]);
		}
		return; 
	};
};

class Cel {
	int stepid;
	double h[3][3];
public:
	Cel(){
		stepid=0;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) h[i][j]=0.0;
	};
	void Readin( FILE * fp ) {
		char line[1024];
		fgets(line, 1023, fp);
		sscanf(line,"%*s %d",&stepid);
		for(int i=0;i<3;i++) {
			fgets(line, 1023, fp);
			sscanf(line,"%lf %lf %lf",&h[i][0],&h[i][1],&h[i][2]);
		}
	};
	int StepID() const { return stepid; };
	void PrintOut( FILE *fp ) const { 
		for(int i=0;i<3;i++) {
			fprintf(fp,"%lf %lf %lf\n",h[i][0],h[i][1],h[i][2]);
		}
		return; 
	};
};

int MyChangeDir( const char * dirname ) {
#ifdef _WIN32
        return _chdir( dirname );
#else
        return chdir( dirname );
#endif
}

char * MyGetPWD() {
#ifdef _WIN32
        return _getcwd(0,0);
#else
        return getcwd(0,0);
#endif
}



void PrintSingleInput( const Pos & p, const Cel & c,  MPI_Comm this_slave_group) {
	char filename[80];
	char dirname[80];
	char currpath[1024];
	if( p.StepID() != c.StepID() ) {
		fprintf( stderr, "WARNING stepid do not match %d %d\n", p.StepID(), c.StepID() );
		return;
	}
	sprintf( filename, "h2o.pw.%d", p.StepID() );
	sprintf( dirname, "RUN%d", p.StepID() );
      	c_mkdir_safe( dirname ); 

	MyChangeDir( dirname );	

	FILE * fp = fopen(filename, "w");
	if( !fp ) {
		fprintf(stderr,"ERROR opening file %s\n", filename );
		exit(1);
	}
	fprintf(fp,"%s",CONTROL_NAMELIST);
	fprintf(fp,"%s",SYSTEM_NAMELIST);
	fprintf(fp,"%s",ELECTRONS_NAMELIST);
	fprintf(fp,"%s",ATOMIC_SPECIES);
	fprintf(fp,"%s",K_POINTS);
	fprintf(fp,"%s",CELL_PARAMETERS);
	c.PrintOut(fp);
	fprintf(fp,"%s",ATOMIC_POSITIONS);
	p.PrintOut(fp);
	fclose(fp);

        int nimage = 1;
        int npots = 1;
        int npool = 1;
        int ntg = 1;
        int nband = 1;
        int ndiag = 1;
        int retval = 0;

        ht_pw_drv( this_slave_group, nimage, npots, npool, ntg, nband, ndiag, &retval, filename);

	MyChangeDir( "../" );	
	
}


void SlaveTask( list<Pos> & positions, list<Cel> & cells, const CommGroup & Group, const CommGroup & World, int * group_leader, int number_of_slave_groups ) {

	bool work_to_do = true;
	int message = 0;
	int message_size = 1;

	while( work_to_do ) {

		if( Group.IamRoot() ) {
			comm_say_i_am_ready( World, &message, message_size );
		}
		MPI_Bcast( &message, message_size, MPI_INT, Group.RootPe(), Group.Comm() );
		
		if( message > 0 ) {
			// send back data to master
		}
		if( Group.IamRoot() ) {
			comm_recv( &message, message_size, World.RootPe(), World);
		}
		MPI_Bcast( &message, message_size, MPI_INT, Group.RootPe(), Group.Comm() );

		if( message > -1 ) {
			list<Pos>::const_iterator i = positions.begin();
			list<Cel>::const_iterator j = cells.begin();
			while( (i != positions.end() ) && ( j != cells.end() ) ) {
				if( i->StepID() == message ) {
					fprintf(stdout,"Slave %d: processing stepid = %d\n", World.MyPe(), i->StepID() );
					PrintSingleInput( *i, *j, Group.Comm() );
				}
				i++;
				j++;
			}
		} else {
			work_to_do = false;
		}
	}

	fprintf(stdout, "SLAVE %d: shutdown\n", World.MyPe());
}


void MasterTask( list<Pos> & positions, list<Cel> & cells, const CommGroup & Group, const CommGroup & World, int * group_leader, int number_of_slave_groups ) {
	list<Pos>::const_iterator i = positions.begin();
	list<Cel>::const_iterator j = cells.begin();
	while( (i != positions.end() ) && ( j != cells.end() ) ) {
		fprintf(stdout,"stepid = %d %d\n", i->StepID(), j->StepID() );
		if( World.NumPe() < 2 ) {
			PrintSingleInput( *i, *j, Group.Comm() );
		} else {
			int message = 0;
			int message_size = 1;
			int slave_id = comm_get_ready_slave(World, &message, message_size );
			if( message > 0 ) {
				// post-process slave work ... if any
			}
			// send work to do
			message = i->StepID();
			fprintf(stdout, "MASTER: assigning to SLAVE %d Task %d\n", slave_id, message );
			comm_send( &message, message_size, slave_id, slave_id, World);
		}

		i++;
		j++;
	}
	for( int ip = 0; ip < number_of_slave_groups; ip++ ) {
		// send shut-down message
		fprintf(stdout, "MASTER: sending to SLAVE %d shutdown signal\n", ip);
		comm_send( &SHUTDOWN_MESSAGE, SHUTDOWN_MESSAGE_SIZE, group_leader[ ip ], ip, World);
		//
	}
}

void ReadPositions( list<Pos> & positions, int nstep ) {
	char pos_filename[]="h2o.pos";
	FILE * in = fopen(pos_filename, "r");
	if( !in ) {
		fprintf(stderr,"ERROR opening file %s\n", pos_filename );
		exit(1);
	}
	for( int i = 0; i < nstep; i++ ) {
		Pos tmp;
		tmp.Readin( in );
		positions.push_back( tmp );
	}
	fclose(in);
}

void ReadCels( list<Cel> & cells, int nstep ) {
	char cel_filename[]="h2o.cel";
	FILE * in = fopen(cel_filename, "r");
	if( !in ) {
		fprintf(stderr,"ERROR opening file %s\n", cel_filename );
		exit(1);
	}
	for( int i = 0; i < nstep; i++ ) {
		Cel tmp;
		tmp.Readin( in );
		cells.push_back( tmp );
	}
	fclose(in);
}



int main(int argc, char ** argv) {

        comm_begin( &argc, &argv );

        CommGroup World( MPI_COMM_WORLD );

        int j = 0;
	int nstep = 0;
	int group_size = 1;
        while ( ++j < argc ) {
                if (!strcmp(argv[j],"-nsys")) {
			// number of systems to be processed
			nstep = atoi( argv[++j] );
                } else if (!strcmp(argv[j],"-nproc")) {
			// number of systems to be processed
			group_size = atoi( argv[++j] );
                }
                else {
			fprintf(stderr,"wrong argument %s\n",argv[j]);
			exit(1);
                }
        }

	if( nstep <= 0 ) {
		fprintf(stderr,"wrong number of systems\n");
		exit(1);
	}

	int role;

	if ( World.IamRoot() ) {
		role = 0;
	} else {
		role = ( World.MyPe() - 1 ) / group_size + 1;
	}

	MPI_Comm intra_comm;
	MPI_Comm_split( MPI_COMM_WORLD, role, World.MyPe(), &intra_comm);
	CommGroup MyGroup( intra_comm );
	int my_fcomm = MPI_Comm_c2f( intra_comm );

	int number_of_slave_groups = ( World.NumPe() - 2 ) / group_size + 1;
	int * group_leader = new int[ number_of_slave_groups ];
	int * group_leader_mask = new int[ World.NumPe() ];

	if( World.IamRoot() )
		fprintf(stdout,"Number of slave groups = %d\n", number_of_slave_groups);

	for( int i = 0; i < World.NumPe(); i++ ) {
		group_leader_mask[ i ] = 0;
	}

	if( MyGroup.IamRoot() && !World.IamRoot() ) {
		group_leader_mask[ World.MyPe() ] = 1;
	}

	MPI_Allreduce(MPI_IN_PLACE, group_leader_mask, World.NumPe(), MPI_INT, MPI_SUM, World.Comm());
	int ind_leader = 0;
	for( int i = 0; i < World.NumPe(); i++ ) {
		if( group_leader_mask[ i ] ) 
			group_leader[ ind_leader++ ] = i;
	}

	if( World.IamRoot() ) {
		for( int i = 0; i < number_of_slave_groups; i++ ) {
			fprintf(stdout,"Group: %d, leader: %d\n", i, group_leader[i] );
		}
	}

	double time_begin = cclock_();

	list<Pos> positions;
	list<Cel> cells;

	ReadPositions(positions, nstep);
	ReadCels(cells,nstep);

	if( World.IamRoot() ) {
		MasterTask( positions, cells, MyGroup, World, group_leader, number_of_slave_groups );
	} else {
		SlaveTask( positions, cells, MyGroup, World, group_leader, number_of_slave_groups );
	}

	double time_end = cclock_();

	if( World.IamRoot() ) {
		fprintf( stdout, "WALLTIME: %lf seconds", time_end - time_begin);
	}

	delete [] group_leader;
	delete [] group_leader_mask;

	comm_end();

	return 0;
}


#if defined(__STANDALONE)
void ht_pw_drv(int lib_comm, int nimage, int npot, int npool, int ntaskgroup,
	int nband, int ndiag, int *exit_status, char *input_file) {
	fprintf(stdout, "DO SOME USEFUL WORK WITH: %s\n", input_file);
}
#endif