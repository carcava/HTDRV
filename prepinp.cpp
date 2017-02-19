#include<stdio.h>
#include<stdlib.h>
#include<list>
#include"pwinput.h"

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



void RunPw( const Pos & p, const Cel & c, const CommGroup * Group) {
	char filename[80];
	char dirname[80];

	if( p.StepID() != c.StepID() ) {
		fprintf( stderr, "WARNING stepid do not match %d %d\n", p.StepID(), c.StepID() );
		return;
	}
	sprintf( filename, "h2o.pw.%d", p.StepID() );
	sprintf( dirname, "RUN%d", p.StepID() );

	if( Group->IamRoot())
		c_mkdir_safe( dirname );
	Group->Sync();

	MyChangeDir( dirname );	

	if (Group->IamRoot()) {
		WritePwInput(filename, p, c);
	}
	Group->Sync();

        int nimage = 1;
        int npots = 1;
        int npool = 1;
        int ntg = 1;
        int nband = 1;
        int ndiag = 1;
        int retval = 0;

        ht_pw_drv( Group->FortranComm(), nimage, npots, npool, ntg, nband, ndiag, &retval, filename);

	MyChangeDir( "../" );	
	
}


void SlaveTask( list<Pos> & positions, list<Cel> & cells, const MasterSlave & MS) {

	bool work_to_do = true;
	int message = 0;
	int message_size = 1;

	while( work_to_do ) {

		MS.SayWeAreReady(&message, message_size);
		
		if( message > 0 ) {
			// send back data to master
		}

		MS.GetWorkFromMaster(&message, message_size);

		if( message > -1 ) {
			// DO WORK
			list<Pos>::const_iterator i = positions.begin();
			list<Cel>::const_iterator j = cells.begin();
			while( (i != positions.end() ) && ( j != cells.end() ) ) {
				if( i->StepID() == message ) {
					fprintf(stdout,"Slave %d: processing stepid = %d\n", MS.MyPe(), i->StepID() );
					RunPw( *i, *j, MS.Slaves() );
				}
				i++;
				j++;
			}
		} else {
			work_to_do = false;
		}
	}

	fprintf(stdout, "SLAVE %d: shutdown\n", MS.MyPe());
}


void MasterTask( list<Pos> & positions, list<Cel> & cells, const MasterSlave & MS ) {
	list<Pos>::const_iterator i = positions.begin();
	list<Cel>::const_iterator j = cells.begin();
	while( (i != positions.end() ) && ( j != cells.end() ) ) {
		fprintf(stdout,"stepid = %d %d\n", i->StepID(), j->StepID() );
		if(MS.World()->NumPe() < 2 ) {
			RunPw( *i, *j, MS.Slaves());
		} else {
			int message = 0;
			int message_size = 1;
			int slave_id = MS.GetReadySlave( &message, message_size );
			if( message > 0 ) {
				// post-process slave work ... if any
			}
			// send work to do
			message = i->StepID();
			fprintf(stdout, "MASTER: assigning to SLAVE %d Task %d\n", slave_id, message );
			MS.SendWorkToSlave( slave_id, &message, message_size );
		}

		i++;
		j++;
	}
	MS.ShutDownSlaves();
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
			// numbero of proocessor for each slave group
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

	MasterSlave MS(World, group_size);

	double time_begin = cclock_();

	list<Pos> positions;
	list<Cel> cells;

	ReadPositions(positions, nstep);
	ReadCels(cells,nstep);

	if( World.IamRoot() ) {
		MasterTask( positions, cells, MS );
	} else {
		SlaveTask( positions, cells, MS );
	}

	double time_end = cclock_();

	if( World.IamRoot() ) {
		fprintf( stdout, "WALLTIME: %lf seconds", time_end - time_begin);
	}

	comm_end();

	return 0;
}


#if defined(__STANDALONE)
void ht_pw_drv(int lib_comm, int nimage, int npot, int npool, int ntaskgroup,
	int nband, int ndiag, int *exit_status, char *input_file) {
	fprintf(stdout, "DO SOME USEFUL WORK WITH: %s\n", input_file);
}
#endif