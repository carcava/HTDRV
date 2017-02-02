#ifndef __COMM_H
#define __COMM_H

#ifdef __USE_MPI
#include "mpi.h"
#else
typedef int MPI_Comm;
const int MPI_COMM_WORLD=0;
#endif

#ifdef __USE_MPI
const bool Parallel = true;
#else
const bool Parallel = false;
#endif

#define ROOT_PE 0
#define NDIGIT_IN_RANK 6

class CommGroup {
	MPI_Comm comm;
	int num_pe;
	int my_pe;
	int root_pe;
	bool I_am_root;
public:
	CommGroup( MPI_Comm comm ) {
		this->comm = comm;
		num_pe = 1;
		my_pe = 0;
		I_am_root = true;
#ifdef __USE_MPI
		MPI_Comm_size( comm, &num_pe );
		MPI_Comm_rank( comm, &my_pe );
#endif
		root_pe = ROOT_PE;
		if (my_pe == root_pe) I_am_root = true;
		else I_am_root = false;
	}

	int NumPe() const {
		return num_pe;
	}
	int MyPe() const {
		return my_pe;
	}
	int RootPe() const {
		return root_pe;
	}
	MPI_Comm Comm()const  {
		return comm;
	}
	bool IamRoot() const {
		return I_am_root;
	}
};

inline void comm_begin( int * argc, char *** argv ) {
#ifdef __USE_MPI
	MPI_Init( argc,  argv);
#endif
}

inline void comm_env( int * NumPe, int * MyPe, int * RootPe ) {
#ifdef __USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, NumPe);
	MPI_Comm_rank(MPI_COMM_WORLD, MyPe);
	*RootPe = ROOT_PE;
#endif
}

inline int comm_rank( ) {
	int MyPe = 0;
#ifdef __USE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &MyPe);
#endif
	return MyPe;
}


inline void comm_end() {
#ifdef __USE_MPI
	MPI_Finalize();
#endif
}

void comm_bcast( int * val, int nval, const CommGroup & G );
void comm_send( const int * val, int nval, int dest, int tag, const CommGroup & G );
void comm_send_string(const char * val, int nval, int dest, int tag, const CommGroup & G);
void comm_recv( int * val, int nval, int sour, const CommGroup & G );
void comm_recv_string(char * val, int nval, int sour, const CommGroup & G);
int comm_get_ready_slave( const CommGroup & G, int * message, int n );
int comm_say_i_am_ready( const CommGroup & G, int * message, int n );
void comm_rootgather(long * gathered_val, long val, const CommGroup & G);

#endif
