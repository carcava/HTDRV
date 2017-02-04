#ifndef __COMM_H
#define __COMM_H

#ifdef __USE_MPI
#include "mpi.h"
#else
typedef int MPI_Comm;
const int MPI_COMM_WORLD=0;
const int MPI_COMM_NULL = 0;
#endif

#ifdef __USE_MPI
const bool Parallel = true;
#else
const bool Parallel = false;
#endif

#define ROOT_PE 0
#define NDIGIT_IN_RANK 6

const int SHUTDOWN_MESSAGE = -1;
const int SHUTDOWN_MESSAGE_SIZE = 1;

class CommGroup {
	MPI_Comm comm;
	int num_pe;
	int my_pe;
	int root_pe;
	bool I_am_root;
public:
	CommGroup() {
		comm = MPI_COMM_NULL;
		num_pe = 1;
		my_pe = 0;
		root_pe = 0;
		I_am_root = true;
	}
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
	void Sync() const {
#ifdef __USE_MPI
		MPI_Barrier(comm);
#endif
	}
	int FortranComm() const {
		return MPI_Comm_c2f(comm);
	}
};


class MasterSlave {
	CommGroup parent;
	CommGroup slaves;
	int group_max_size;
	int number_of_slave_groups;
	int * group_leader;
public:
	MasterSlave() {
		number_of_slave_groups = 0;
		group_leader = 0;
	}
	MasterSlave(CommGroup parent_comm, int group_max_size); 
	~MasterSlave() {
		if (group_leader) {
			delete[] group_leader;
		}
		group_leader = 0;
		number_of_slave_groups = 0;
	}
	bool IamSlaveRoot() const {
		return slaves.IamRoot();
	}
	int SlaveRootPe() const {
		return slaves.RootPe();
	}
	MPI_Comm SlaveComm() const {
		return slaves.Comm();
	}
	const CommGroup * World() const {
		return &parent;
	}
	const CommGroup * Slaves() const {
		return &slaves;
	}
	int NumberOfSlaveGroups() const {
		return number_of_slave_groups;
	}
	int GroupLeader(int ind) const {
		return group_leader[ind];
	}
	int MyPe() const {
		return parent.MyPe();
	}
	int SayIamReady(int * message, int n) const;
	int BcastToOtherSlaves(int * message, int n) const;
	int SayWeAreReady(int * message, int n)const;
	int GetWorkFromMaster(int * message, int n)const;
	int ShutDownSlaves() const;
	int GetReadySlave(int * message, int n)const;
	int SendWorkToSlave(int slave_id, int * message, int n )const;
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
