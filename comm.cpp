#include "comm.h"
#include "stdio.h"



MasterSlave::MasterSlave(CommGroup parent_comm, int group_max_size) {
	this->group_max_size = group_max_size;
	this->parent = parent_comm;

	int role;
	if (this->parent.IamRoot()) {
		role = 0;
	}
	else {
		role = (this->parent.MyPe() - 1) / this->group_max_size + 1;
	}

	MPI_Comm intra_comm;
	MPI_Comm_split( this->parent.Comm(), role, this->parent.MyPe(), &intra_comm );
	CommGroup _tmp( intra_comm );
	this->slaves = _tmp;

	this->number_of_slave_groups = (this->parent.NumPe() - 2) / this->group_max_size + 1;
	this->group_leader = new int[ this->number_of_slave_groups ];
	int * group_leader_mask = new int[ this->parent.NumPe() ];

	if (this->parent.IamRoot())
		fprintf(stdout, "Number of slave groups = %d\n", this->number_of_slave_groups);

	for (int i = 0; i < this->parent.NumPe(); i++) {
		group_leader_mask[i] = 0;
	}

	if (this->slaves.IamRoot() && !this->parent.IamRoot()) {
		group_leader_mask[this->parent.MyPe()] = 1;
	}

	MPI_Allreduce(MPI_IN_PLACE, group_leader_mask, this->parent.NumPe(), MPI_INT, MPI_SUM, this->parent.Comm());
	int ind_leader = 0;
	for (int i = 0; i < this->parent.NumPe(); i++) {
		if (group_leader_mask[i])
			this->group_leader[ind_leader++] = i;
	}

	if (this->parent.IamRoot()) {
		for (int i = 0; i < this->number_of_slave_groups; i++) {
			fprintf(stdout, "Group: %d, leader: %d\n", i, this->group_leader[i]);
		}
	}
	delete[] group_leader_mask;
}

void comm_bcast( int * val, int nval, const CommGroup & G ){
#ifdef __USE_MPI
	MPI_Bcast( (void *)val, nval, MPI_INT, ROOT_PE, G.Comm() );
#endif
}



void comm_send(const int * val, int nval, int dest, int tag, const CommGroup & G ){
#ifdef __USE_MPI
	MPI_Send( (const void *)val, nval, MPI_INT, dest, tag, G.Comm() );
#endif
}

void comm_send_string(const char * val, int nval, int dest, int tag, const CommGroup & G) {
#ifdef __USE_MPI
	MPI_Send((const void *)val, nval, MPI_CHAR, dest, tag, G.Comm());
#endif
}

void comm_recv( int * val, int nval, int sour, const CommGroup & G ) {
#ifdef __USE_MPI
	MPI_Status status;
	MPI_Recv( (void *)val, nval, MPI_INT, sour, MPI_ANY_TAG, G.Comm(), &status );
#endif
}

void comm_recv_string(char * val, int nval, int sour, const CommGroup & G) {
#ifdef __USE_MPI
	MPI_Status status;
	MPI_Recv((void *)val, nval, MPI_CHAR, sour, MPI_ANY_TAG, G.Comm(), &status);
#endif
}

void comm_recv( int * val, int nval, const CommGroup & G ) {
#ifdef __USE_MPI
	MPI_Status status;
	MPI_Recv( (void *)val, nval, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, G.Comm(), &status );
#endif
}

void comm_rootgather(long * gathered_val, long val, const CommGroup & G) {
#ifdef __USE_MPI
	MPI_Gather(&val, 1, MPI_LONG, gathered_val, 1, MPI_LONG, G.RootPe(), G.Comm());
#endif
}

int comm_get_ready_slave( const CommGroup & G, int * message, int n ) {
	int whois = 0;
#ifdef __USE_MPI
	MPI_Status status;
	MPI_Recv( (void *) message, n, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, G.Comm(), &status );
	whois = status.MPI_SOURCE;
#endif
	return whois;
}

int comm_say_i_am_ready( const CommGroup & G, int * message, int n ) {
	comm_send( message, n, G.RootPe(), G.MyPe(), G );
	return G.MyPe();
}
