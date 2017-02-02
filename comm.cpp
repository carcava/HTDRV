#include "comm.h"


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
