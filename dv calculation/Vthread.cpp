#include<iostream>
#include<mpi.h>
//child process 

int main(int argc, char** argv) {
	int rank,*rec_row,rec_row_N,dv=0;
	//rec_row is pointer to the received row, have to allocate space for the received row as well
	MPI_Comm parentcomm;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status status;
	//MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_get_parent(&parentcomm);
	if (parentcomm==MPI_COMM_NULL)
		std::cout<<"ERROR!, No parent for the child"<<std::endl;
	else
		{
		//Don't know the size of row to be received
		//Two ways to implement: 1. send the size from the parent process first
		// 2. use MPI_Probe because Only one to send messages to this process is the parent
		MPI_Probe(0,0,parentcomm,&status);
		MPI_Get_count(&status,MPI_INT,&rec_row_N);
		rec_row=new int[rec_row_N];		
		MPI_Recv(rec_row,rec_row_N,MPI_INT,0,0,parentcomm,MPI_STATUS_IGNORE);
		std::cout<<"I am the child number "<<rank<<" who received a row "<<std::endl;
		//sequential calculation of dv (it maybe optimised)
		for (int i=0;i<rec_row_N;++i)
			dv+=rec_row[i];
		//sending the dv data to parent process
		MPI_Send(&dv,1,MPI_INT,0,2,parentcomm);
		}	
	
	MPI_Finalize();
	return 0;
}
