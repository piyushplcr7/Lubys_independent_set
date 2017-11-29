#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int main()
{
	FILE*file;
	file = fopen("./out.txt","r");
	int vertix = 10;

	int* matrix = malloc(vertix*vertix*sizeof(int));
	int row=1;
	int column;

	while(fscanf(file,"%d", &row)==1)
		{
			fscanf(file,"%d", &column);
			printf("row:%d, column:%d\n",row,column );
			matrix[row*vertix+column] = 1;
		}

	
	for (int i = 0; i < vertix; ++i)
	{
		for (int j = 0; j < vertix; j++)
		{
			printf("%d ", matrix[i*vertix+j]);
		}
		printf("\n");
	}

	return 0;
}