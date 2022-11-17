#include <mpi.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include "bmp.h"

using namespace std;

//Defines the number of smoothing operations
#define NSmooth 1000

/****************************************************************/
/*variable declaration：                                        */
/*  bmpHeader    ： BMP file header                             */
/*  bmpInfo      ： BMP file information                        */
/*  **BMPSaveData： Store pixel data to be written              */
/*  **BMPData    ： Temporarily store pixel data to be written  */
/****************************************************************/
BMPHEADER bmpHeader;
BMPINFO bmpInfo;
RGBTRIPLE **BMPSaveData = NULL;
RGBTRIPLE **BMPData = NULL;

/*********************************************************************************/
/*function declaration：                                                         */
/*  readBMP       ： Read the image file and store the pixel data in BMPSaveData */
/*  saveBMP       ： Write the image file and write the pixel data BMPSaveData   */
/*  swap          ： Swap two indicators                                         */
/*  **alloc_memory： Dynamically allocate a Y*X matrix                           */
/*********************************************************************************/
int readBMP( char* fileName);		//read file
int saveBMP( char* fileName);		//save file
void swap(RGBTRIPLE* a, RGBTRIPLE* b);
RGBTRIPLE** alloc_memory( int Y, int X );		//allocate memory

int main(int argc, char* argv[])
{
/*********************************************************/
/*variable declaration：                                 */
/*  *infileName  ： read filename                        */
/*  *outfileName ： write filename                       */
/*  startwtime   ： record start time                    */
/*  endwtime     ： record end time                      */
/*********************************************************/
		char *infileName = "input.bmp";
		char *outfileName = "output.bmp";
		double startwtime = 0.0,endwtime=0.0;

    	int id, numprocs;    			// process id and process size 

		int *sendNum;					// number of pcs to send
		int *sendDisplc;                // displacement (starting addr) for each proc 

		// init
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &id);

    	//read file
        readBMP( infileName);

		MPI_Datatype MPI_RGBTRIPLE;		// new MPI datatype for RGB
		int elements[3];				// num of elements per set	
		MPI_Datatype type[3] ;			// element's type
		MPI_Aint displc[3];				// displacement of each elm
		
		for (int i = 0; i < 3; i++) {
			elements[i] = 1;
			type[i] = MPI_UNSIGNED_CHAR;
			displc[i] = i; 
		}     
		MPI_Type_create_struct(3, elements, displc, type, &MPI_RGBTRIPLE);
		MPI_Type_commit(&MPI_RGBTRIPLE);

		//record start time
		startwtime = MPI_Wtime();
		MPI_Barrier(MPI_COMM_WORLD);
    
		// allocate mem
		RGBTRIPLE** upperBuff = alloc_memory(1, bmpInfo.biWidth);
		RGBTRIPLE** lowerBuff = alloc_memory(1, bmpInfo.biWidth);
		sendNum = (int *) malloc(sizeof(int) * numprocs);
		sendDisplc = (int *) malloc(sizeof(int) * numprocs);                    
     
		// distrib remaining pieces
		int rem = bmpInfo.biHeight % numprocs;
		int sum = 0;                                    
		for (int i = 0; i < numprocs; i++) {
			sendNum[i] = bmpInfo.biHeight * bmpInfo.biWidth / numprocs;
			if (rem > 0) {
				sendNum[i] += bmpInfo.biWidth;                             
				rem--; 
			}
			sendDisplc[i] = sum;                                          
			sum += sendNum[i];                                         
		}

		//Dynamically allocate memory to scratch space
		if (id == 0)
        	BMPData = alloc_memory( bmpInfo.biHeight, bmpInfo.biWidth);
		else 
        	BMPData = alloc_memory( sendNum[id] / bmpInfo.biWidth, bmpInfo.biWidth);

		// scatter data to other processes
		MPI_Scatterv(BMPSaveData[0], sendNum, sendDisplc, MPI_RGBTRIPLE, BMPData[0], sendNum[0], MPI_RGBTRIPLE, 0, MPI_COMM_WORLD);

		//Perform multiple smoothing operations
		for (int count = 0 ; count < NSmooth ; count ++) {
			if (count > 0)
				//Swap pixel data with staging metrics
				swap(BMPSaveData,BMPData);

			// send and obtain upper and lower bound
			int index = sendNum[id] / bmpInfo.biWidth - 1;
			int partner;
			// odd 
			if (id % 2 != 0) {
				// check if the next process is proc 0
				partner = id < numprocs - 1 ? id + 1 : 0;
				
				MPI_Send(BMPData[0], bmpInfo.biWidth, MPI_RGBTRIPLE, id - 1, 0, MPI_COMM_WORLD);
				MPI_Send(BMPData[index], bmpInfo.biWidth, MPI_RGBTRIPLE, partner, 0, MPI_COMM_WORLD);

				MPI_Recv(upperBuff[0], bmpInfo.biWidth, MPI_RGBTRIPLE, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(lowerBuff[0], bmpInfo.biWidth, MPI_RGBTRIPLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			} 
			// even
			else {
				// check if currently proc 0
				partner = id ? id - 1 : numprocs - 1;
				
				MPI_Send(BMPData[0], bmpInfo.biWidth, MPI_RGBTRIPLE, partner, 0, MPI_COMM_WORLD);
				MPI_Send(BMPData[index], bmpInfo.biWidth, MPI_RGBTRIPLE, id + 1, 0, MPI_COMM_WORLD);
				
				MPI_Recv(upperBuff[0], bmpInfo.biWidth, MPI_RGBTRIPLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(lowerBuff[0], bmpInfo.biWidth, MPI_RGBTRIPLE, id + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			//perform smoothing
			for (int i = 0; i < (sendNum[id] / bmpInfo.biWidth); i++) {
				for (int j = 0 ; j < bmpInfo.biWidth; j++) {
					/**************************************************************/
					/*Set the position of the upper, lower, left and right pixels */
					/**************************************************************/
					int Top = i>0 ? i-1 : bmpInfo.biHeight-1;
					int Down = i<bmpInfo.biHeight-1 ? i+1 : 0;
					int Left = j>0 ? j-1 : bmpInfo.biWidth-1;
					int Right = j<bmpInfo.biWidth-1 ? j+1 : 0;

					RGBTRIPLE t, d, l, r, tl, tr, dl, dr;
					// get data from all direction
					if (i == (sendNum[id] / bmpInfo.biWidth) - 1) {
						t = BMPData[Top][j];
						d = lowerBuff[0][j];
						l = BMPData[i][Left];
						r = BMPData[i][Right];
						tl = BMPData[Top][Left];
						tr = BMPData[Top][Right];
						dl = lowerBuff[0][Left];
						dr = lowerBuff[0][Right];
					}
					else if (i == 0) {
						t = upperBuff[0][j];
						d = BMPData[Down][j];
						l = BMPData[i][Left];
						r = BMPData[i][Right];
						tl = upperBuff[0][Left];
						tr = upperBuff[0][Right];
						dl = BMPData[Down][Left];
						dr = BMPData[Down][Right];
					}  
					else {
						t = BMPData[Top][j];
						d = BMPData[Down][j];
						l = BMPData[i][Left];
						r = BMPData[i][Right];
						tl = BMPData[Top][Left];
						tr = BMPData[Top][Right];
						dl = BMPData[Down][Left];
						dr = BMPData[Down][Right];
					}
					/****************************************************************/
					/*Average with top, bottom, left, and right pixels and round up */
					/****************************************************************/
					BMPSaveData[i][j].rgbBlue = (double)(BMPData[i][j].rgbBlue + t.rgbBlue + d.rgbBlue + l.rgbBlue + r.rgbBlue + tl.rgbBlue + tr.rgbBlue + dl.rgbBlue + dr.rgbBlue) / 9 + 0.5;
					BMPSaveData[i][j].rgbGreen = (double)(BMPData[i][j].rgbGreen + t.rgbGreen + d.rgbGreen + l.rgbGreen + r.rgbGreen + tl.rgbGreen + tr.rgbGreen + dl.rgbGreen + dr.rgbGreen) / 9 + 0.5;
					BMPSaveData[i][j].rgbRed = (double)(BMPData[i][j].rgbRed + t.rgbRed + d.rgbRed + l.rgbRed + r.rgbRed + tl.rgbRed + tr.rgbRed + dl.rgbRed + dr.rgbRed) / 9 + 0.5;
				}
			}
		}

		// gather each results back
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gatherv(BMPSaveData[0], sendNum[0], MPI_RGBTRIPLE, BMPData[0], sendNum, sendDisplc, MPI_RGBTRIPLE, 0, MPI_COMM_WORLD);
		
		//write file
		if (id == 0) { 
			if ( saveBMP( outfileName ) )
					cout << "Save file successfully!!" << endl;
			else
					cout << "Save file fails!!" << endl;
			
			//Get the end time and print the execution time
			endwtime = MPI_Wtime();
			cout << "The execution time = "<< endwtime-startwtime <<endl ;	
		}

		free(BMPData[0]);
		free(BMPSaveData[0]);
		free(BMPData);
		free(BMPSaveData);

		free(upperBuff[0]);
		free(lowerBuff[0]);
		free(upperBuff);
		free(lowerBuff);

		free(sendNum);
		free(sendDisplc);

		MPI_Type_free(&MPI_RGBTRIPLE);

		MPI_Finalize();

		return 0;
}

/*********************************************************/
/* read image                                            */
/*********************************************************/
int readBMP(char *fileName)
{
		//Create input file object
        ifstream bmpFile( fileName, ios::in | ios::binary );

        //File cannot be opened
        if ( !bmpFile ){
                cout << "It can't open file!!" << endl;
                return 0;
        }

        //Read the header data of the BMP file
    	bmpFile.read( ( char* ) &bmpHeader, sizeof( BMPHEADER ) );

        //Determine whether it is a BMP image file
        if( bmpHeader.bfType != 0x4d42 ){
                cout << "This file is not .BMP!!" << endl ;
                return 0;
        }

        //Read BMP information
        bmpFile.read( ( char* ) &bmpInfo, sizeof( BMPINFO ) );

        //Determine if the bit depth is 24 bits
        if ( bmpInfo.biBitCount != 24 ){
                cout << "The file is not 24 bits!!" << endl;
                return 0;
        }

        //Correct the width of the image to be a multiple of 4
        while( bmpInfo.biWidth % 4 != 0 )
        	bmpInfo.biWidth++;

        //Dynamically allocate memory
        BMPSaveData = alloc_memory( bmpInfo.biHeight, bmpInfo.biWidth);

        //read pixel data
    	//for(int i = 0; i < bmpInfo.biHeight; i++)
        //	bmpFile.read( (char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE));
	    bmpFile.read( (char* )BMPSaveData[0], bmpInfo.biWidth*sizeof(RGBTRIPLE)*bmpInfo.biHeight);

        //close file
        bmpFile.close();

        return 1;

}
/*********************************************************/
/* save image                                            */
/*********************************************************/
int saveBMP( char *fileName)
{
 		//Determine whether it is a BMP image file
        if( bmpHeader.bfType != 0x4d42 ){
                cout << "This file is not .BMP!!" << endl ;
                return 0;
        }

 		//Create output file object
        ofstream newFile( fileName,  ios:: out | ios::binary );

        //File could not be created
        if ( !newFile ){
                cout << "The File can't create!!" << endl;
                return 0;
        }

        //Write the header data of the BMP file
        newFile.write( ( char* )&bmpHeader, sizeof( BMPHEADER ) );

		//Information written to BMP
        newFile.write( ( char* )&bmpInfo, sizeof( BMPINFO ) );

        //write pixel data
        //for( int i = 0; i < bmpInfo.biHeight; i++ )
        //        newFile.write( ( char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE) );
		newFile.write((char*)BMPData[0], bmpInfo.biWidth * sizeof(RGBTRIPLE) * bmpInfo.biHeight);
		
		//write file
		newFile.close();

		return 1;
}

/*********************************************************/
/* Allocate memory: return a matrix of Y*X               */
/*********************************************************/
RGBTRIPLE **alloc_memory(int Y, int X )
{
		//Create an index array of length Y
        RGBTRIPLE **temp = new RGBTRIPLE *[ Y ];
		RGBTRIPLE *temp2 = new RGBTRIPLE [ Y * X ];
        memset( temp, 0, sizeof( RGBTRIPLE ) * Y);
        memset( temp2, 0, sizeof( RGBTRIPLE ) * Y * X );

		//declares an array of length X for each indicator in the indicator array
        for( int i = 0; i < Y; i++) {
                temp[ i ] = &temp2[i*X];
        }

        return temp;

}
/*********************************************************/
/* Swap two indicators                                   */
/*********************************************************/
void swap(RGBTRIPLE *a, RGBTRIPLE *b)
{
		RGBTRIPLE *temp;
		temp = a;
		a = b;
		b = temp;
}
