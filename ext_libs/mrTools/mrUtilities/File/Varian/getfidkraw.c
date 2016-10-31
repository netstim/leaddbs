#ifdef documentation
=========================================================================

    program: getfidkraw.c
         by: justin gardner
    purpose: reads fid data from VNMR system - does not try to do
             any conversion of the data - this is left for matlab
             programs to do. Based on getfidk.c
       date: 08/10/2011
    compile: mex getfidkraw.c
=========================================================================
#endif

////////////////////
// include section
////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "mex.h"
#include "vnmrdata.h"

///////////////////
// define section
///////////////////
#define STRSIZE 2048
#define FALSE 0
#define TRUE 1
#define INT16 int16
#define INT32 int32
#define FLOAT float

///////////////////
// function decls
///////////////////
void usageError();
int fileLength(char *,FILE *);
void errorExit(mxArray *[]);
void swapBytes(void * p, size_t s);

////////////////////
// global variables
////////////////////
FILE *ffid = NULL; // declared global so we can close files on error

//////////////////////////////////////////
// function mexFunction called by matlab
//////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char filename[STRSIZE];
  struct datafilehead header;
  struct datablockhead blockheader;
  INT16 *data_int16,*datai_int16,*raw_int16;
  INT32 *data_int32,*datai_int32,*raw_int32;
  FLOAT *data_float,*datai_float,*raw_float;
  double *data,*datai;
  INT16 *block;
  mxArray *mxdata,*mxdatar,*mxdatai;
  mxClassID dataType;
  int i,j,k,n;
  int isLittleEndianPlatform, swapFlag;
  int verbose,headerOnly;
  
  // check input arguments
  if ((nrhs < 1) || (nrhs > 3)){
    usageError();
    return;
  }
  // check output arguments
  else if (nlhs > 1) {
    usageError();
    return;
  }

  // set whether to verbosely report errors
  if (nrhs >= 2) {
    verbose = (int)*mxGetPr(prhs[1]);
  }
  else {
    verbose = FALSE;
  }

  // set whether to verbosely report errors
  if (nrhs >= 3) {
    headerOnly = (int)*mxGetPr(prhs[2]);
  }
  else {
    headerOnly = FALSE;
  }

  // get filename
  mxGetString(prhs[0], filename, mxGetN(prhs[0])+1);

  // create the output structure 
  const char *fieldNames[] = {"name","nblocks","ntraces","np","ebytes","tbytes","bbytes","vers_id","status","nbheaders","real","imag"};
  int dims[2] = {1, 1};
  plhs[0] = mxCreateStructArray(2, dims, 12, fieldNames);

  // put filename in output structure
  mxSetField(plhs[0],0,"name",mxCreateString(filename)); 
 
  // open fid file
  if ((ffid = fopen(filename,"r")) == NULL) {
    mexPrintf("(getfidkraw) Could not open file %s\n",filename);
    return;
  }
  if (verbose) mexPrintf("(getfidkraw) Opened %s\n", filename);

  // get the fid file length
  if ((n = fileLength(filename,ffid)) == -1) {
    errorExit(plhs);
    return;
  }

  // read the fid file header
  if (fread(&header, sizeof(struct datafilehead), 1, ffid) == 0) {
    mexPrintf("(getfidkraw) ERROR: reading header from file %s\n",filename);
    errorExit(plhs);
    return;
  }
  
  // try to see which endian platform we are running on 
  isLittleEndianPlatform = 0;
  int one=1;
  if( *(char*)&one == 1 ) 
    isLittleEndianPlatform = 1;
  
  swapFlag = 0;
  // then make sure fid is in big endian
  if ( isLittleEndianPlatform & header.nbheaders > 9 ) {
    swapFlag = 1;
    if (verbose)
      mexPrintf("(getfidkraw) Running on little endian platform and data is big endian (default), will swap bytes\n");
  }
  else
    swapFlag = 0;
  
  if ( swapFlag ) {
    swapBytes(&header.nblocks, 4); // long
    swapBytes(&header.ntraces, 4);
    swapBytes(&header.np, 4);
    swapBytes(&header.ebytes, 4);
    swapBytes(&header.tbytes, 4);
    swapBytes(&header.bbytes, 4);
    swapBytes(&header.vers_id, 2);  // short
    swapBytes(&header.status, 2);  // short
    swapBytes(&header.nbheaders, 4);  // long
  }

  // set header info in output structure
  mxSetField(plhs[0],0,"nblocks",mxCreateDoubleScalar(header.nblocks)); 
  mxSetField(plhs[0],0,"ntraces",mxCreateDoubleScalar(header.ntraces));
  mxSetField(plhs[0],0,"np",mxCreateDoubleScalar(header.np));
  mxSetField(plhs[0],0,"ebytes",mxCreateDoubleScalar(header.ebytes));
  mxSetField(plhs[0],0,"tbytes",mxCreateDoubleScalar(header.tbytes));
  mxSetField(plhs[0],0,"bbytes",mxCreateDoubleScalar(header.bbytes));
  mxSetField(plhs[0],0,"vers_id",mxCreateDoubleScalar(header.vers_id));
  mxSetField(plhs[0],0,"nbheaders",mxCreateDoubleScalar(header.nbheaders));

  const char *statusFieldNames[] = {"data","spec","is32","float","complex","hypercomplex","acqpar","secnd","transf","is3d","np","nf","ni","ni2"};
  mxSetField(plhs[0],0,"status",mxCreateStructArray(2, dims, 14, statusFieldNames));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"data",mxCreateDoubleScalar((header.status&S_DATA)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"spec",mxCreateDoubleScalar((header.status&S_SPEC)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"is32",mxCreateDoubleScalar((header.status&S_32)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"float",mxCreateDoubleScalar((header.status&S_FLOAT)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"complex",mxCreateDoubleScalar((header.status&S_COMPLEX)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"hypercomplex",mxCreateDoubleScalar((header.status&S_HYPERCOMPLEX)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"acqpar",mxCreateDoubleScalar((header.status&S_ACQPAR)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"secnd",mxCreateDoubleScalar((header.status&S_SECND)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"transf",mxCreateDoubleScalar((header.status&S_TRANSF)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"is3d",mxCreateDoubleScalar((header.status&S_3D)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"np",mxCreateDoubleScalar((header.status&S_NP)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"nf",mxCreateDoubleScalar((header.status&S_NF)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"ni",mxCreateDoubleScalar((header.status&S_NI)!=0));
  mxSetField(mxGetField(plhs[0],0,"status"),0,"ni2",mxCreateDoubleScalar((header.status&S_NI2)!=0));


  // display fid file header
  if (verbose) {
    mexPrintf("(getfidkraw) FID HEADER\n");
    mexPrintf("  nblocks   = %li\n",header.nblocks);
    mexPrintf("  ntraces   = %li\n",header.ntraces);
    mexPrintf("  np        = %li\n",header.np);
    mexPrintf("  ebytes    = %li\n",header.ebytes);
    mexPrintf("  tbytes    = %li\n",header.tbytes);
    mexPrintf("  bbytes    = %li\n",header.bbytes);
    mexPrintf("  vers_id   = %i\n",header.vers_id);
    mexPrintf("  status    = %i\n",header.status);
    mexPrintf("  nbheaders = %li\n",header.nbheaders);
    mexPrintf("(getfidkraw)  STATUS DECODE:\n");
    mexPrintf("    S_DATA         = %i\n",(header.status&S_DATA)!=0);
    mexPrintf("    S_SPEC         = %i\n",(header.status&S_SPEC)!=0);
    mexPrintf("    S_32           = %i\n",(header.status&S_32)!=0);
    mexPrintf("    S_FLOAT        = %i\n",(header.status&S_FLOAT)!=0);
    mexPrintf("    S_COMPLEX      = %i\n",(header.status&S_COMPLEX)!=0);
    mexPrintf("    S_HYPERCOMPLEX = %i\n",(header.status&S_HYPERCOMPLEX)!=0);
    mexPrintf("    S_ACQPAR       = %i\n",(header.status&S_ACQPAR)!=0);
    mexPrintf("    S_SECND        = %i\n",(header.status&S_SECND)!=0);
    mexPrintf("    S_TRANSF       = %i\n",(header.status&S_TRANSF)!=0);
    mexPrintf("    S_3D           = %i\n",(header.status&S_3D)!=0);
    mexPrintf("    S_NP           = %i\n",(header.status&S_NP)!=0);
    mexPrintf("    S_NF           = %i\n",(header.status&S_NF)!=0);
    mexPrintf("    S_NI           = %i\n",(header.status&S_NI)!=0);
    mexPrintf("    S_NI2          = %i\n",(header.status&S_NI2)!=0);
  }

  // just return header if called for
  if (headerOnly) return;

  // figure out data class
  if (header.status & S_FLOAT) {
    if (verbose) mexPrintf("(getfidkraw) datatype is FLOAT\n");
    // set data type
    dataType = mxSINGLE_CLASS;
  }
  else {
    if (header.status & S_32) {
      dataType = mxINT32_CLASS;
      if (verbose) mexPrintf("(getfidkraw) datatype is int32\n");
    }
    else {
      dataType = mxINT16_CLASS;
      if (verbose) mexPrintf("(getfidkraw) datatype is int16\n");
    }
  }

  // make sure we have enough filelength
  if ((sizeof(struct datafilehead)+header.bbytes*header.nblocks) != n) {
    mexPrintf("(getfidkraw) ERROR: Expected %i bytes, but found %i\n",sizeof(struct datafilehead)+header.bbytes*header.nblocks,n);
    errorExit(plhs);
    return;
  }

  // set space for real data in its original format
  if (verbose)
    mexPrintf("(getfidkraw) Allocating space for %ix%i data\n",(int)(header.nblocks),(int)(header.np/2));

  int dataDims[2] = {header.ntraces*header.np/2,header.nblocks};
  if ((mxdatar = mxCreateNumericArray(2,dataDims,dataType,mxREAL)) == NULL) {
    errorExit(plhs);
    return;
  }
  mxSetField(plhs[0],0,"real",mxdatar);
  // set space for imaginary data in its original format
  if ((mxdatai = mxCreateNumericArray(2,dataDims,dataType,mxREAL)) == NULL) {
    errorExit(plhs);
    return;
  }
  mxSetField(plhs[0],0,"imag",mxdatai);

  // get space to hold a block of data
  block = (INT16 *)malloc(header.tbytes*header.ntraces);

  // get pointers to places to store data
  switch (dataType) {
    case mxSINGLE_CLASS:
      data_float = (FLOAT*)mxGetPr(mxGetField(plhs[0],0,"real"));
      datai_float = (FLOAT*)mxGetPr(mxGetField(plhs[0],0,"imag"));
      break;
    case mxINT16_CLASS:
      data_int16 = (INT16 *)mxGetPr(mxGetField(plhs[0],0,"real"));
      datai_int16 = (INT16 *)mxGetPr(mxGetField(plhs[0],0,"imag"));
      break;
    case mxINT32_CLASS:
      data_int32 = (INT32 *)mxGetPr(mxGetField(plhs[0],0,"real"));
      datai_int32 = (INT32 *)mxGetPr(mxGetField(plhs[0],0,"imag"));
      break;
  }

  // read each line of k-space until the end of file
  for (i=0;i<header.nblocks;i++) {
    // read block header if nbheaders is 1 or 9 (9 is a code set by epibsi meaning
    // the same as 1, except specifiying the epibsi processing has been run)
    if ((header.nbheaders == 1) || (header.nbheaders == 9)) {
      if (fread(&blockheader, sizeof(struct datablockhead), 1, ffid) == 0) {
	mexPrintf("(getfidkraw) ERROR: Could not read block header from file %s\n",filename);
	errorExit(plhs);
	free(block);
	return;
      }
    
      //swap bytes  
      if ( swapFlag ) {
	swapBytes(&blockheader.scale, 2); // short
	swapBytes(&blockheader.status, 2);
	swapBytes(&blockheader.index, 2);
	swapBytes(&blockheader.mode, 2);
	swapBytes(&blockheader.ctcount, 4); // long
	swapBytes(&blockheader.lpval, 4);   // float
	swapBytes(&blockheader.rpval, 4);  // float
	swapBytes(&blockheader.lvl, 4);  // float
	swapBytes(&blockheader.tlt, 4);  // float
      }
      
      // print out block headers
      if (verbose > 1) {
	mexPrintf("(getfidkraw)  BLOCK %i\n",i);
	mexPrintf("    scale   = %i\n",blockheader.scale);
	mexPrintf("    status  = %i\n",blockheader.status);
	mexPrintf("    index   = %i\n",blockheader.index);
	mexPrintf("    mode    = %i\n",blockheader.mode);
	mexPrintf("    ctcount = %li\n",blockheader.ctcount);
	mexPrintf("    lpval   = %f\n",blockheader.lpval);
	mexPrintf("    rpval   = %f\n",blockheader.rpval);
	mexPrintf("    lvl     = %f\n",blockheader.lvl);
	mexPrintf("    tlt     = %f\n",blockheader.tlt);
      }
    }
    else if (header.nbheaders != 0){
      mexPrintf("(getfidkraw) Don't know how to handle number of block headers setting to: %i (should be 0 or 1)\n",header.nbheaders);
      errorExit(plhs);
      free(block);
      return;
    }

    // read in the block of data
    if (fread(block, 1, header.tbytes*header.ntraces, ffid) == 0) {
      mexPrintf("(getfidkraw) ERROR: reading data from file %s\n",filename);
      errorExit(plhs);
      free(block);
      return;
    }
    // move block into data and datai pointer approriately
    // note that the real and imaginary parts of the data
    // are stored sequentially in the fid file. 
    // They are stored in its native format into two arrays
    for (k=0; k<header.np*header.ntraces; k+=2) {
      switch(dataType) {
	// data is stored as a float
        case mxSINGLE_CLASS:
	  // convert endian
	  if (swapFlag) {
	    swapBytes(((FLOAT*)block)+k,4); 
	    swapBytes(((FLOAT*)block)+k+1,4); 
	  }
	  // save data in ouput array
	  *data_float++ = ((FLOAT*)block)[k];
	  *datai_float++ = ((FLOAT*)block)[k+1];
	  break;
	  // data is stored as two byte integers
	case mxINT16_CLASS:
	  // convert endian
	  if (swapFlag) {
	    swapBytes(((INT16*)block)+k,2); 
	    swapBytes(((INT16*)block)+k+1,2); 
	  }
	  // save data in ouput array
	  *data_int16++ = ((INT16*)block)[k];
	  *datai_int16++ = ((INT16*)block)[k+1];
	  break;
	  // data is stored as four byte integers
	case mxINT32_CLASS:
	  // convert endian
	  if (swapFlag) {
	    swapBytes(((INT32*)block)+k,4); 
	    swapBytes(((INT32*)block)+k+1,4); 
          }  
	  // save data in ouput array
	  *data_int32++ = ((INT32*)block)[k];
	  *datai_int32++ = ((INT32*)block)[k+1];
	  break;
      }
    }
  }
  // free memory for block
  free(block);
  // close file handles
  fclose(ffid);
}

////////////////////////
// function usageError
////////////////////////
void usageError()
{
  mxArray *callInput[] = {mxCreateString("getfidkraw")};
  mexCallMATLAB(0,NULL,1,callInput,"help");
}

///////////////////////
// function errorexit
///////////////////////
void errorExit(mxArray *plhs[])
{
  if (ffid) fclose(ffid);
  mxSetField(plhs[0],0,"data",NULL);
}

////////////////////////
// function fileLength
////////////////////////
int fileLength(char *filename,FILE *filepointer)
{
  fpos_t filelen;
  // get size of file: seek to end of file, read pos, seek to begin of file
  if (fseek(filepointer,0,SEEK_END)) {
    mexPrintf("(getfidkraw) ERROR: Could not seek in file %s",filename);
    return -1;
  }

  // find the position, i.e. end of file
  // this may not work with all compilers since fpos_t
  // is not required to be the position in bytes.
  if (fgetpos(filepointer,&filelen)) {
    mexPrintf("(getfidkraw) ERROR: Could not seek in file %s",filename);
    return -1;
  }
  
  // seek back to the beginning of the file
  if (fseek(filepointer,0,SEEK_SET)) {
    mexPrintf("(getfidkraw) ERROR: Could not seek in file %s",filename);
    return -1;
  }
  return (int)filelen;
}  


//////////////////////////////////////////
// function swap bytes
//////////////////////////////////////////
void swapBytes(void * p, size_t s)
{
  unsigned char tmp, *a, *b ;
 
  a = (unsigned char*)p ;
  b = a + s ;

  while (a<b) {
    tmp = *a ;
    *a++ = *--b ;
    *b = tmp ;
  }
}


