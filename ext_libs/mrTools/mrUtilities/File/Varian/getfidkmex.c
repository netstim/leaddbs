#ifdef documentation
=========================================================================

    program: getfidkmex.c
         by: justin gardner
    purpose: reads fid data from VNMR system
       date: 05/06/03
    compile: mex getfidk.c
    history: 2008/12/01, Pei, add endian detect and byte swap part.
             Default format of fid is big endian, as we acquaire data
             on VNMR Sun Solaris, also epirri5 keep the same endian.
             Note: SWAP data firstly and then DO typecast.

             This is an obsolete function that has been replaced by
             getfidkraw.c and getfidk.m
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
#define LINESIZE 4096
#define FALSE 0
#define TRUE 1
#define INT16 int16
#define INT32 int32
#define FLOAT float

////////////////////
// data structures
////////////////////
typedef struct {
  int numslices;                  // number of slices
  int numvolumes;                 // number of volumes
  int numlines;                   // number of k-space lines in each image
  int numreceivers;               // number of receivers
  int numshots;                   // number of shots
  int linelen;                    // number of voxels in each line
  int navechoes;                  // number of naviagtor echoes
  int accFactor;                   // acceleration factor
  int ndim;                       // number of dimensions
  int dims[5];                    // size of those dimensions
  int numimages;                  // number of images = slices*volumes*receivers
  int imagesize;                  // number of voxels in each image
  mxClassID datatype;             // matlab type of data
  int totalsize;                  // total size of array
} datainfo;

///////////////////
// function decls
///////////////////
void usageError();
void fileNotFound(char *);
int fileLength(char *,FILE *);
void errorExit(mxArray *[]);
// functions of bytes swap
void swap_bytes(void * p, size_t s);

////////////////////
// global variables
////////////////////
int verbose;
FILE *ffid = NULL,*fprocpar = NULL; // declared global so we can close files on error

/////////////////////////
// for output structure
/////////////////////////
const char *fieldNames[] = {"name","filepath","data"};
#define NUMFIELDS 3
const char *originalFormatFieldNames[] = {"name","filepath","real","imag"};
#define ORIGINAL_FORMAT_NUMFIELDS 4
int dims[2] = {1, 1};

//////////////////////////////////////////
// function mexFunction called by matlab
//////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char filename[STRSIZE],inputname[STRSIZE],
    pathname[STRSIZE],procparname[STRSIZE],
    fullpath[STRSIZE],line[LINESIZE],*token;
  struct datafilehead header;
  struct datablockhead blockheader;
  datainfo info;
  INT16 *data_int16,*datai_int16,*raw_int16, tokenint16;
  INT32 *data_int32,*datai_int32,*raw_int32, tokenint32;
  FLOAT *data_float,*datai_float,*raw_float, tokenfloat;
  double *data,*datai;
  INT16 *block;
  mxArray *mxdata,*mxdatar,*mxdatai;
  int originalformat;
  int i,j,k,n;
  int isLittleEndianPlatform, swapFlag, gotNavechoes;
  
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

  // set whether to return data in original format or not
  if (nrhs >= 3) {
    originalformat = (int)*mxGetPr(prhs[2]);
  }
  else {
    originalformat = FALSE;
  }

  // get filename
  mxGetString(prhs[0], inputname, mxGetN(prhs[0])+1);

  // create the output structure 
  if (originalformat) 
    plhs[0] = mxCreateStructArray(2, dims, ORIGINAL_FORMAT_NUMFIELDS, originalFormatFieldNames);
  else
    plhs[0] = mxCreateStructArray(2, dims, NUMFIELDS, fieldNames);

  // put filename in output structure
  mxSetField(plhs[0],0,"name",mxCreateString(inputname)); 
 
  // figure out where the fid file is
  // did the user input the full directory?
  sprintf(pathname,"%s",inputname);
  sprintf(filename,"%s/fid",pathname);
  if (verbose) mexPrintf("(getfidkmex) Trying filename %s\n",filename);
  if ((ffid = fopen(filename,"r")) == NULL) {
    // maybe the user gave the directory without .fid
    sprintf(pathname,"%s.fid",inputname);
    sprintf(filename,"%s/fid",pathname);
    if (verbose) mexPrintf("(getfidkmex) Trying filename %s\n",filename);
    if ((ffid = fopen(filename,"r")) == NULL) {
      // nope, maybe it is simply the filename
      sprintf(pathname,".");
      sprintf(filename,"%s/%s",pathname,inputname);
      if (verbose) mexPrintf("(getfidkmex) Trying filename %s\n",filename);
      if ((ffid = fopen(filename,"r")) == NULL) {
	// nothing found, give up.
	fileNotFound(inputname);
	errorExit(plhs);
	return;
      }
    }
  }
  if (verbose) mexPrintf("(getfidkmex) Opened %s\n", filename);

  // get the full path
  getcwd(fullpath,STRSIZE);
  strcat(fullpath,"/");
  if (!strcmp(pathname,".")) {
    if (!strcmp(inputname,"."))
      strcat(fullpath,"fid");
    else
      strcat(fullpath,inputname);
  }
  else
    strcat(fullpath,filename);

  // save the full qualified path in the output structure
  mxSetField(plhs[0],0,"filepath",mxCreateString(fullpath)); 

  // open the procpar
  sprintf(procparname,"%s/procpar",pathname);
  if ((fprocpar = fopen(procparname,"r")) == NULL) {
    // no procpar
    fileNotFound(filename);
    errorExit(plhs);
    return;
  }

  info.navechoes = 0;
  gotNavechoes = 0;
  info.numshots = 1;
  info.accFactor = 1;
  info.numvolumes = 1;
  info.numreceivers = 1;
  // scan through procpar looking for parameters needed for parsing
  while (fgets(line,LINESIZE,fprocpar)) {
    // get the first token of the line
    token = (char*)strtok(line," ");
    // get number of slices
    if (!strcmp(token,"pss")) {
      info.numslices = atoi(fgets(line,LINESIZE,fprocpar));
      if (verbose) mexPrintf("(getfidkmex) numslices = %i\n",info.numslices);
    }
    // get number of volumes
    if (!strcmp(token,"cntr")) {
      info.numvolumes = atoi(fgets(line,LINESIZE,fprocpar));
      if (verbose) mexPrintf("(getfidkmex) numvolumes = %i\n",info.numvolumes);
    }
    // get number of voxels per line
    if (!strcmp(token,"np")) {
      fgets(line,LINESIZE,fprocpar);
      token = (char*)strtok(line," ");
      token = (char*)strtok(NULL," ");
      // that should be number of slices
      info.linelen = atoi(token)/2;
      if (verbose) mexPrintf("(getfidkmex) linelen = %i\n",info.linelen);
    }
    // get the number of receivers
    if (!strcmp(token,"rcvrs")) {
      fgets(line,LINESIZE,fprocpar);
      token = (char*)strtok(line," ");
      token = (char*)strtok(NULL," ");
      info.numreceivers = 0;
      // Count the number of 'y's that is the number of receivers
      for (i = 0;i < strlen(token);i++)
	if (token[i] == 'y') info.numreceivers++;
      // make sure we have atleast 1 receiver
      if (info.numreceivers <= 0) info.numreceivers = 1;
      if (verbose) mexPrintf("(getfidkmex) Number of receivers = %i\n",info.numreceivers);
    }
    // get number of lines per image
    if (!strcmp(token,"ni")) {
      fgets(line,LINESIZE,fprocpar);
      token = (char*)strtok(line," ");
      token = (char*)strtok(NULL," ");
      info.numlines = atoi(token);
      if (verbose) mexPrintf("(getfidkmex) numlines = %i\n",info.numlines);
    }
    // get number of navechoes
    if (!strcmp(token,"nv")) {
      fgets(line,LINESIZE,fprocpar);
      token = (char*)strtok(line," ");
      token = (char*)strtok(NULL," ");
      // nv should be the number of lines of k-space including the nav echoes.
      // later we subtract of the lines of k-space found from ni to get the number of nav echoes
      info.navechoes = atoi(token);
      if (verbose) mexPrintf("(getfidkmex) nv (from procpar) = %i\n",info.navechoes);
    }
    // for epi images, there are navigator echos, which
    // should be subtracted from the number of lines.
    // this can be known from the name of the petable
    if (!strcmp(token,"petable")) {
      fgets(line,LINESIZE,fprocpar);
      token = (char*)strtok(line," ");
      token = (char*)strtok(NULL," ");
      if (verbose) mexPrintf("(getfidkmex) petable name: %s",token);
      // the petable name should be something like
      // "epi132alt8k". We want the second number
      if (!strncmp(token+1,"epi",3)) {
	j = 0;
	// go past any intial characters
	while((j < strlen(token)) && !isdigit(token[j])) j++;
	// then skip the numbers
	while((j < strlen(token)) && isdigit(token[j])) j++;
	// and skip the next characters
	while((j < strlen(token)) && !isdigit(token[j])) j++;
	// this is the number of shots
	info.numshots = atoi(token+j);
	if (verbose) mexPrintf("(getfidkmex) Number of shots (from petable name) = %i\n",info.numshots);
	// get the accelaration factor, should be like epi132alt8k2r
	while((j < strlen(token)) && isdigit(token[j])) j++;
	while((j < strlen(token)) && !isdigit(token[j])) j++;
	if (j < strlen(token)) {
	  info.accFactor = atoi(token+j);
	  if (verbose) mexPrintf("(getfidkmex) Acceleration factor (from petable name) = %i\n",info.accFactor);
	}
      }
    }
  }

  // now compute the number of navigator echoes total
  info.navechoes = info.navechoes-info.numlines;
  if (verbose) mexPrintf("(getfik) Number of navigator echoes: %i\n",info.navechoes);
  
  
  // set up other parts of info
  // set the dimensions for our array
  info.ndim = 5;
  info.dims[0] = info.linelen;
  info.dims[1] = info.numlines;
  info.dims[2] = info.numslices;
  info.dims[3] = info.numreceivers;
  info.dims[4] = info.numvolumes;
  // calculate some handy parameters
  info.totalsize = info.linelen*info.numlines*info.numslices*info.numvolumes*info.numreceivers;
  info.numimages = info.numslices*info.numvolumes*info.numreceivers;
  info.imagesize = info.linelen*info.numlines;
  // display parameters if verbose mode
  if (verbose) mexPrintf("(getfidkmex) (%ix%ix%ix%ix%i) numimages=%i,imagesize=%i\n",info.linelen,info.numlines,info.numslices,info.numvolumes,info.numreceivers,info.numimages,info.imagesize);    
  if (verbose) mexPrintf("(getfidkmex) totalsize = %i\n",info.totalsize);

  // get the fid file length and number of blocks
  if ((n = fileLength(filename,ffid)) == -1) {
    errorExit(plhs);
    return;
  }

  // read the fid file header
  if (fread(&header, sizeof(struct datafilehead), 1, ffid) == 0) {
    mexPrintf("(getfidkmex) ERROR: reading header from file %s\n",filename);
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
      mexPrintf("(getfidkmex) Running on little endian platform and data is big endian (default), will swap bytes\n");
  }
  else
    swapFlag = 0;
  
  if ( swapFlag ) {
      
    swap_bytes(&header.nblocks, 4); // long
    swap_bytes(&header.ntraces, 4);
    swap_bytes(&header.np, 4);
    swap_bytes(&header.ebytes, 4);
    swap_bytes(&header.tbytes, 4);
    swap_bytes(&header.bbytes, 4);
    swap_bytes(&header.vers_id, 2);  // short
    swap_bytes(&header.status, 2);  // short
    swap_bytes(&header.nbheaders, 4);  // long
  }
  
  // display fid file header
  if (verbose) {
    mexPrintf("(getfidkmex) FID HEADER\n");
    mexPrintf("  nblocks   = %li\n",header.nblocks);
    mexPrintf("  ntraces   = %li\n",header.ntraces);
    mexPrintf("  np        = %li\n",header.np);
    mexPrintf("  ebytes    = %li\n",header.ebytes);
    mexPrintf("  tbytes    = %li\n",header.tbytes);
    mexPrintf("  bbytes    = %li\n",header.bbytes);
    mexPrintf("  vers_id   = %i\n",header.vers_id);
    mexPrintf("  status    = %i\n",header.status);
    mexPrintf("  nbheaders = %li\n",header.nbheaders);
    if (verbose > 1) {
      mexPrintf("(getfidkmex)  STATUS DECODE:\n");
      mexPrintf("    S_DATA         = %i\n",header.status&S_DATA);
      mexPrintf("    S_SPEC         = %i\n",header.status&S_SPEC);
      mexPrintf("    S_32           = %i\n",header.status&S_32);
      mexPrintf("    S_FLOAT        = %i\n",header.status&S_FLOAT);
      mexPrintf("    S_COMPLEX      = %i\n",header.status&S_COMPLEX);
      mexPrintf("    S_HYPERCOMPLEX = %i\n",header.status&S_HYPERCOMPLEX);
      mexPrintf("    S_ACQPAR       = %i\n",header.status&S_ACQPAR);
      mexPrintf("    S_SECND        = %i\n",header.status&S_SECND);
      mexPrintf("    S_TRANSF       = %i\n",header.status&S_TRANSF);
      mexPrintf("    S_3D           = %i\n",header.status&S_3D);
      mexPrintf("    S_NP           = %i\n",header.status&S_NP);
      mexPrintf("    S_NF           = %i\n",header.status&S_NF);
      mexPrintf("    S_NI           = %i\n",header.status&S_NI);
      mexPrintf("    S_NI2          = %i\n",header.status&S_NI2);
    }
  }

  if ((header.ntraces > 1) && (header.nbheaders != 9)) {
    mexPrintf("(getfidkmex) Epirri5 has not been run on this fid file (ntraces = %i and nbheader=%i)\n",header.ntraces,header.nbheaders);
    errorExit(plhs);
    return;
  }

  // figure out data class
  if (header.status & S_FLOAT) {
    if (verbose) mexPrintf("(getfidkmex) datatype is FLOAT\n");
    // set data type
    info.datatype = mxSINGLE_CLASS;
  }
  else {
    if (header.status & S_32) {
      info.datatype = mxINT32_CLASS;
      if (verbose) mexPrintf("(getfidkmex) datatype is int32\n");
    }
    else {
      info.datatype = mxINT16_CLASS;
      if (verbose) mexPrintf("(getfidkmex) datatype is int16\n");
    }
  }

  // make sure we have enough filelength
  if ((sizeof(struct datafilehead)+header.bbytes*header.nblocks) != n) {
    mexPrintf("(getfidkmex) ERROR: Expected %i bytes, but found %i\n",sizeof(struct datafilehead)+header.bbytes*header.nblocks,n);
    errorExit(plhs);
    return;
  }

  if (originalformat) {
    // set space for real data in its original format
    if ((mxdatar = mxCreateNumericArray(info.ndim,info.dims,info.datatype,mxREAL)) == NULL) {
      errorExit(plhs);
      return;
    }
    mxSetField(plhs[0],0,"real",mxdatar);
    // set space for imaginary data in its original format
    if ((mxdatai = mxCreateNumericArray(info.ndim,info.dims,info.datatype,mxREAL)) == NULL) {
      errorExit(plhs);
      return;
    }
    mxSetField(plhs[0],0,"imag",mxdatai);
  }
  else {
    // set space for data as complex doubles
    if ((mxdata = mxCreateNumericArray(info.ndim,info.dims,mxDOUBLE_CLASS,mxCOMPLEX)) == NULL) {
      errorExit(plhs);
      return;
    }
    mxSetField(plhs[0],0,"data",mxdata);
  }
  // get space to hold a block of data
  block = (INT16 *)malloc(header.tbytes);

  // get pointers to places to store data
  if (originalformat) {
    switch (info.datatype) {
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
  }
  else {
    // set pointer to data
    data = (double*)mxGetPr(mxGetField(plhs[0],0,"data"));
    datai = (double*)mxGetPi(mxGetField(plhs[0],0,"data"));
  }

  // if the number of slices and volumes doesn't divide nicely
  // into the number of blocks, then something is wrong.
  if (((double)header.nblocks/info.numimages) != floor((double)header.nblocks/info.numimages)) {
    mexPrintf("(getfidkmex) ERROR: numimages (%i) does not divide evenly into num blocks %i\n",info.numimages,header.nblocks);
    errorExit(plhs);
    free(block);
    return;
  }

  // check that the number of blocks we have corresponds to how much data we think
  // we should have. This might be wrong if not all the nav echos have been correctly removed
  // from the data.
  if (((double)header.nblocks/info.numimages) != info.numlines) {
    mexPrintf("(getfidkmex) ERROR: File contains %0.0f k-space lines when %i are expected. Have you run the correct version of epirri to removes all %i navigator echoes?\n",(double)header.nblocks/info.numimages,info.numlines,info.navechoes);
    errorExit(plhs);
    free(block);
    return;
  }

  // the data are stored such that we read one k-space line
  // for each image, cyclying through all slices and volumes
  // then we read the next k-space line etc.
  for (i=0;i<info.numlines;i++) {
    if (verbose >= 2)
      mexPrintf("(getfidkmex) Reading line %i/%i\n",i,info.numlines);
    for (j=0;j<info.numimages;j++) {

      // read block header
      if ((header.ntraces == 1) || ((j % header.ntraces) == 0)) {
	if (fread(&blockheader, sizeof(struct datablockhead), 1, ffid) == 0) {
	  mexPrintf("(getfidkmex) ERROR: Could not read block header from file %s\n",filename);
	  errorExit(plhs);
	  free(block);
	  return;
	}
      }
    
      //swap bytes  
      if ( swapFlag ) {
	swap_bytes(&blockheader.scale, 2); // short
	swap_bytes(&blockheader.status, 2);
	swap_bytes(&blockheader.index, 2);
	swap_bytes(&blockheader.mode, 2);
	swap_bytes(&blockheader.ctcount, 4); // long
	swap_bytes(&blockheader.lpval, 4);   // float
	swap_bytes(&blockheader.rpval, 4);  // float
	swap_bytes(&blockheader.lvl, 4);  // float
	swap_bytes(&blockheader.tlt, 4);  // float
      }
      
      // print out block headers
      if (verbose > 2) {
	mexPrintf("(getfidkmex)  BLOCK (%i,%i)\n",i,j);
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
      // read in the block of data
      if (fread(block, 1, header.tbytes, ffid) == 0) {
	mexPrintf("(getfidkmex) ERROR: reading data from file %s\n",filename);
	errorExit(plhs);
	free(block);
	return;
      }
      // move block into data and datai pointer approriately
      // note that the real and imaginary parts of the data
      // are stored sequentially in the fid file. 
      // these are stored either as a double if originalformat is not set.
      // otherwise, it is stored in its native format into two arrays
      // in originalformat is set.
      for (k=0; k<header.np; k+=2) {
	switch(info.datatype) {
	  // data is stored as a float
	case mxSINGLE_CLASS:
	  // store as floats
	  if (originalformat) {
	    if (swapFlag) {
	      tokenfloat = (((FLOAT*)block)[k]); 
	      swap_bytes(&tokenfloat,4); 
	      data_float[j*info.imagesize+i*info.linelen+k/2] = tokenfloat;
	      tokenfloat = (((FLOAT*)block)[k+1]); 
	      swap_bytes(&tokenfloat,4);
	      datai_float[j*info.imagesize+i*info.linelen+k/2] = tokenfloat;
	    }
	    else {
	      data_float[j*info.imagesize+i*info.linelen+k/2] = ((FLOAT*)block)[k];
	      datai_float[j*info.imagesize+i*info.linelen+k/2] = ((FLOAT*)block)[k+1];
	    }
	  }
	  // or store as doubles (default)
	  else {
	    if (swapFlag) {
	      tokenfloat = (((FLOAT*)block)[k]); 
	      swap_bytes(&tokenfloat,4); 
	      data[j*info.imagesize+i*info.linelen+k/2] = (double)tokenfloat;
	      tokenfloat = (((FLOAT*)block)[k+1]); 
	      swap_bytes(&tokenfloat,4);
	      datai[j*info.imagesize+i*info.linelen+k/2] = (double)tokenfloat;
	    }
	    else {
	      data[j*info.imagesize+i*info.linelen+k/2] = (double)((FLOAT*)block)[k];
	      datai[j*info.imagesize+i*info.linelen+k/2] = (double)((FLOAT*)block)[k+1];
	    }
	  }
	  break;
	  // data is stored as two byte integers
	case mxINT16_CLASS:
	  // store as two byte integers
	  if (originalformat) {
	    if (swapFlag) {
	      tokenint16 = (((INT16*)block)[k]); 
	      swap_bytes(&tokenint16,2); 
	      data_int16[j*info.imagesize+i*info.linelen+k/2] = tokenint16;
	      tokenint16 = (((INT16*)block)[k+1]); 
	      swap_bytes(&tokenint16,2);
	      datai_int16[j*info.imagesize+i*info.linelen+k/2] = tokenint16;
	    }
	    else {
	      data_int16[j*info.imagesize+i*info.linelen+k/2] = ((INT16*)block)[k];
	      datai_int16[j*info.imagesize+i*info.linelen+k/2] = ((INT16*)block)[k+1];
	    }
	  }
	  // or store as doubles (default)
	  else {
	    if (swapFlag) {
	      tokenint16 = (((INT16*)block)[k]); 
	      swap_bytes(&tokenint16,2); 
	      data[j*info.imagesize+i*info.linelen+k/2] = (double)tokenint16;
	      tokenint16 = (((INT16*)block)[k+1]); 
	      swap_bytes(&tokenint16,2);
	      datai[j*info.imagesize+i*info.linelen+k/2] = (double)tokenint16;
	    }
	    else {
	      data[j*info.imagesize+i*info.linelen+k/2] = (double)((INT16*)block)[k];
	      datai[j*info.imagesize+i*info.linelen+k/2] = (double)((INT16*)block)[k+1];
	    }
	  }
	  break;
	  // data is stored as four byte integers
	case mxINT32_CLASS:
	  // save as four byte integer
	  if (originalformat) {
	    if (swapFlag) {
	      tokenint32 = (((INT32*)block)[k]); 
	      swap_bytes(&tokenint32,4); 
	      data_int32[j*info.imagesize+i*info.linelen+k/2] = tokenint32;
	      tokenint32 = (((INT32*)block)[k+1]); 
	      swap_bytes(&tokenint32,4);
	      datai_int32[j*info.imagesize+i*info.linelen+k/2] = tokenint32;
	    }
	    else {
	      data_int32[j*info.imagesize+i*info.linelen+k/2] = ((INT32*)block)[k];
	      datai_int32[j*info.imagesize+i*info.linelen+k/2] = ((INT32*)block)[k+1];
	    }
	  }
	  // or store as doubles (default)
	  else {
	    if (swapFlag) {
	      tokenint32 = (((INT32*)block)[k]); 
	      swap_bytes(&tokenint32,4);
	      data[j*info.imagesize+i*info.linelen+k/2] = (double)tokenint32;
	      tokenint32 = (((INT32*)block)[k+1]); 
	      swap_bytes(&tokenint32,4);
	      datai[j*info.imagesize+i*info.linelen+k/2] = (double)tokenint32;
	    }
	    else {
	      data[j*info.imagesize+i*info.linelen+k/2] = (double)(((INT32*)block)[k]);
	      datai[j*info.imagesize+i*info.linelen+k/2] = (double)(((INT32*)block)[k+1]);
	    }
	  }
	  break;
	}
      }
    }
  }
  // free memory for block
  free(block);
  // close file handles
  fclose(ffid);
  fclose(fprocpar);
}

////////////////////////
// function usageError
////////////////////////
void usageError()
{

  mexPrintf("    program: getfidk.c\n");
  mexPrintf("         by: justin gardner\n");
  mexPrintf("    purpose: reads fid data from VNMR system\n");
  mexPrintf("       date: 05/06/03\n");
  mexPrintf("    compile: mex getfidk.c\n");
  mexPrintf("    history: 2008/12/01, Pei, add endian detect and byte swap part.\n");
  mexPrintf("             Default format of fid is big endian, as we acquaire data\n");
  mexPrintf("             on VNMR Sun Solaris, also epirri5 keep the same endian.\n");
  mexPrintf("             Note: SWAP data firstly and then DO typecast.\n");
  mexPrintf("      USAGE: d = getfidk('filedir',verbose,uint16);\n");
  mexPrintf("                 filename = name of fid directory\n");
  mexPrintf("                 verbose  = (0,1,23) for reporting verbose information\n");
  mexPrintf("                            defaults to no reporting (0)\n");
  mexPrintf("                 uint16   = (1 or 0) to return data in native format or double\n");
  mexPrintf("                            defaults to double (0)\n");
  mexPrintf("                 Retruns data as x,y,slices,receivers,volumes\n\n");
  mexPrintf("                 Note that for 2D data with multiple receivers, because of Varian\n");
  mexPrintf("                 silliness, the receivers dimension and the slices dimensions will be swapped\n");
  mexPrintf("                 getfid fixes this\n");
}

///////////////////////
// function errorexit
///////////////////////
void errorExit(mxArray *plhs[])
{
  if (ffid) fclose(ffid);
  if (fprocpar) fclose(fprocpar);
  mxSetField(plhs[0],0,"data",NULL);

}

//////////////////////////
// function fileNotFound
//////////////////////////
void fileNotFound(char *filename)
{
  mexPrintf("(getfidkmex) ERROR: Could not open %s\n", filename);
}

////////////////////////
// function fileLength
////////////////////////
int fileLength(char *filename,FILE *filepointer)
{
  fpos_t filelen;
  // get size of file: seek to end of file, read pos, seek to begin of file
  if (fseek(filepointer,0,SEEK_END)) {
    mexPrintf("(getfidkmex) ERROR: Could not seek in file %s",filename);
    return -1;
  }

  // find the position, i.e. end of file
  // this may not work with all compilers since fpos_t
  // is not required to be the position in bytes.
  if (fgetpos(filepointer,&filelen)) {
    mexPrintf("(getfidkmex) ERROR: Could not seek in file %s",filename);
    return -1;
  }
  
  // seek back to the beginning of the file
  if (fseek(filepointer,0,SEEK_SET)) {
    mexPrintf("(getfidkmex) ERROR: Could not seek in file %s",filename);
    return -1;
  }
  return (int)filelen;
}  


//////////////////////////////////////////
// function swap bytes
//////////////////////////////////////////
void swap_bytes(void * p, size_t s)
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


