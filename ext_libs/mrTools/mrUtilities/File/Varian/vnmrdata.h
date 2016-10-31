/* THIS FILE IS TAKEN FROM vnmr/psg/data.h */
/**/
/**/
/* static char SCCSid[] = "@(#)data.h 14.2 04/25/00 Copyright (c) 1991-2000 Varian, Inc. All Rights Reserved"; */
/*
* Varian, Inc. All Rights Reserved.
* This software contains proprietary and confidential
* information of Varian, Inc. and its contributors.
* Use, disclosure and reproduction is prohibited without
* prior consent.
*/
/*---------------------------------------
|					|
|   data.h   -   data file handler 	|
|					|
+--------------------------------------*/

/*-------------------------------
|				|
|  Rene Richarz 12-85		|
|  Sandy Farmer  4-89		|
|  Sandy Farmer 10-89		|
|				|
+------------------------------*/

/*---------------------------------------------------------------
0xc0								|
|  The data file handler allows processing of NMR data files.   |
|  The 1D NMR data files have the following structure:          |
|								|
|  filehead  blockhead  blockdata  blockhead  blockdata ...  	|
|								|
|								|
|  The 2D NMR data files have the following structure:		|
|								|
|  filehead  blockhead  blockhead2  blockdata ...		|
|								|
|								|
|  The 3D NMR data files have the following structure:		|
|								|
|  filehead  blockhead  blockhead2  blockhead2  blockdata ..	|
|								|
|								|
|  All blocks within one file have exactly the same size.       |
|  Two files can be processed at the same time.  The first is	|
|  always the file "data" of the current experiment. It is	|
|  left open and can be used by any module.  The second one	|
|  can be used for "acqfil/fid" or other files and is only	|
|  used within one module.					|
|								|
+--------------------------------------------------------------*/

#define VERSION_1

#define MAX_DATA_BUFFERS   16	/* maximum nubmer of data buffers per exp	 */
#define D_MAXFILES	3	/* maximum number of open files		 */
#define D_DATAFILE	0	/* curexp/data file index 		 */
#define D_PHASFILE	1	/* curexp/phasefile file index		 */
#define D_USERFILE	2	/* any user file used in one module	 */

/*  The file headers are defines as follows:                     */

#ifdef __x86_64__
#define int32 int
#define int16 short
#else
#define int32 long
#define int16 short
#endif

/*****************/
struct datafilehead
/*****************/
/* Used at the beginning of each data file (fid's, spectra, 2D)  */
{
   int32    nblocks;       /* number of blocks in file			*/
   int32    ntraces;       /* number of traces per block			*/
   int32    np;            /* number of elements per trace		*/
   int32    ebytes;        /* number of bytes per element		*/
   int32    tbytes;        /* number of bytes per trace			*/
   int32    bbytes;        /* number of bytes per block			*/
   int16   vers_id;       /* software version and file_id status bits	*/
   int16   status;        /* status of whole file			*/
   int32	   nbheaders;	  /* number of block headers			*/
};

/*******************/
struct datablockhead
/*******************/
/* Each file block contains the following header        */
{
   int16 	  scale;	/* scaling factor                   */
   int16 	  status;	/* status of data in block          */
   int16 	  index;	/* block index                      */
   int16	  mode;		/* mode of data in block	    */
   int32		  ctcount;	/* ct value for FID		    */
   float 	  lpval;	/* F2 left phase in phasefile       */
   float 	  rpval;	/* F2 right phase in phasefile      */
   float 	  lvl;		/* F2 level drift correction        */
   float 	  tlt;		/* F2 tilt drift correction         */
};

/*********************/
struct hypercmplxbhead
/*********************/
/* Additional datafile block header for hypercomplex 2D data */
{
   int16	s_spare1;	/* int16 word:  spare		*/
   int16	status;		/* status word for block header	*/
   int16	s_spare2;	/* int16 word:  spare		*/
   int16	s_spare3;	/* int16 word:  spare		*/
   int32	l_spare1;	/* int32 word:  spare		*/
   float 	lpval1;		/* additional left phase	*/
   float 	rpval1;		/* additional right phase 	*/
   float 	f_spare1;	/* float word:  spare		*/
   float 	f_spare2;	/* float word:  spare		*/
};

/*****************/
struct datapointers
/*****************/
{
   struct datablockhead	*head;	/* pointer to the head of block     */
   float		*data;	/* pointer to the data of the block */
};


typedef struct datafilehead	dfilehead;
typedef struct datablockhead	dblockhead;
typedef struct datapointers	dpointers;
typedef struct hypercmplxbhead	hycmplxhead;


/*---------------------------------------------------------------
|								|
| The data file handler provides a set of subroutines, which    |
| allow to buffer several blocks of a file in memory at a time. |
| The number of blocks buffered in memory is limited by the     |
| factor (D_MAXBYTES*bufferscale) in the data file handler.     |
| All data file handler routines return an error code.		|
|								|
+--------------------------------------------------------------*/

#define D_IS_OPEN  -101  /* there is a open file currently       */
#define D_NOTOPEN  -102  /* there is no open file currently      */
#define D_READERR  -103  /* there has been a file read error     */
#define D_WRITERR  -104  /* there has been a file write error    */
#define D_INCONSI  -105  /* inconsistent data head               */
#define D_NOALLOC  -106  /* not allocated, out of memory         */
#define D_NOBLOCK  -107  /* too many blocks in memory            */
#define D_NOTFOUND -108  /* entry with index not found           */
#define D_SEEKERR  -109  /* seek in file not successful          */
#define D_BLINDEX  -110  /* illegal block index                  */
#define D_FLINDEX  -111  /* illegal file index			 */
#define D_BUFERR   -112  /* internal buffer too small		 */
#define D_TRERR    -113  /* cannot transpose block               */
#define D_READONLY -114	 /* no write permission			 */
#define D_INVNBHDR -115  /* invalid number of block headers	 */
#define D_INVOP    -116	 /* invalid file operation		 */
#define D_TRUNCERR -117  /* a file truncation error has occurred */
#define D_OPENERR  -118  /* error opening a data file            */


/*-----------------------------------------------
|						|
|     Data file header status bits  (0-9) 	|
|						|
+----------------------------------------------*/

/* Bits 0-6:  file header and block header status bits (bit 6 = unused) */
#define S_DATA		0x1	/* 0 = no data      1 = data       */
#define S_SPEC		0x2	/* 0 = fid	    1 = spectrum   */
#define S_32		0x4	/* 0 = 16 bit	    1 = 32 bit	   */
#define S_FLOAT 	0x8	/* 0 = integer      1 = fltng pt   */
#define S_COMPLEX	0x10	/* 0 = real	    1 = complex	   */
#define S_HYPERCOMPLEX	0x20	/* 1 = hypercomplex		   */

/* Bits 7-10:  file header status bits (bit 10 = unused) */
#define S_ACQPAR	0x80    /* 0 = not acqpar   1 = acqpar     */
#define S_SECND		0x100	/* 0 = first ft     1 = second ft  */
#define S_TRANSF	0x200	/* 0 = regular      1 = transposed */
#define S_3D		0x400	/* 1 = 3D data			   */

/* Bits 11-14:  status bits for processed dimensions (file header) */
#define S_NP		0x800	/* 1 = np  dimension is active	*/
#define S_NF		0x1000	/* 1 = nf  dimension is active	*/
#define S_NI		0x2000	/* 1 = ni  dimension is active	*/
#define S_NI2		0x4000	/* 1 = ni2 dimension is active	*/


/*-------------------------------------------------------
|							|
|     Main data block header status bits  (7-15)	|
|							|
+------------------------------------------------------*/

/* Bit 7 */
#define MORE_BLOCKS	0x80	/* 0 = absent	    1 = present	   */

/* Bits 8-11:  bits 12-14 are currently unused */
#define NP_CMPLX	0x100	/* 1 = complex     0 = real	*/
#define NF_CMPLX	0x200	/* 1 = complex     0 = real	*/
#define NI_CMPLX	0x400	/* 1 = complex     0 = real	*/
#define NI2_CMPLX	0x800   /* 1 = complex     0 = real	*/


/*-----------------------------------------------
|						|
|  Main data block header mode bits  (0-13)	|
|						|
+----------------------------------------------*/

/* Bits 0-3:  bit 3 is used for pamode*/
#define NP_PHMODE	0x1	/* 1 = ph mode  */
#define NP_AVMODE	0x2	/* 1 = av mode	*/
#define NP_PWRMODE	0x4	/* 1 = pwr mode */
#define NP_PAMODE	0x8	/* 1 = pa mode */

/* Bits 4-7:  bit 7 is used for pamode */
#define NF_PHMODE	0x10	/* 1 = ph mode	*/
#define NF_AVMODE	0x20	/* 1 = av mode	*/
#define NF_PWRMODE	0x40	/* 1 = pwr mode	*/
#define NF_PAMODE	0x80	/* 1 = pa mode */

/* Bits 8-10 */
#define NI_PHMODE	0x100	/* 1 = ph mode	*/
#define NI_AVMODE	0x200	/* 1 = av mode	*/
#define NI_PWRMODE	0x400	/* 1 = pwr mode	*/

/* Bits 11-13:  bits 14 and 15 are used for pamode */
#define NI2_PHMODE	0x800	/* 1 = ph mode	*/
#define NI2_AVMODE	0x1000	/* 1 = av mode	*/
#define NI2_PWRMODE	0x2000	/* 1 = pwr mode	*/

/* Bits 14-15 */
#define NI_PAMODE	0x4000	/* 1 = pa mode  */
#define NI2_PAMODE	0x8000	/* 1 = pa mode  */

#define NP_DSPLY	(NP_PHMODE  |  NP_AVMODE  |  NP_PWRMODE  |  NP_PAMODE)
#define NI_DSPLY	(NI_PHMODE  |  NI_AVMODE  |  NI_PWRMODE  |  NI_PAMODE)
#define NI2_DSPLY	(NI2_PHMODE |  NI2_AVMODE |  NI2_PWRMODE |  NI2_PAMODE)


/*-----------------------------------------------
|						|
|    Usage Bits for Additional Block Headers	|
|						|
+----------------------------------------------*/

#define U_HYPERCOMPLEX	0x2	/* 1 = hypercomplex block structure */


/*-------------------------------
|				|
|	File Header Types	|
|				|
+------------------------------*/

#define FILEHEAD	0
#define BLOCKHEAD	1


/*-------------------------------
|				|
|     Processing Dimensions	|
|				|
+------------------------------*/

#define F2_DIM		0	/* first dimension:   np       */
#define F1_DIM		1	/* second dimension:  ni or nf */
#define F3_DIM		2	/* second dimension:  ni2      */
#define F1F2_DIM	101
#define F3F2_DIM	102


/*-------------------------------
|				|
|       Data Definitions        |
|				|
+------------------------------*/

#define REAL		1
#define COMPLEX		2
#define HYPERCOMPLEX	4


/*-------------------------------
|				|
|	 I/O Definitions	|
|				|
+------------------------------*/

#define READ_DATA	0
#define WRITE_DATA	1


/*-------------------------------
|				|
|    Software Version Status 	|
|				|
+------------------------------*/

/* Bits 0-5:  31 different software versions; 32-63 are for 3D */
#ifdef VERSION_1
#define VERSION         1
#endif

#ifdef VERSION_2
#define VERSION         2
#endif

#define VNMR_VERSION	VERSION		/* for Ft3d */


/*-------------------------------
|				|
|	 File ID Status		|
|				|
+------------------------------*/

/* Bits 6-10:  31 different file types (64-1984 in steps of 64) */
#define FID_FILE	0x40
#define DATA_FILE	0x80
#define PHAS_FILE	0xc0
#define DATA3D_FILE	0x100		/* for Ft3d */


/*-------------------------------
|				|
|	 Vendor ID Status	|
|				|
+------------------------------*/

/* Bits 11-14 */
#define S_CM		0x800   /* 1 = Chem. Magnetics data   	*/ 
#define S_GE		0x1000  /* 1 = GE data			*/ 
#define S_JEOL		0x2000  /* 1 = JEOL data		*/ 
#define S_BRU		0x4000  /* 1 = Bruker data		*/ 

#define P_VERS		0x003F	/* preserves version number	*/
#define P_FILE_ID	0x07C0	/* preserves file ID status	*/
#define P_VENDOR_ID	0x7800	/* preserves vendor ID status	*/


/*-------------------------------------------------------
|							|
|  Some data file handler routines return pointers to	|
|  the data in the following structure.			|
|							|
+------------------------------------------------------*/

#if !defined( FT3D ) && !defined( PLEXTRACT )
/******************/
extern int D_init();
/******************/
/*---------------------------------------
|   Initialize the data file handler;	|
|   should be called at bootup		|
+--------------------------------------*/


/*******************************************/
extern int D_open(/* fileindex,file path */);
/*******************************************/
/*-------------------------------------------------------
|    Open a file for processing using the data file	|
|    handler:						|
|							|
|    int fileindex;    defines which file to use	|
|    char *filepath;   string defining the file path 	|
|		       of file				|
+------------------------------------------------------*/


/**********************************/
extern int D_close(/* fileindex */);
/**********************************/
/*-------------------------------------------------------
|    Close the file currently processed by the data	|
|    file handler; store all updated data currently	|
|    in the buffer:					|
|							|
|    int fileindex;      defines which file to use	|
+------------------------------------------------------*/


/**********************/
extern int D_compress();
/**********************/
/*-------------------------------------------------------
|    Compress the data file to real data		|
+------------------------------------------------------*/


/**********************************/
extern int D_flush(/* fileindex */);
/**********************************/
/*-------------------------------------------------------
|   Remove all released blocks from the current buffer	|
|   and store the updated and released blocks on disk:	|
|							|
|   int fileindex;      defines which file to use	|
+------------------------------------------------------*/


/*******************************************************/
extern int D_getfilepath(/* fileindex, path, basepath*/);
/*******************************************************/
/*-------------------------------------------------------
|   Get the file path for the particular data file:	|
|							|
|   char *path;		file path			|
|   char *basepath;	base for the file path		|
|   int  fileindex;	defines which file		|
+------------------------------------------------------*/


/********************8***********************************/
extern int D_getparfilepath(/* partree, path, basepath*/);
/***********************8********************************/ 
/*------------------------------------------------------- 
|   Get the file path for the particular data file:     | 
|                                                       |
|   char *path;         file path                       | 
|   char *basepath;     base for the file path          |
|   int  partree;       defines which parameter tree    |
+------------------------------------------------------*/


/***********************************************/
extern int D_writeintoallblocks(/* fileindex */);
/***********************************************/
/*-------------------------------------------------------
|   Write block header into all data blocks on disk	|
|							|
|   int fileindex;	defines which file to use	|
+------------------------------------------------------*/


/***************************************************/
extern int D_allocbuf(/* fileindex,index,&bpntrs */);
/***************************************************/
/*-------------------------------------------------------
|    Allocate memory for the block index of the file. 	|
|    If this block is subsequently marked as updated,	|
|    all data in the disk file of this block will be	|
|    overwritten.  It returns a pointer to the block	|
|    head and the data of this block:			|
|							|
|    int fileindex;                 defines which file	|
|				    to use		|
|    int index;                     block index		|
|    struct datapointers *bpntrs;   block pointers	|
+------------------------------------------------------*/


/**********************************************************/
extern int D_getbuf(/* fileindex,blocks,index,&bpntrs */);
/**********************************************************/
/*-------------------------------------------------------
|    Allocate memory for the indexed block; then read	|
|    this block in from the disk file.  A pointer to	|
|    the block head the the data of that block is	|
|    returned:						|
|							|
|    int fileindex;      	     defines which file	|
|				     to use		|
|    int blocks;		     number of blocks	|
|    int index;		 	     block index	|
|    struct datapointers *bpntrs;    block pointers	|
+------------------------------------------------------*/


/******************************************/
extern int D_release(/* fileindex,index */);
/******************************************/
/*-------------------------------------------------------
|   Release the block with the specified index so that	|
|   memory can be used again to store other blocks:	|
|							|
|   int fileindex;      defines which file to use	|
|   int index;		block index			|
+------------------------------------------------------*/


/**********************************************/
extern int D_markupdated(/* fileindex,index */);
/**********************************************/
/*-----------------------------------------------
|   Mark the indexed block to be updated:	|
|						|
|   int fileindex;   defines which file to use	|
|   int index;	     block index		|
+----------------------------------------------*/


/*******************************************/
extern int D_getblhead(/* fileindex,index,&blhead */);
/*******************************************/
/*-------------------------------------------------------
|    Gives a program a copy of the data file head:      |
|                                                       |
|    int fileindex;                 defines which file  |
|                                   to use              |
|    int index;                     block index		|
|    struct datafilehead *blhead;   pointer to the      |
|                                   block head          |
+------------------------------------------------------*/


/*******************************************/
extern int D_gethead(/* fileindex,&fhead */);
/*******************************************/
/*-------------------------------------------------------
|    Gives a program a copy of the data file head: 	|
|							|
|    int fileindex;		    defines which file	|
|				    to use		|
|    struct datafilehead *fhead;    pointer to the	|
|				    file head		|
+------------------------------------------------------*/


/****************************************************/
extern int D_newhead(/* fileindex,filepath,&fhead */);
/****************************************************/
/*---------------------------------------------------------------
|    Open a disk file with a new file header; the currrent	|
|    contents of the file are lost:				|
|								|
|    int fileindex;               defines which file to use	|
|    char *filepath;              string defining the file path	|
|			          of file			|
|    struct datafilehead *fhead;  pointer to new file head	|
+--------------------------------------------------------------*/


/****************************************************/
extern int D_foldt(/* fileindex,data,seekmode,seekbytes,iomode,iobytes */);
/****************************************************/
/*---------------------------------------------------------------
|    Perform specialized I/O operations associated with the	|
|    FOLDT algorithm:						|
|								|
|    int fileindex;	defines which file to use		|
|    float *data;	pointer to data				|
|    int seekmode;	mode for lseek() function call		|
|    int seekbytes;	number of bytes to seek past		|
|    int iomode;	mode for I/O operation (read or write)	|
|    int iobytes;	number of bytes in I/O operation	|
+--------------------------------------------------------------*/


/**************************/
extern int D_error(/* e */);
/**************************/
/*-----------------------------------------------
|  Display an error message for a data handler	|
|  error:					|
|						|
|  int e;	error code			|
+----------------------------------------------*/


/*************************************/
extern void setfilepaths(/* datano */);
/*************************************/
#endif

#undef VERSION_1
