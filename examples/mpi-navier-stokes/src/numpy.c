/* Dump routines for the NumPy binary format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define BUFLEN 256

void dump_mkdir(char *dname, int dlen)
{
  dname[dlen] = 0;
  mkdir(dname, 0755);
}


void dump_numpy(char *dname, char *fname, char endian[4], int dim, int *shape, int nvars, double *array)
{
  unsigned short  i, len, pad;
  FILE*           fp;
  char            errmsg[BUFLEN];
  char            header[256*256];
  char            buf[BUFLEN], shp[BUFLEN];

  /*
   * open output file
   */

  snprintf(buf, BUFLEN, "%s/%s", dname, fname);
  fp = fopen(buf, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create npy file (%s)", buf);
    perror(errmsg);
    return;
  }

  /*
   * build shape string
   */
  shp[0] = 0;
  for (i=0; i<dim; i++) {
    snprintf(shp + strlen(shp), BUFLEN-strlen(shp)-1,
             "%d,", shape[i]);
  }

  /* 
   * build numpy header 
   */
  snprintf(header, 256*256, 
           "{'descr': '%s', 'fortran_order': True, 'shape': (%s), }", 
           endian, shp);

  len = strlen(header);
  if (len > 256*254) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create npy file (%s)", fname);
    fprintf(stderr, "%s: NumPy header too long.\n", errmsg);
    return;
  }

  pad = 16 - (8 + len + 1) % 16;
  for (i=len; i<len+pad; i++)
    header[i] = ' ';
  header[len+pad+1] = '\n';
  header[len+pad+2] = '\0';

  /* 
   * write npy header, v1.0 
   */
  len = strlen(header);
  
  fwrite("\x93NUMPY\x01\x00", 1, 8, fp);
  fwrite(&len, 1, 2, fp);
  fwrite(header, 1, len, fp);

  /*
   * write data and close 
   */
  fwrite(array, sizeof(double), nvars, fp);
  fclose(fp);
}
