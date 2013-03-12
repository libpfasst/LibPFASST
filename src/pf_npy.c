/* Dump routines for the NumPy binary format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define BUFLEN 256

void pf_dnpy_mkdir(char *dname, int dlen)
{
  dname[dlen] = 0;
  mkdir(dname, 0755);
}

void pf_dnpy_npy(char *fname, int flen, char endian[4], double *arr, int nvars)
{
  char errmsg[BUFLEN];
  char header[256*256];
  unsigned short i, len, pad;
  FILE *fp;

  fname[flen] = 0;
  
  fp = fopen(fname, "wb");
  if (fp == NULL) {
    snprintf(errmsg, BUFLEN, "WARNING: Unable to create npy file (%s)", fname);
    perror(errmsg);
    return;
  }

  /* build numpy header */
  snprintf(header, 256*256, 
           "{'descr': '%s', 'fortran_order': True, 'shape': (%d,), }", 
           endian, nvars);
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

  len = strlen(header);
  
  /* write npy header, v1.0 */
  fwrite("\x93NUMPY\x01\x00", 1, 8, fp);
  fwrite(&len, 1, 2, fp);
  fwrite(header, 1, len, fp);

  /* write data and close */
  fwrite(arr, sizeof(double), nvars, fp);
  fclose(fp);
}

void pf_dnpy_solution_npy(char *dirname, int dlen, char endian[4],
                          double *q, int nvars, 
                          int level, int step, int cycle, int iter)
{
  char fname[BUFLEN];

  dirname[dlen] = 0;
  snprintf(fname, BUFLEN, 
           "%s/level=%d_step=%d_iter=%d_cycle=%d.npy", 
           dirname, level, step, iter, cycle);

  pf_dnpy_npy(fname, strlen(fname), endian, q, nvars);
}
