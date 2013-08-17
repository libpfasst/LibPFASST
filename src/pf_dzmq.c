/* DZMQ related debugging routines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zmq.h>

typedef struct {
  double t0, dt;
  int nsteps, block, cycle, step, iter, level, hook, proc;
  int status, pstatus, nmoved, first, last;
  double res;
} pf_state_t;

typedef struct {
  void *context;
  void *socket;
} zmq_ctx_t;

void dzmq_free(void *data, void *hint)
{
  free(data);
}

void *dzmq_connect()
{
  zmq_ctx_t *ctx;
  int rc;

  ctx = (zmq_ctx_t *) malloc(sizeof(zmq_ctx_t));
  if (ctx == NULL) {
    perror("Error: zmq allocate:");
    return NULL;
  }

  ctx->context = zmq_init(1);
  if (ctx->context == NULL) {
    fprintf(stderr, "Error: zmq create: %s\n",
  	    zmq_strerror(errno));
  }

  ctx->socket = zmq_socket(ctx->context, ZMQ_PUSH);
  rc = zmq_connect(ctx->socket, "tcp://127.0.0.1:31415");

  if (rc != 0) {
    fprintf(stderr, "Error: zmq connect: %s\n",
  	    zmq_strerror(errno));
    // XXX: tear down properly
    return NULL;
  }

  return ctx;
}

void dzmq_close(void *ptr)
{
  zmq_ctx_t *ctx;
  ctx = (zmq_ctx_t *) ptr;

  zmq_close(ctx->socket);
  zmq_term(ctx->context);
}

void dzmq_status(void *ptr, pf_state_t *state, char *where, int wlen)
{
  zmq_ctx_t *ctx;
  zmq_msg_t msg;
  int size;
  void *qp;
  char *buf;

  ctx = (zmq_ctx_t *) ptr;

  size = 2*sizeof(int) + sizeof(pf_state_t) + wlen;
  qp = malloc(size);

  buf = qp;
  *((int*) buf) = (int) sizeof(pf_state_t); buf += sizeof(int);
  *((int*) buf) = wlen;                     buf += sizeof(int);
  memcpy(buf, state, sizeof(pf_state_t));   buf += sizeof(pf_state_t);
  memcpy(buf, where, wlen);                 buf += wlen;

  zmq_msg_init_data(&msg, qp, size, dzmq_free, NULL);
  zmq_send(ctx->socket, &msg, 0);
}
