/*
!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!
*/

/* C routines to deal with pthreads mutexes, conditions etc 
 *
 * The PTH version of PFASST creates instances of the pf_pth_t struct
 * for each processor and each level.  This struct contains the
 * mutexes and conditions used to syncronize PFASST communications.
 *
 * To send (per level):
 * 
 *   1) The sending processor waits until its "received" condition is
 *      cleared (set to 0).  This signals that the next processor has
 *      received the previous message sent and that the send buffer
 *      can be overwritten.
 *
 *   2) The sending processor locks its send buffer and copies qend to
 *      it.
 *
 *   3) The sending processor sets its send and recv conditions to the
 *      given tag.
 *
 * To receive (per level):
 *
 *   1) The receiving processor looks up the sending processors
 *      pf_pth_t struct.
 * 
 *   2) The receiving processor waits until the senders "send"
 *      condition matches the given tag.
 *
 *   3) The receiving processor locks the sending processors send
 *      buffer and copies it to q0.
 *
 *   4) The receiving processor clears the senders "received"
 *      condition to signal that it has received the message and the
 *      sending processor is free to use its send buffer again.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/*
 * PFASSTH Pthread "pth" struct.
 */
typedef struct {
  
  int send;
  int recv;

  pthread_cond_t  send_cond;
  pthread_mutex_t send_mtx;

  pthread_cond_t  recv_cond;
  pthread_mutex_t recv_mtx;
  
} pf_pth_t;


pf_pth_t *pf_pth_create(void) {
  pf_pth_t *pth;

  pth = (pf_pth_t *) malloc(sizeof(pf_pth_t));

  pthread_cond_init(&pth->send_cond, NULL);
  pthread_mutex_init(&pth->send_mtx, NULL);

  pthread_cond_init(&pth->recv_cond, NULL);
  pthread_mutex_init(&pth->recv_mtx, NULL);

  pth->send = -1;
  pth->recv = 0;

  return pth;
}

void pf_pth_destroy(pf_pth_t *pth) {
  pthread_cond_destroy(&pth->send_cond);
  pthread_mutex_destroy(&pth->send_mtx);

  pthread_cond_destroy(&pth->recv_cond);
  pthread_mutex_destroy(&pth->recv_mtx);

  free(pth);
}

void pf_pth_wait_send(pf_pth_t *pth, int tag) {
  pthread_mutex_lock(&pth->send_mtx);
  while (pth->send != tag) {
    pthread_cond_wait(&pth->send_cond, &pth->send_mtx);
  }
  pthread_mutex_unlock(&pth->send_mtx);
}

void pf_pth_wait_recv(pf_pth_t *pth, int tag) {
  pthread_mutex_lock(&pth->recv_mtx);
  while (pth->recv != tag) {
    pthread_cond_wait(&pth->recv_cond, &pth->recv_mtx);
  }
  pthread_mutex_unlock(&pth->recv_mtx);
}

void pf_pth_set_send(pf_pth_t *pth, int tag) {
  pthread_mutex_lock(&pth->send_mtx);
  pth->send = tag;
  pthread_cond_broadcast(&pth->send_cond);
  pthread_mutex_unlock(&pth->send_mtx);
}

void pf_pth_set_recv(pf_pth_t *pth, int tag) {
  pthread_mutex_lock(&pth->recv_mtx);
  pth->recv = tag;
  pthread_cond_broadcast(&pth->recv_cond);
  pthread_mutex_unlock(&pth->recv_mtx);
}

void pf_pth_lock(pf_pth_t *pth) {
  pthread_mutex_lock(&pth->send_mtx);
}

void pf_pth_unlock(pf_pth_t *pth) {
  pthread_mutex_unlock(&pth->send_mtx);
}
