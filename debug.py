
import zmq
import ctypes
import urwid
import threading
import time
import struct

from ctypes import c_double, c_int

maxprocs = 16

class pf_state_t(ctypes.Structure):
    _fields_ = [ ('t0', c_double),
                 ('dt', c_double),
                 ('nsteps', c_int),
                 ('block', c_int),
                 ('cycle', c_int),
                 ('step', c_int),
                 ('iter', c_int),
                 ('level', c_int),
                 ('hook', c_int),
                 ('proc', c_int),
                 ('status', c_int),
                 ('pstatus', c_int),
                 ('nmoved', c_int),
                 ('first', c_int),
                 ('last', c_int),
                 ('res', c_double) ]

    def from_bytes(self, bytes):
        fit = min(len(bytes), ctypes.sizeof(self))
        ctypes.memmove(ctypes.addressof(self), bytes, fit)


class Acquire(threading.Thread):
    """Listen for data and call 'render_scene' when appropriate.

    This thread creates a ZMQ socket (bound to port 31415) and listens
    for data.  When data is received, the GUI thread's 'render_scene'
    method is called.

    A non-blocking ZMQ receive is used so that we can periodically
    check our 'stop' event and exit if appropriate.  The 'stop' event
    is set by the main thread when the user closes the main window.

    """

    def __init__(self):
        super(Acquire, self).__init__()
        self.stop = threading.Event()
        self.viewer = None
        

    def recv(self, flag=0):
        try:
            msg = self.socket.recv(flag)
        except zmq.ZMQError:
            return None, None

        size, wlen = struct.unpack_from("i i", msg)
        size, wlen, sbuf, where = struct.unpack("i i %ds %ds" % (size, wlen), msg)

        state = pf_state_t()
        state.from_bytes(sbuf)

        return state, where


    def run(self):

        self.context = zmq.Context()
        self.socket  = self.context.socket(zmq.PULL)
        self.socket.bind("tcp://*:31415")

        while True:
            if self.stop.is_set():
                return

            state, where = self.recv(zmq.NOBLOCK)

            if state is None:
                time.sleep(0.005)
                continue

            if self.viewer is not None:
                self.viewer.render(state, where)

log = []

class MainWindow():
    """Build main window and render scene when called.

    When the user closes the main window the Acquire thread's stop
    event is set.
    """

    tmpl  = "{proc:2d} {step:3d} {where:>16s} {first:2d} {last:2d} {nmoved:2d} {status:>2s} {pstatus:>2s} {res:e}"
    title = " P   S                W  F  L NM  S PS"

    def render(self, state, where):
        stat = { 1: 'I', 2: 'C', 3: 'P' }

        s = self.tmpl.format(proc=state.proc, step=state.step, where=where,
                             first=state.first+1, last=state.last+1, nmoved=state.nmoved,
                             status=stat[state.status], pstatus=stat[state.pstatus], res=state.res)
        log.append(s)
        self.pile.contents[state.proc] = (urwid.Text(s), self.pile.options())
        self.loop.draw_screen()

    def exit(self, key):
        if key in ('q', 'Q'):
            self.acquire.stop.set()
            self.acquire.join()
            raise urwid.ExitMainLoop()
        

    def __init__(self):
        self.acquire = Acquire()
        self.acquire.viewer = self
        self.acquire.start()

        pile = urwid.Pile([ urwid.Text(self.title) ] + [ urwid.Text('-') for i in range(maxprocs) ])
        # pile = urwid.Pile([ urwid.Text('-') for i in range(maxprocs+1) ])
        loop = urwid.MainLoop(urwid.Filler(pile), unhandled_input=self.exit)
        
        self.pile = pile
        self.loop = loop


if __name__ == '__main__':
    main = MainWindow()
    main.loop.run()
    print '\n'.join(log)


    

