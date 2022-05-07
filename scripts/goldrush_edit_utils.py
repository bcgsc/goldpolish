import base64
import uuid
import threading
import os
from os.path import join
import signal
import time
import fcntl

INIT_PID = 1
PARENT_QUERY_PERIOD = 1  # In seconds
THREADS_IN_USE_FILE = "threads_in_use"


def get_random_name():
    """Return a unique string of characters."""
    return (
        base64.urlsafe_b64encode(uuid.uuid4().bytes)
        .rstrip(b"=")
        .replace(b"-", b"")
        .decode("ascii")
    )


def process_terminate():
    """SIGTERM the calling process."""
    os.kill(os.getpid(), signal.SIGTERM)


def watch_process(process):
    """SIGTERM the calling process if the observed process dies with an error."""

    def _watch_process(process):
        process.wait()
        if process.returncode != 0:
            print(f"{process.args} failed!")
            process_terminate()

    threading.Thread(target=_watch_process, args=(process,), daemon=True).start()


def wait_till_parent_ends():
    """Periodically poll for parent's end."""
    while os.getppid() != INIT_PID:
        time.sleep(PARENT_QUERY_PERIOD)


def bind_to_parent():
    """SIGTERM the calling process if the parent process ends."""

    def _bind_to_parent():
        wait_till_parent_ends()
        process_terminate()

    threading.Thread(target=_bind_to_parent, daemon=True).start()


def create_threads_in_use(path, prefix):
    fpath = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    with open(fpath, "w") as f:
        print(0, file=f)


def get_threads_in_use(path, prefix):
    fpath = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    with open(fpath, "r") as f:
        fcntl.flock(f, fcntl.LOCK_SH)
        threads = int(f.read())
        fcntl.flock(f, fcntl.LOCK_UN)
    return threads


def update_threads_in_use(path, prefix, threads_delta):
    fpath = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    with open(fpath, "r+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        threads = int(f.read())
        os.ftruncate(f.fileno(), 0)
        f.seek(0)
        print(threads + threads_delta, file=f)
        os.fsync(f)
        fcntl.flock(f, fcntl.LOCK_UN)
