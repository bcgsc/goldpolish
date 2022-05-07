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
THREADS_IN_USE_LOCK_FILE = "threads_in_use.lock"


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
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    lock_path = join(path, f"{prefix}-{THREADS_IN_USE_LOCK_FILE}")

    with open(threads_in_use_path, "w") as f:
        print(0, file=f)
    with open(lock_path, "w") as f:
        pass


def get_threads_in_use(path, prefix):
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")

    lock_path = join(path, f"{prefix}-{THREADS_IN_USE_LOCK_FILE}")
    lock = os.open(lock_path, os.O_RDONLY)

    fcntl.flock(lock, fcntl.LOCK_SH)
    with open(threads_in_use_path, "r") as threads_in_use:
        threads = int(threads_in_use.read())
    fcntl.flock(lock, fcntl.LOCK_UN)

    os.close(lock)

    return threads


def update_threads_in_use(path, prefix, threads_delta):
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")

    lock_path = join(path, f"{prefix}-{THREADS_IN_USE_LOCK_FILE}")
    lock = os.open(lock_path, os.O_RDONLY)

    fcntl.flock(lock, fcntl.LOCK_EX)
    with open(threads_in_use_path, "r") as threads_in_use:
        threads = int(threads_in_use.read())
    with open(threads_in_use_path, "w") as threads_in_use:
        print(threads + threads_delta, file=threads_in_use)
    fcntl.flock(lock, fcntl.LOCK_UN)

    os.close(lock)
