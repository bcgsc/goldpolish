import base64
import uuid
import threading
import os
import signal
import time

INIT_PID = 1
PARENT_QUERY_PERIOD = 1  # In seconds


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
