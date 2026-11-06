"""Execute the shipped notebook in this Python environment, with a ten-minute limit."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import tempfile
import time


def worker(notebook, output):
    import nbformat
    from nbclient import NotebookClient
    from jupyter_client.kernelspec import KernelSpecManager

    output.mkdir(parents=True, exist_ok=True)
    nb = nbformat.read(notebook, as_version=4)
    # Temporary kernel specification selects this interpreter without changing
    # the user's global Jupyter configuration or relying on a named venv.
    with tempfile.TemporaryDirectory(prefix="pase-kernel-") as temp:
        resource = Path(temp) / "kernels" / "pase-release-check"
        resource.mkdir(parents=True)
        (resource / "kernel.json").write_text(json.dumps({
            "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
            "display_name": "PASE verification", "language": "python"
        }), encoding="utf-8")
        manager = KernelSpecManager(kernel_dirs=[str(resource.parent)])
        client = NotebookClient(nb, timeout=120, kernel_name=resource.name,
                                resources={"metadata": {"path": str(output)}},
                                allow_errors=False)
        from jupyter_client import KernelManager
        client.km = KernelManager(kernel_name=resource.name, kernel_spec_manager=manager)
        try:
            client.execute()
        finally:
            nbformat.write(nb, output / notebook.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", type=Path, default=Path(__file__).resolve().parents[1] / "examples" / "gse234933_quickstart.ipynb")
    parser.add_argument("--output", type=Path, default=Path("tutorial-output"))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    notebook, output = args.notebook.resolve(), args.output.resolve()
    if args.worker:
        worker(notebook, output)
        return
    output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    started = datetime.now(timezone.utc).isoformat()
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--notebook", str(notebook), "--output", str(output)]
    timed_out = False
    with (output / "execution.log").open("w", encoding="utf-8") as log:
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=os.name != "nt")
        try:
            exit_code = proc.wait(timeout=600)
        except subprocess.TimeoutExpired:
            timed_out = True
            if os.name == "nt":
                subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)],
                               stdout=log, stderr=subprocess.STDOUT, check=False)
            else:
                os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            exit_code = 124
    elapsed = time.perf_counter() - start
    report = dict(started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=elapsed, limit_seconds=600, timed_out=timed_out,
                  exit_code=exit_code, passed=exit_code == 0 and elapsed <= 600,
                  python=platform.python_version(), platform=platform.platform(),
                  scope="Notebook execution after installation; no raw-source reconstruction")
    (output / "runtime.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(exit_code or 1)


if __name__ == "__main__":
    main()

