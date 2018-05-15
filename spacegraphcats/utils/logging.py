"""Log the history of commands that have been run."""
import os
import subprocess as sp

LOGFILE = "commands.log"
GIT_COMMIT = ["git", "rev-parse", "--short",  "HEAD"]


def log(project, args):
    """Log a run."""
    path = os.path.join(project, LOGFILE)
    output = sp.run(GIT_COMMIT, stdout=sp.PIPE)
    git_hash = output.stdout.strip().decode("utf-8")
    with open(path, 'a') as logfile:
        call = " ".join(args)
        logfile.write("{} {}\n".format(git_hash, call))
